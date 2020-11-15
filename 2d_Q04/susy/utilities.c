// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force

// #define DET_DIST prints out all determinants for plotting distribution
// CAUTION: Do not run DET_DIST with MPI!

//#define DET_DIST
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute at each site all NUMLINK * (NUMLINK - 1) plaquette determinants
// counting both orientations, and saving ZW* = plaqdet (plaqdet - 1)^*
// Use Tr_Uinv as temporary storage
void compute_plaqdet() {
  register int i;
  register site *s;
  char **local_pt[2][2];
  int a, b, flip = 0;
  complex tc;
  msg_tag *tag0[2], *tag1[2];

#ifdef DET_DIST
  if (this_node != 0) {
    printf("compute_plaqdet: don't run DET_DIST in parallel\n");
    fflush(stdout);
    terminate(1);
  }
#endif

  for (a = 0; a < 2; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[2 + a];
  }

  // Gather determinants rather than the full matrices
  // Recall det[Udag] = (det[U])^*
  FORALLSITES(i, s) {
    FORALLDIR(a)
      Tr_Uinv[a][i] = find_det(&(s->link[a]));
  }

  // Start first set of gathers (a = 0 and b = 1)
  // local_pt[0][0] is det[U_1(x+0)], local_pt[0][1] is det[U_0(x+1)]
  tag0[0] = start_gather_field(Tr_Uinv[1], sizeof(complex),
                               goffset[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(Tr_Uinv[0], sizeof(complex),
                               goffset[1], EVENANDODD, local_pt[0][1]);

  // Main loop
  FORALLDIR(a) {
    FORALLDIR(b) {
      if (a == b)
        continue;

      if (a == 0) {         // Start other set of gathers (a = 1 and b = 0)
        tag0[1] = start_gather_field(Tr_Uinv[0], sizeof(complex),
                                     goffset[1], EVENANDODD, local_pt[1][0]);
        tag1[1] = start_gather_field(Tr_Uinv[1], sizeof(complex),
                                     goffset[0], EVENANDODD, local_pt[1][1]);
      }

      // Initialize plaqdet[a][b] with det[U_b(x)] det[Udag_a(x)]
      FORALLSITES(i, s)
        CMULJ_(Tr_Uinv[a][i], Tr_Uinv[b][i], plaqdet[a][b][i]);

      // Now put it all together
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        // local_pt[flip][0] is det[U_b(x+a)]
        // Conjugate it to get det[Udag_b(x+a)]
        CMUL_J(plaqdet[a][b][i], *((complex *)(local_pt[flip][0][i])), tc);
        // local_pt[flip][1] is det[U_a(x+b)]
        CMUL(*((complex *)(local_pt[flip][1][i])), tc, plaqdet[a][b][i]);

        // ZWstar = plaqdet (plaqdet - 1)^*
        CADD(plaqdet[a][b][i], minus1, tc);
        CMUL_J(plaqdet[a][b][i], tc, ZWstar[a][b][i]);
#ifdef DET_DIST
        if (a < b) {
          printf("DET_DIST %d %d %d %d %.4g %.4g %.4g\n",
                 s->x, s->t, a, b,
                 plaqdet[a][b][i].real, plaqdet[a][b][i].imag, cabs_sq(&tc1));
        }
#endif
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip++;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Save U_a^{-1} and Udag_a^{-1} = (U_a^{-1})^dag at each site
void compute_Uinv() {
  register int i, mu;
  register site *s;

  FORALLSITES(i, s) {
    FORALLDIR(mu) {
      invert(&(s->link[mu]), &(Uinv[mu][i]));
      adjoint(&(Uinv[mu][i]), &(Udag_inv[mu][i]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Separate routines for each term in the fermion operator
// All called by fermion_op at the bottom of the file
#ifdef VP
void Dplus(matrix *src[NUMLINK], matrix *dest) {
  register int i;
  register site *s;
  msg_tag *tag[4];

  // Only have one set of gathers (mu = 0 and nu = 1)
  tag[0] = start_gather_field(src[1], sizeof(matrix),
                              goffset[0], EVENANDODD, gen_pt[0]);

  tag[1] = start_gather_site(F_OFFSET(link[0]), sizeof(matrix),
                             goffset[1], EVENANDODD, gen_pt[1]);

  tag[2] = start_gather_field(src[0], sizeof(matrix),
                              goffset[1], EVENANDODD, gen_pt[2]);

  tag[3] = start_gather_site(F_OFFSET(link[1]), sizeof(matrix),
                             goffset[0], EVENANDODD, gen_pt[3]);

  wait_gather(tag[0]);
  wait_gather(tag[1]);
  wait_gather(tag[2]);
  wait_gather(tag[3]);
  FORALLSITES(i, s) {
    // Initialize dest[i]
    scalar_mult_nn(&(s->link[0]), (matrix *)(gen_pt[0][i]),
                   s->bc[0], &(plaq_dest[i]));

    // Add or subtract the other three terms
    mult_nn_dif(&(src[1][i]), (matrix *)(gen_pt[1][i]), &(plaq_dest[i]));
    scalar_mult_nn_dif(&(s->link[1]), (matrix *)(gen_pt[2][i]),
                       s->bc[1], &(plaq_dest[i]));

    mult_nn_sum(&(src[0][i]), (matrix *)(gen_pt[3][i]), &(plaq_dest[i]));
  }
  cleanup_gather(tag[0]);
  cleanup_gather(tag[1]);
  cleanup_gather(tag[2]);
  cleanup_gather(tag[3]);
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use tempmat and tempmat2 for temporary storage
#ifdef VP
void Dminus(matrix *src, matrix *dest[NUMLINK]) {
  register int i;
  register site *s;
  char **local_pt[2][2];
  int mu, nu, flip = 0, opp_mu;
  msg_tag *tag0[2], *tag1[2];

  for (mu = 0; mu < 2; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[2 + mu];
  }

  // Start first set of gathers (nu = 0 and mu = 1)
  tag0[0] = start_gather_site(F_OFFSET(link[1]), sizeof(matrix),
                              goffset[0], EVENANDODD, local_pt[0][0]);

  FORALLSITES(i, s) {   // mu = 1 > nu = 0
    scalar_mult_nn(&(src[i]), &(s->link[1]), -1.0, &(tempmat[i]));
    FORALLDIR(mu)
      clear_mat(&(dest[mu][i]));        // Initialize
  }
  tag1[0] = start_gather_field(tempmat, sizeof(matrix),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]);

  // Main loop
  FORALLDIR(nu) {
    FORALLDIR(mu) {
      if (mu == nu)
        continue;

      if (nu == 0) {        // Start other set of gathers (nu = 1 and mu = 0)
        tag0[1] = start_gather_site(F_OFFSET(link[0]), sizeof(matrix),
                                    goffset[1], EVENANDODD, local_pt[1][0]);

        FORALLSITES(i, s)   // mu = 0 < nu = 1
          mult_nn(&(src[i]), &(s->link[0]), &(tempmat2[i]));
        tag1[1] = start_gather_field(tempmat2, sizeof(matrix),
                                     goffset[0] + 1, EVENANDODD, local_pt[1][1]);
      }

      opp_mu = OPP_LDIR(mu);
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu)      // src is anti-symmetric under mu <--> nu
          mult_nn_dif((matrix *)(local_pt[flip][0][i]), &(src[i]),
                      &(dest[nu][i]));
        else
          mult_nn_sum((matrix *)(local_pt[flip][0][i]), &(src[i]),
                      &(dest[nu][i]));

        scalar_mult_dif_matrix((matrix *)(local_pt[flip][1][i]),
                               s->bc[opp_mu], &(dest[nu][i]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip++;
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Term in action connecting site fermion to the link fermions
// bc[mu](x) on psi_mu(x) eta(x + mu)
// Add to dest instead of overwriting; note factor of 1/2
#ifdef SV
void DbplusStoL(matrix *src, matrix *dest[NUMLINK]) {
  register int i;
  register site *s;
  int mu;
  msg_tag *tag[NUMLINK];
  matrix tmat;

  tag[0] = start_gather_field(src, sizeof(matrix), goffset[0],
                              EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu == 0)     // Start other gather
      tag[1] = start_gather_field(src, sizeof(matrix), goffset[1],
                                  EVENANDODD, gen_pt[1]);

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      mult_na((matrix *)(gen_pt[mu][i]), &(s->link[mu]), &tmat);
      scalar_mult_matrix(&tmat, s->bc[mu], &tmat);
      mult_an_dif(&(s->link[mu]), &(src[i]), &tmat);
      scalar_mult_sum_matrix(&tmat, 0.5, &(dest[mu][i]));
    }
    cleanup_gather(tag[mu]);
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant coupling from site source to link destination
//   U^{-1}[a](x) * sum_b {D[b][a](x) + D[a][b](x - b)}
// D is Tr[eta] * plaqdet, Tr[eta] = i sqrt(N) eta^D
// T is Tr[U^{-1} Lambda]
// Assume compute_plaqdet() has already been run
// Use tr_dest and tempdet for temporary storage
// bc[b](x - b) = bc[-b](x) on eta(x - b) psi_a(x)
// Add negative to dest instead of overwriting
// Negative sign is due to anti-commuting eta past psi
#ifdef SV
void detStoL(matrix *dest[NUMLINK]) {
  register int i;
  register site *s;
  int a, b, opp_b;
  Real localG = -0.5 * C2 * G;
  complex tc;
  msg_tag *tag[NUMLINK];

  // Save Tr[eta(x)] plaqdet[a][b](x)
  //   or Tr[eta(x)] ZWstar[a][b](x) in tempdet[a][b]
  FORALLSITES(i, s) {
    CMUL(tr_eta[i], plaqdet[0][1][i], tempdet[0][1][i]);
    CMUL(tr_eta[i], plaqdet[1][0][i], tempdet[1][0][i]);
  }

  // Now we gather tempdet in both cases
  // Start first gather for (a, b) = (0, 1)
  tag[1] = start_gather_field(tempdet[0][1], sizeof(complex),
                              goffset[1] + 1, EVENANDODD, gen_pt[1]);

  FORALLDIR(a) {
    // Initialize accumulator for sum over b
    FORALLSITES(i, s)
      tr_dest[i] = cmplx(0.0, 0.0);

    FORALLDIR(b) {
      if (a == b)
        continue;

      if (a == 0) {      // Start other gather (a, b) = (1, 0)
        tag[0] = start_gather_field(tempdet[1][0], sizeof(complex),
                                    goffset[0] + 1, EVENANDODD, gen_pt[0]);
      }

      // Accumulate tempdet[b][a](x) + tempdet[a][b](x - b)
      opp_b = OPP_LDIR(b);
      wait_gather(tag[b]);
      FORALLSITES(i, s) {
        tc = *((complex *)(gen_pt[b][i]));
        tr_dest[i].real += s->bc[opp_b] * tc.real;
        tr_dest[i].imag += s->bc[opp_b] * tc.imag;
        CSUM(tr_dest[i], tempdet[b][a][i]);
      }
      cleanup_gather(tag[b]);
    }

    // Multiply U_a^{-1} by sum, add to dest[a][i]
    FORALLSITES(i, s) {
      CMULREAL(tr_dest[i], localG, tc);
      c_scalar_mult_sum_mat(&(Uinv[a][i]), &tc, &(dest[a][i]));
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Term in action connecting the link fermions to the site fermion
// Given src psi_a, dest is Dbar_a psi_a (Eq. 63 in the arXiv:1108.1503)
// Use tempmat and tempmat2 for temporary storage
// bc[OPP_LDIR(mu)](x) on eta(x - mu) psi_mu(x - mu)
// Initialize dest; note factor of 1/2
#ifdef SV
void DbminusLtoS(matrix *src[NUMLINK], matrix *dest) {
  register int i, mu, opp_mu;
  register site *s;
  msg_tag *tag[NUMLINK];

  FORALLSITES(i, s) {         // Set up first gather
    clear_mat(&(dest[i]));     // Initialize
    mult_an(&(s->link[0]), &(src[0][i]), &(tempmat[i]));
  }
  tag[0] = start_gather_field(tempmat, sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);

  FORALLDIR(mu) {
    if (mu == 0) {                  // Start other gather
      FORALLSITES(i, s)
        mult_an(&(s->link[1]), &(src[1][i]), &(tempmat2[i]));
      tag[1] = start_gather_field(tempmat2, sizeof(matrix),
                                  goffset[1] + 1, EVENANDODD, gen_pt[1]);
    }

    opp_mu = OPP_LDIR(mu);
    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_dif_matrix((matrix *)(gen_pt[mu][i]), s->bc[opp_mu],
                             &(dest[i]));
      mult_na_sum(&(src[mu][i]), &(s->link[mu]), &(dest[i]));
    }
    cleanup_gather(tag[mu]);
  }

  // Overall factor of 1/2
  FORALLSITES(i, s)
    scalar_mult_matrix(&(dest[i]), 0.5, &(dest[i]));
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant coupling from link source to site destination
//   sum_{a, b} D[a][b](x) * {T[b](x) +  T[a](x + b)}
// D is plaqdet and T is Tr[U^{-1} psi]
// Assume compute_plaqdet() has already been run
// bc[b](x) on eta(x) psi_a(x + b)
// Use Tr_Uinv and tr_dest for temporary storage
// Add to dest instead of overwriting
// Has same sign as DbminusLtoS (negative comes from generator normalization)
#ifdef SV
void detLtoS(matrix *src[NUMLINK], matrix *dest) {
  register int i;
  register site *s;
  int a, b;
  Real localG = 0.5 * C2 * G * sqrt((Real)NCOL);
  complex tc, tc2;
  msg_tag *tag[NUMLINK];

  // Prepare Tr[U_a^{-1} psi_a] = sum_j Tr[U_a^{-1} Lambda^j] psi_a^j
  // and save in Tr_Uinv[a]
  FORALLSITES(i, s) {
    tr_dest[i] = cmplx(0.0, 0.0);   // Initialize
    FORALLDIR(a)
      Tr_Uinv[a][i] = complextrace_nn(&(Uinv[a][i]), &(src[a][i]));
  }

  // Start first gather of Tr[U_a^{-1} psi_a] from x + b for (0, 1)
  tag[1] = start_gather_field(Tr_Uinv[0], sizeof(complex),
                              goffset[1], EVENANDODD, gen_pt[1]);

  // Main loop
  FORALLDIR(a) {
    FORALLDIR(b) {
      if (a == b)
        continue;

      if (a == 0) {         // Start other gather (1, 0)
        tag[0] = start_gather_field(Tr_Uinv[1], sizeof(complex),
                                    goffset[0], EVENANDODD, gen_pt[0]);
      }

      // Accumulate D[a][b](x) {T[b](x) + T[a](x + b)} in tr_dest
      wait_gather(tag[b]);
      FORALLSITES(i, s) {
        tc = *((complex *)(gen_pt[b][i]));
        tc2.real = Tr_Uinv[b][i].real + s->bc[b] * tc.real;
        tc2.imag = Tr_Uinv[b][i].imag + s->bc[b] * tc.imag;
        CMUL(plaqdet[a][b][i], tc2, tc);
        // localG is purely imaginary...
        tr_dest[i].real -= tc.imag * localG;
        tr_dest[i].imag += tc.real * localG;
      }
      cleanup_gather(tag[b]);
    }
  }

  // Add to dest (negative comes from generator normalization)
  FORALLSITES(i, s)
    c_scalar_mult_dif_mat(&(Lambda[DIMF - 1]), &(tr_dest[i]), &(dest[i]));
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Twist_Fermion matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
void fermion_op(Twist_Fermion *src, Twist_Fermion *dest, int sign) {
  register int i, mu;
  register site *s;

  // Copy src TwistFermion into fieldwise site, link and plaq fermions,
  // overwriting all of the latter
  if (sign == 1) {
    FORALLSITES(i, s) {
      mat_copy(&(src[i].Fsite), &(site_src[i]));
      FORALLDIR(mu)
        mat_copy(&(src[i].Flink[mu]), &(link_src[mu][i]));
      mat_copy(&(src[i].Fplaq), &(plaq_src[i]));
    }
  }
  else if (sign == -1) {
    FORALLSITES(i, s) {
      adjoint(&(src[i].Fsite), &(site_src[i]));
      FORALLDIR(mu)
        adjoint(&(src[i].Flink[mu]), &(link_src[mu][i]));
      adjoint(&(src[i].Fplaq), &(plaq_src[i]));
    }
  }
  else {
    node0_printf("Error: incorrect sign in fermion_op: %d\n", sign);
    terminate(1);
  }
  FORALLSITES(i, s)
    tr_eta[i] = trace(&(site_src[i]));

  // Assemble separate routines for each term in the fermion operator
#ifdef VP
  Dplus(link_src, plaq_dest);             // Overwrites plaq_dest
  Dminus(plaq_src, link_dest);            // Overwrites link_dest
#endif

#ifdef SV
  DbplusStoL(site_src, link_dest);        // Adds to link_dest

  // Site-to-link plaquette determinant contribution if G is non-zero
  // Only depends on Tr[eta(x)]
  if (doG)
    detStoL(link_dest);                   // Adds to link_dest

  DbminusLtoS(link_src, site_dest);       // Overwrites site_dest

  // Link-to-site plaquette determinant contribution if G is non-zero
  if (doG)
    detLtoS(link_src, site_dest);         // Adds to site_dest
#endif

  // Copy local plaquette, link and site fermions into dest TwistFermion
  if (sign == 1) {
    FORALLSITES(i, s) {
      mat_copy(&(site_dest[i]), &(dest[i].Fsite));
      FORALLDIR(mu)
        mat_copy(&(link_dest[mu][i]), &(dest[i].Flink[mu]));
      mat_copy(&(plaq_dest[i]), &(dest[i].Fplaq));
    }
  }
  else if (sign == -1) {    // Both negate and conjugate
    FORALLSITES(i, s) {
      neg_adjoint(&(site_dest[i]), &(dest[i].Fsite));
      FORALLDIR(mu)
        neg_adjoint(&(link_dest[mu][i]), &(dest[i].Flink[mu]));
      neg_adjoint(&(plaq_dest[i]), &(dest[i].Fplaq));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Squared Twist_Fermion matrix--vector operation
//   dest = (D^2 + fmass^2).src
// Use tempTF for temporary storage
void DSq(Twist_Fermion *src, Twist_Fermion *dest) {
  register int i;
  register site *s;

  fermion_op(src, tempTF, PLUS);
  fermion_op(tempTF, dest, MINUS);
  if (fmass > IMAG_TOL) {           // Assume fmass non-negative
    Real fmass2 = fmass * fmass;
    FORALLSITES(i, s)
      scalar_mult_sum_TF(&(src[i]), fmass2, &(dest[i]));
  }
}
// -----------------------------------------------------------------
