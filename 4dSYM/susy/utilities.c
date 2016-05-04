// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force

// #define DET_DIST prints out all determinants for plotting distribution
// #define TR_DIST prints out all (Tr[U.Udag]/N-1)^2 for plotting distribution
// CAUTION: Do not run DET_DIST or TR_DIST with MPI!

//#define DET_DIST
//#define TR_DIST
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
  int a, b, gather, flip = 0, mu, nu;
  complex tc;
  msg_tag *tag0[2], *tag1[2];

#ifdef DET_DIST
  if (this_node != 0) {
    printf("compute_plaqdet: don't run DET_DIST in parallel\n");
    fflush(stdout);
    terminate(1);
  }
#endif

  for (mu = 0; mu < 2; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[2 + mu];
  }

  // Gather determinants rather than the full matrices
  // Recall det[Udag] = (det[U])^*
  FORALLSITES(i, s) {
    FORALLDIR(a)
      Tr_Uinv[a][i] = find_det(&(s->linkf[a]));
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

      gather = (flip + 1) % 2;
      if (a < NUMLINK - 1 || b < NUMLINK - 2) { // Start next set of gathers
        if (b == NUMLINK - 1) {
          mu = a + 1;
          nu = 0;
        }
        else if (b == a - 1) {
          mu = a;
          nu = b + 2;
        }
        else {
          mu = a;
          nu = b + 1;
        }
        tag0[gather] = start_gather_field(Tr_Uinv[nu], sizeof(complex),
                                          goffset[mu], EVENANDODD,
                                          local_pt[gather][0]);
        tag1[gather] = start_gather_field(Tr_Uinv[mu], sizeof(complex),
                                          goffset[nu], EVENANDODD,
                                          local_pt[gather][1]);
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
          printf("DET_DIST %d %d %d %d %d %d %.4g %.4g %.4g\n",
                 s->x, s->y, s->z, s->t, a, b,
                 plaqdet[a][b][i].real, plaqdet[a][b][i].imag, cabs_sq(&tc1));
        }
#endif
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Save U_a^{-1} at each site
void compute_Uinv() {
  register int i, mu;
  register site *s;

  FORALLSITES(i, s) {
    FORALLDIR(mu)
      invert(&(s->linkf[mu]), &(Uinv[mu][i]));
  }
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Separate routines for each term in the fermion operator
// All called by fermion_op at the bottom of the file
#ifdef VP
void Dplus(matrix_f *src[NUMLINK], matrix_f *dest[NPLAQ]) {
  register int i;
  register site *s;
  char **local_pt[2][4];
  int mu, nu, index, gather, flip = 0, a, b;
  msg_tag *tag0[2], *tag1[2], *tag2[2], *tag3[2];

  for (mu = 0; mu < 4; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[4 + mu];
  }

  // Start first set of gathers (mu = 0 and nu = 1)
  tag0[0] = start_gather_field(src[1], sizeof(matrix_f),
                               goffset[0], EVENANDODD, local_pt[0][0]);

  tag1[0] = start_gather_site(F_OFFSET(linkf[0]), sizeof(matrix_f),
                              goffset[1], EVENANDODD, local_pt[0][1]);

  tag2[0] = start_gather_field(src[0], sizeof(matrix_f),
                               goffset[1], EVENANDODD, local_pt[0][2]);

  tag3[0] = start_gather_site(F_OFFSET(linkf[1]), sizeof(matrix_f),
                              goffset[0], EVENANDODD, local_pt[0][3]);

  // Main loop
  FORALLDIR(mu) {
    for (nu = mu + 1; nu < NUMLINK; nu++) {
      index = plaq_index[mu][nu];
      gather = (flip + 1) % 2;
      if (index < NPLAQ - 1) {               // Start next set of gathers
        if (nu == NUMLINK - 1) {
          a = mu + 1;
          b = a + 1;
        }
        else {
          a = mu;
          b = nu + 1;
        }
        tag0[gather] = start_gather_field(src[b], sizeof(matrix_f), goffset[a],
                                          EVENANDODD, local_pt[gather][0]);

        tag1[gather] = start_gather_site(F_OFFSET(linkf[a]), sizeof(matrix_f),
                                         goffset[b], EVENANDODD,
                                         local_pt[gather][1]);

        tag2[gather] = start_gather_field(src[a], sizeof(matrix_f), goffset[b],
                                          EVENANDODD, local_pt[gather][2]);

        tag3[gather] = start_gather_site(F_OFFSET(linkf[b]), sizeof(matrix_f),
                                         goffset[a], EVENANDODD,
                                         local_pt[gather][3]);
      }

      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      wait_gather(tag2[flip]);
      wait_gather(tag3[flip]);
      FORALLSITES(i, s) {
        // Initialize dest[index][i]
        scalar_mult_nn_f(&(s->linkf[mu]),
                         (matrix_f *)(local_pt[flip][0][i]),
                         s->bc1[mu], &(plaq_pmat[index][i]));

        // Add or subtract the other three terms
        mult_nn_dif_f(&(src[nu][i]), (matrix_f *)(local_pt[flip][1][i]),
                      &(plaq_pmat[index][i]));

        scalar_mult_nn_dif_f(&(s->linkf[nu]),
                             (matrix_f *)(local_pt[flip][2][i]),
                             s->bc1[nu], &(plaq_pmat[index][i]));

        mult_nn_sum_f(&(src[mu][i]), (matrix_f *)(local_pt[flip][3][i]),
                      &(plaq_pmat[index][i]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      cleanup_gather(tag2[flip]);
      cleanup_gather(tag3[flip]);
      flip = gather;
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use tempmat and tempmat2 for temporary storage
#ifdef VP
void Dminus(matrix_f *src[NPLAQ], matrix_f *dest[NUMLINK]) {
  register int i;
  register site *s;
  char **local_pt[2][2];
  int mu, nu, index, gather, flip = 0, a, b, next, opp_mu;
  matrix_f *mat[2];
  msg_tag *tag0[2], *tag1[2];

  for (mu = 0; mu < 2; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[2 + mu];
  }
  mat[0] = tempmat;
  mat[1] = tempmat2;

  // Start first set of gathers (mu = 1 and nu = 0)
  index = plaq_index[1][0];
  tag0[0] = start_gather_site(F_OFFSET(linkf[1]), sizeof(matrix_f),
                              goffset[0], EVENANDODD, local_pt[0][0]);

  FORALLSITES(i, s) {   // mu = 1 > nu = 0
    scalar_mult_nn_f(&(src[index][i]), &(s->linkf[1]), -1.0, &(mat[0][i]));
    FORALLDIR(mu)
      clear_mat_f(&(dest[mu][i]));        // Initialize
  }
  tag1[0] = start_gather_field(mat[0], sizeof(matrix_f),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]);

  // Main loop
  FORALLDIR(nu) {
    FORALLDIR(mu) {
      if (mu == nu)
        continue;

      gather = (flip + 1) % 2;
      if (nu < NUMLINK - 1 || mu < NUMLINK - 2) { // Start next set of gathers
        if (mu == NUMLINK - 1) {
          a = 0;
          b = nu + 1;
        }
        else if (mu == nu - 1) {
          a = mu + 2;
          b = nu;
        }
        else {
          a = mu + 1;
          b = nu;
        }
        next = plaq_index[a][b];
        tag0[gather] = start_gather_site(F_OFFSET(linkf[a]), sizeof(matrix_f),
                                         goffset[b], EVENANDODD,
                                         local_pt[gather][0]);

        FORALLSITES(i, s) {
          if (a > b) {      // src is anti-symmetric under a <--> b
            scalar_mult_nn_f(&(src[next][i]), &(s->linkf[a]), -1.0,
                             &(mat[gather][i]));
          }
          else {
            mult_nn_f(&(src[next][i]), &(s->linkf[a]), &(mat[gather][i]));
          }
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix_f),
                                          goffset[a] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      opp_mu = OPP_LDIR(mu);
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu)      // src is anti-symmetric under mu <--> nu
          mult_nn_dif_f((matrix_f *)(local_pt[flip][0][i]), &(src[index][i]),
                        &(dest[nu][i]));
        else
          mult_nn_sum_f((matrix_f *)(local_pt[flip][0][i]), &(src[index][i]),
                        &(dest[nu][i]));

        scalar_mult_dif_matrix_f((matrix_f *)(local_pt[flip][1][i]),
                                 s->bc1[opp_mu], &(dest[nu][i]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Add to dest instead of overwriting; note factor of 1/2
#ifdef QCLOSED
void DbplusPtoP(matrix_f *src[NPLAQ], matrix_f *dest[NPLAQ]) {
  register int i;
  register site *s;
  char **local_pt[2][4];
  int a, b, c, d, e, j, gather, next, flip = 0, i_ab, i_de;
  Real tr;
  msg_tag *tag0[2], *tag1[2], *tag2[2], *tag3[2];

  for (a = 0; a < 4; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[4 + a];
  }

  // Start first set of gathers
  // From setup_lambda.c, we see b > a and e > d
  c = DbplusPtoP_lookup[0][2];
  d = DbplusPtoP_lookup[0][3];
  e = DbplusPtoP_lookup[0][4];
  i_de = plaq_index[d][e];
  tag0[0] = start_gather_site(F_OFFSET(linkf[c]), sizeof(matrix_f),
                              DbpP_d1[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(src[i_de], sizeof(matrix_f),
                               DbpP_d2[0], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(src[i_de], sizeof(matrix_f),
                               DbpP_d1[0], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_site(F_OFFSET(linkf[c]), sizeof(matrix_f),
                              goffset[c] + 1, EVENANDODD, local_pt[0][3]);

  // Loop over lookup table
  for (j = 0; j < NTERMS; j++) {
    gather = (flip + 1) % 2;
    if (j < NTERMS - 1) {               // Start next set of gathers
      next = j + 1;
      c = DbplusPtoP_lookup[next][2];
      d = DbplusPtoP_lookup[next][3];
      e = DbplusPtoP_lookup[next][4];
      i_de = plaq_index[d][e];

      tag0[gather] = start_gather_site(F_OFFSET(linkf[c]), sizeof(matrix_f),
                                       DbpP_d1[next], EVENANDODD,
                                       local_pt[gather][0]);
      tag1[gather] = start_gather_field(src[i_de], sizeof(matrix_f),
                                        DbpP_d2[next], EVENANDODD,
                                        local_pt[gather][1]);
      tag2[gather] = start_gather_field(src[i_de], sizeof(matrix_f),
                                        DbpP_d1[next], EVENANDODD,
                                        local_pt[gather][2]);
      tag3[gather] = start_gather_site(F_OFFSET(linkf[c]), sizeof(matrix_f),
                                       goffset[c] + 1, EVENANDODD,
                                       local_pt[gather][3]);
    }

    // Do this set of computations while next set of gathers runs
    a = DbplusPtoP_lookup[j][0];
    b = DbplusPtoP_lookup[j][1];
    c = DbplusPtoP_lookup[j][2];
    d = DbplusPtoP_lookup[j][3];
    e = DbplusPtoP_lookup[j][4];
    tr = 0.5 * perm[a][b][c][d][e];
    i_ab = plaq_index[a][b];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);
    wait_gather(tag3[flip]);
    FORALLSITES(i, s) {
      scalar_mult_na_sum_f((matrix_f *)(local_pt[flip][1][i]),
                           (matrix_f *)(local_pt[flip][0][i]),
                           tr * s->bc3[a][b][c], &(dest[i_ab][i]));

      scalar_mult_an_dif_f((matrix_f *)(local_pt[flip][3][i]),
                           (matrix_f *)(local_pt[flip][2][i]),
                           tr * s->bc2[a][b], &(dest[i_ab][i]));
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    cleanup_gather(tag3[flip]);
    flip = gather;
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Add to dest instead of overwriting; note factor of 1/2
#ifdef QCLOSED
void DbminusPtoP(matrix_f *src[NPLAQ], matrix_f *dest[NPLAQ]) {
  register int i, opp_a, opp_b, opp_c;
  register site *s;
  char **local_pt[2][4];
  int a, b, c, d, e, j, gather, next, flip = 0, i_ab, i_de;
  Real tr;
  msg_tag *tag0[2], *tag1[2], *tag2[2], *tag3[2];

  for (a = 0; a < 4; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[4 + a];
  }

  // Start first set of gathers
  // From setup_lambda.c, we see b > a and e > d
  a = DbminusPtoP_lookup[0][0];
  b = DbminusPtoP_lookup[0][1];
  c = DbminusPtoP_lookup[0][2];
  i_ab = plaq_index[a][b];

  tag0[0] = start_gather_site(F_OFFSET(linkf[c]), sizeof(matrix_f),
                              DbmP_d1[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(src[i_ab], sizeof(matrix_f),
                               DbmP_d2[0], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(src[i_ab], sizeof(matrix_f),
                               DbmP_d1[0], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_site(F_OFFSET(linkf[c]), sizeof(matrix_f),
                              goffset[c] + 1, EVENANDODD, local_pt[0][3]);

  // Loop over lookup table
  for (j = 0; j < NTERMS; j++) {
    gather = (flip + 1) % 2;
    if (j < NTERMS - 1) {               // Start next set of gathers
      next = j + 1;
      a = DbminusPtoP_lookup[next][0];
      b = DbminusPtoP_lookup[next][1];
      c = DbminusPtoP_lookup[next][2];
      i_ab = plaq_index[a][b];
      tag0[gather] = start_gather_site(F_OFFSET(linkf[c]), sizeof(matrix_f),
                                       DbmP_d1[next], EVENANDODD,
                                       local_pt[gather][0]);
      tag1[gather] = start_gather_field(src[i_ab], sizeof(matrix_f),
                                        DbmP_d2[next], EVENANDODD,
                                        local_pt[gather][1]);
      tag2[gather] = start_gather_field(src[i_ab], sizeof(matrix_f),
                                        DbmP_d1[next], EVENANDODD,
                                        local_pt[gather][2]);
      tag3[gather] = start_gather_site(F_OFFSET(linkf[c]), sizeof(matrix_f),
                                       goffset[c] + 1, EVENANDODD,
                                       local_pt[gather][3]);
    }

    // Do this set of computations while next set of gathers runs
    a = DbminusPtoP_lookup[j][0];
    b = DbminusPtoP_lookup[j][1];
    c = DbminusPtoP_lookup[j][2];
    d = DbminusPtoP_lookup[j][3];
    e = DbminusPtoP_lookup[j][4];
    tr = 0.5 * perm[a][b][c][d][e];
    i_de = plaq_index[d][e];

    opp_a = OPP_LDIR(a);
    opp_b = OPP_LDIR(b);
    opp_c = OPP_LDIR(c);
    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);
    wait_gather(tag3[flip]);
    FORALLSITES(i, s) {
      scalar_mult_na_sum_f((matrix_f *)(local_pt[flip][1][i]),
                           (matrix_f *)(local_pt[flip][0][i]),
                           tr * s->bc2[opp_a][opp_b], &(dest[i_de][i]));

      scalar_mult_an_dif_f((matrix_f *)(local_pt[flip][3][i]),
                           (matrix_f *)(local_pt[flip][2][i]),
                           tr * s->bc3[opp_a][opp_b][opp_c], &(dest[i_de][i]));
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    cleanup_gather(tag3[flip]);
    flip = gather;
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Term in action connecting site fermion to the link fermions
// bc1[mu](x) on psi_mu(x) eta(x + mu)
// Add to dest instead of overwriting; note factor of 1/2
#ifdef SV
void DbplusStoL(matrix_f *src, matrix_f *dest[NUMLINK]) {
  register int i;
  register site *s;
  int mu;
  msg_tag *tag[NUMLINK];
  matrix_f tmat;

  tag[0] = start_gather_field(src, sizeof(matrix_f), goffset[0],
                              EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_field(src, sizeof(matrix_f), goffset[mu + 1],
                                       EVENANDODD, gen_pt[mu + 1]);

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      mult_na_f((matrix_f *)(gen_pt[mu][i]), &(s->linkf[mu]), &tmat);
      scalar_mult_matrix_f(&tmat, s->bc1[mu], &tmat);
      mult_an_dif_f(&(s->linkf[mu]), &(src[i]), &tmat);
      scalar_mult_sum_matrix_f(&tmat, 0.5, &(dest[mu][i]));
    }
    cleanup_gather(tag[mu]);
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant coupling from site source to link destination
//   U^{-1}[a](x) * sum_b {D[b][a](x) + D[a][b](x - b)}
// In the global case D is Tr[eta] * plaqdet, Tr[eta] = i sqrt(N) eta^D
// In the local case D is Tr[eta] plaqdet (plaqdet - 1)^*
// In both cases T is Tr[U^{-1} Lambda] and
// Assume compute_plaqdet() has already been run
// Use tr_dest and tempdet for temporary storage
// bc1[b](x - b) = bc1[-b](x) on eta(x - b) psi_a(x)
// Add negative to dest instead of overwriting
// Negative sign is due to anti-commuting eta past psi
#ifdef SV
void detStoL(matrix_f *dest[NUMLINK]) {
  register int i;
  register site *s;
  int a, b, opp_b, next;
  complex tc, localGc;
#ifdef LINEAR_DET
  CMULREAL(Gc, -0.5, localGc);            // Since not squared
#else
  CMULREAL(Gc, -1.0, localGc);
#endif
  msg_tag *tag[NUMLINK];

  // Save Tr[eta(x)] plaqdet[a][b](x)
  //   or Tr[eta(x)] ZWstar[a][b](x) in tempdet[a][b]
  FORALLDIR(a) {
    for (b = a + 1; b < NUMLINK; b++) {
      FORALLSITES(i, s) {
#ifdef LINEAR_DET
        CMUL(tr_eta[i], plaqdet[a][b][i], tempdet[a][b][i]);
        CMUL(tr_eta[i], plaqdet[b][a][i], tempdet[b][a][i]);
#else
        CMUL(tr_eta[i], ZWstar[a][b][i], tempdet[a][b][i]);
        CMUL(tr_eta[i], ZWstar[b][a][i], tempdet[b][a][i]);
#endif
      }
    }
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

      // Start next gather unless we're doing the last (a=4, b=3)
      next = b + 1;
      if (next < NUMLINK && a + b < 2 * NUMLINK - 3) {
        if (next == a)              // Next gather is actually (a, b + 2)
          next++;

        tag[next] = start_gather_field(tempdet[a][next], sizeof(complex),
                                       goffset[next] + 1, EVENANDODD,
                                       gen_pt[next]);
      }
      else if (next == NUMLINK) {   // Start next gather (a + 1, 0)
        tag[0] = start_gather_field(tempdet[a + 1][0], sizeof(complex),
                                    goffset[0] + 1, EVENANDODD, gen_pt[0]);
      }

      // Accumulate tempdet[b][a](x) + tempdet[a][b](x - b)
      opp_b = OPP_LDIR(b);
      wait_gather(tag[b]);
      FORALLSITES(i, s) {
        tc = *((complex *)(gen_pt[b][i]));
        tr_dest[i].real += s->bc1[opp_b] * tc.real;
        tr_dest[i].imag += s->bc1[opp_b] * tc.imag;
        CSUM(tr_dest[i], tempdet[b][a][i]);
      }
      cleanup_gather(tag[b]);
    }

    // Multiply U_a^{-1} by sum, add to dest[a][i]
    FORALLSITES(i, s) {
      CMUL(tr_dest[i], localGc, tc);
      c_scalar_mult_sum_mat_f(&(Uinv[a][i]), &tc, &(dest[a][i]));
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar potential coupling from site source to link destination
//   Udag[a](x) Tr[eta(x)] (Tr[U_a(x) Udag_a(x)] / N - 1)
// Add negative to dest instead of overwriting
// Negative sign is due to anti-commuting eta past psi
#ifdef SV
void potStoL(matrix_f *dest[NUMLINK]) {
  register int i, a;
  register site *s;
  Real tr;
  complex tc;

  FORALLSITES(i, s) {
    FORALLDIR(a) {
      tr = 1.0 - one_ov_N * realtrace_f(&(s->linkf[a]), &(s->linkf[a]));
      CMUL(tr_eta[i], Bc, tc);
      CMULREAL(tc, tr, tc);
      c_scalar_mult_sum_mat_adj_f(&(s->linkf[a]), &tc, &(dest[a][i]));
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Term in action connecting the link fermions to the site fermion
// Given src psi_a, dest is Dbar_a psi_a (Eq. 63 in the arXiv:1108.1503)
// Use tempmat and tempmat2 for temporary storage
// bc1[OPP_LDIR(mu)](x) on eta(x - mu) psi_mu(x - mu)
// Initialize dest; note factor of 1/2
#ifdef SV
void DbminusLtoS(matrix_f *src[NUMLINK], matrix_f *dest) {
  register int i, mu, nu, opp_mu;
  register site *s;
  int gather = 1, flip = 0;
  matrix_f *mat[2];
  msg_tag *tag[NUMLINK];

  mat[0] = tempmat;
  mat[1] = tempmat2;

  FORALLSITES(i, s) {         // Set up first gather
    clear_mat_f(&(dest[i]));     // Initialize
    mult_an_f(&(s->linkf[0]), &(src[0][i]), &(mat[0][i]));
  }
  tag[0] = start_gather_field(mat[0], sizeof(matrix_f),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);

  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {   // Start next gather
      nu = mu + 1;
      gather = (flip + 1) % 2;
      FORALLSITES(i, s)
        mult_an_f(&(s->linkf[nu]), &(src[nu][i]), &(mat[gather][i]));
      tag[nu] = start_gather_field(mat[gather], sizeof(matrix_f),
                                   goffset[nu] + 1, EVENANDODD, gen_pt[nu]);
    }

    opp_mu = OPP_LDIR(mu);
    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_dif_matrix_f((matrix_f *)(gen_pt[mu][i]), s->bc1[opp_mu],
                               &(dest[i]));
      mult_na_sum_f(&(src[mu][i]), &(s->linkf[mu]), &(dest[i]));
    }
    cleanup_gather(tag[mu]);
    flip = gather;
  }

  // Overall factor of 1/2
  FORALLSITES(i, s)
    scalar_mult_matrix_f(&(dest[i]), 0.5, &(dest[i]));
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant coupling from link source to site destination
//   sum_{a, b} D[a][b](x) * {T[b](x) +  T[a](x + b)}
// In the global case D is plaqdet
// In the local case D is plaqdet (plaqdet - 1)^*
// In both cases T is Tr[U^{-1} psi]
// Assume compute_plaqdet() has already been run
// bc1[b](x) on eta(x) psi_a(x + b)
// Use Tr_Uinv and tr_dest for temporary storage
// Add to dest instead of overwriting
// Has same sign as DbminusLtoS (negative comes from generator normalization)
#ifdef SV
void detLtoS(matrix_f *src[NUMLINK], matrix_f *dest) {
  register int i;
  register site *s;
  int a, b, next;
  complex tc, tc2, localGc;
#ifdef LINEAR_DET
  CMULREAL(Gc, 0.5, localGc);                  // Since not squared
#endif
  msg_tag *tag[NUMLINK];

  // Prepare Tr[U_a^{-1} psi_a] = sum_j Tr[U_a^{-1} Lambda^j] psi_a^j
  // and save in Tr_Uinv[a]
  FORALLSITES(i, s) {
    tr_dest[i] = cmplx(0.0, 0.0);   // Initialize
    FORALLDIR(a)
      Tr_Uinv[a][i] = complextrace_nn_f(&(Uinv[a][i]), &(src[a][i]));
  }

  // Start first gather of Tr[U_a^{-1} psi_a] from x + b for (0, 1)
  tag[1] = start_gather_field(Tr_Uinv[0], sizeof(complex),
                              goffset[1], EVENANDODD, gen_pt[1]);

  // Main loop
  FORALLDIR(a) {
    FORALLDIR(b) {
      if (a == b)
        continue;

      // Start next gather unless we're doing the last (a=4, b=3)
      next = b + 1;
      if (next < NUMLINK && a + b < 2 * NUMLINK - 3) {
        if (next == a)              // Next gather is actually (a, b + 2)
          next++;
        tag[next] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                       goffset[next], EVENANDODD,
                                       gen_pt[next]);
      }
      else if (next == NUMLINK) {   // Start next gather (a + 1, 0)
        tag[0] = start_gather_field(Tr_Uinv[a + 1], sizeof(complex),
                                    goffset[0], EVENANDODD, gen_pt[0]);
      }

      // Accumulate D[a][b](x) {T[b](x) + T[a](x + b)} in tr_dest
      wait_gather(tag[b]);
      FORALLSITES(i, s) {
        tc = *((complex *)(gen_pt[b][i]));
        tc2.real = Tr_Uinv[b][i].real + s->bc1[b] * tc.real;
        tc2.imag = Tr_Uinv[b][i].imag + s->bc1[b] * tc.imag;
#ifdef LINEAR_DET
        CMUL(plaqdet[a][b][i], tc2, tc);
#else
        CMUL(ZWstar[a][b][i], tc2, tc);
#endif
        tr_dest[i].real += tc.real * localGc.real - tc.imag * localGc.imag;
        tr_dest[i].imag += tc.imag * localGc.real + tc.real * localGc.imag;
      }
      cleanup_gather(tag[b]);
    }
  }

  // Add to dest (negative comes from generator normalization)
  FORALLSITES(i, s)
    c_scalar_mult_dif_mat_f(&(Lambda[DIMF - 1]), &(tr_dest[i]), &(dest[i]));
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar potential coupling from link source to site destination
//   sum_a (Tr[U_a(x) Udag_a(x)] / N - 1)^2 psi(x) Udag[a](x)
// Use tr_dest for temporary storage
// Add to dest instead of overwriting
// Has same sign as DbminusLtoS (negative comes from generator normalization)
#ifdef SV
void potLtoS(matrix_f *src[NUMLINK], matrix_f *dest) {
  register int i, a;
  register site *s;
  Real tr;
  complex tc;

  FORALLSITES(i, s) {
    // Initialize tr_dest
    tr = one_ov_N * realtrace_f(&(s->linkf[0]), &(s->linkf[0])) - 1.0;
    tr_dest[i] = complextrace_na_f(&(src[0][i]), &(s->linkf[0]));
    CMULREAL(tr_dest[i], tr, tr_dest[i]);
    for (a = 1; a < NUMLINK; a++) {
      tr = one_ov_N * realtrace_f(&(s->linkf[a]), &(s->linkf[a])) - 1.0;
      tc = complextrace_na_f(&(src[a][i]), &(s->linkf[a]));
      tr_dest[i].real += tr * tc.real;
      tr_dest[i].imag += tr * tc.imag;
    }

    // Add to dest (negative comes from generator normalization)
    CMUL(tr_dest[i], Bc, tc);
    c_scalar_mult_dif_mat_f(&(Lambda[DIMF - 1]), &tc, &(dest[i]));
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Twist_Fermion matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
void fermion_op(Twist_Fermion *src, Twist_Fermion *dest, int sign) {
  register int i, j, mu, nu, index;
  register site *s;

  // Copy src TwistFermion into fieldwise site, link and plaq fermions
  // All of the latter are overwritten -- don't need to clear explicitly
  if (sign == 1) {
    FORALLSITES(i, s) {
      reconstruct(&(src[i].Fsite), &(site_mat[i]));
      tr_eta[i] = src[i].Fsite.c[DIMF - 1];
      for (j = 0; j < DIMF; j++)
        site_src[i].c[j] = src[i].Fsite.c[j];
      FORALLDIR(mu) {
        reconstruct(&(src[i].Flink[mu]), &(link_mat[mu][i]));
        for (j = 0; j < DIMF; j++)
          link_src[mu][i].c[j] = src[i].Flink[mu].c[j];
      }
      for (mu = 0; mu < NPLAQ; mu++) {
        reconstruct(&(src[i].Fplaq[mu]), &(plaq_mat[mu][i]));
        for (j = 0; j < DIMF; j++)
          plaq_src[mu][i].c[j] = src[i].Fplaq[mu].c[j];
      }
    }
  }
  else if (sign == -1) {
    FORALLSITES(i, s) {
      reconstruct_star(&(src[i].Fsite), &(site_mat[i]));
      CONJG(src[i].Fsite.c[DIMF - 1], tr_eta[i]);
      for (j = 0; j < DIMF; j++)
        CONJG(src[i].Fsite.c[j], site_src[i].c[j]);
      FORALLDIR(mu) {
        reconstruct_star(&(src[i].Flink[mu]), &(link_mat[mu][i]));
        for (j = 0; j < DIMF; j++)
          CONJG(src[i].Flink[mu].c[j], link_src[mu][i].c[j]);
      }
      for (mu = 0; mu < NPLAQ; mu++) {
        reconstruct_star(&(src[i].Fplaq[mu]), &(plaq_mat[mu][i]));
        for (j = 0; j < DIMF; j++)
          CONJG(src[i].Fplaq[mu].c[j], plaq_src[mu][i].c[j]);
      }
    }
  }
  else {
    node0_printf("Error: incorrect sign in fermion_op: %d\n", sign);
    terminate(1);
  }

  // Assemble separate routines for each term in the fermion operator
#ifdef VP
  Dplus(link_mat, plaq_pmat);             // Initializes plaq_pmat
  Dminus(plaq_mat, link_pmat);            // Initializes link_pmat
#endif

#ifdef SV
  DbplusStoL(site_mat, link_pmat);        // Adds to link_pmat

  // Site-to-link plaquette determinant contribution if G is non-zero
  // Only depends on Tr[eta(x)]
  if (doG)
    detStoL(link_pmat);                   // Adds to link_pmat

  // Site-to-link scalar potential contribution if B is non-zero
  // Only depends on Tr[eta(x)]
  if (doB)
    potStoL(link_pmat);                   // Adds to link_pmat

  DbminusLtoS(link_mat, site_pmat);       // Initializes site_pmat

  // Link-to-site plaquette determinant contribution if G is non-zero
  if (doG)
    detLtoS(link_mat, site_pmat);         // Adds to site_dest

  // Link-to-site scalar potential contribution if B is non-zero
  if (doB)
    potLtoS(link_mat, site_pmat);         // Adds to site_dest
#endif

#ifdef QCLOSED
  DbminusPtoP(plaq_mat, plaq_pmat);       // Adds to plaq_dest
  DbplusPtoP(plaq_mat, plaq_pmat);        // Adds to plaq_dest
#endif
  FORALLSITES(i, s) {
    deconstruct(&(site_pmat[i]), &(site_dest[i]));
    FORALLDIR(mu) {
      deconstruct(&(link_pmat[mu][i]), &(link_dest[mu][i]));
      for (nu = mu + 1; nu < NUMLINK; nu++) {
        index = plaq_index[mu][nu];
        deconstruct(&(plaq_pmat[index][i]), &(plaq_dest[index][i]));
      }
    }
  }
  // Copy local plaquette, link and site fermions into dest TwistFermion
  if (sign == 1) {
    FORALLSITES(i, s) {
      for (j = 0; j < DIMF; j++)
        dest[i].Fsite.c[j] = site_dest[i].c[j];
      FORALLDIR(mu) {
        for (j = 0; j < DIMF; j++)
          dest[i].Flink[mu].c[j] = link_dest[mu][i].c[j];
      }
      for (mu = 0; mu < NPLAQ; mu++) {
        for (j = 0; j < DIMF; j++)
          dest[i].Fplaq[mu].c[j] = plaq_dest[mu][i].c[j];
      }
    }
  }
  else if (sign == -1) {    // Both negate and conjugate
    FORALLSITES(i, s) {
      for (j = 0; j < DIMF; j++) {
        dest[i].Fsite.c[j].real = -site_dest[i].c[j].real;
        dest[i].Fsite.c[j].imag = site_dest[i].c[j].imag;
      }
      FORALLDIR(mu) {
        for (j = 0; j < DIMF; j++) {
          dest[i].Flink[mu].c[j].real = -link_dest[mu][i].c[j].real;
          dest[i].Flink[mu].c[j].imag = link_dest[mu][i].c[j].imag;
        }
      }
      for (mu = 0; mu < NPLAQ; mu++) {
        for (j = 0; j < DIMF; j++) {
          dest[i].Fplaq[mu].c[j].real = -plaq_dest[mu][i].c[j].real;
          dest[i].Fplaq[mu].c[j].imag = plaq_dest[mu][i].c[j].imag;
        }
      }
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
  if (fmass > IMAG_TOL) {
    Real fmass2 = fmass * fmass;
    FORALLSITES(i, s)
      scalar_mult_sum_TF(&(src[i]), fmass2, &(dest[i]));
  }
}
// -----------------------------------------------------------------
