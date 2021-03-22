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
// This term connects the link fermions to the plaquette fermions
// Given src psi_a,
//   dest_ab(n) = bc[a] U_a(n) psi_b(n+a) - psi_b(n) U_a(n+b)
//              - bc[b] U_b(n) psi_a(n+b) + psi_a(n) U_b(n+a)
// Initialize dest; absorb factor of 1/2 by keeping b > a
#ifdef VP
void Dplus(matrix *src[NUMLINK], matrix *dest[NPLAQ]) {
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
  tag0[0] = start_gather_field(src[1], sizeof(matrix),
                               goffset[0], EVENANDODD, local_pt[0][0]);

  tag1[0] = start_gather_site(F_OFFSET(link[0]), sizeof(matrix),
                              goffset[1], EVENANDODD, local_pt[0][1]);

  tag2[0] = start_gather_field(src[0], sizeof(matrix),
                               goffset[1], EVENANDODD, local_pt[0][2]);

  tag3[0] = start_gather_site(F_OFFSET(link[1]), sizeof(matrix),
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
        tag0[gather] = start_gather_field(src[b], sizeof(matrix), goffset[a],
                                          EVENANDODD, local_pt[gather][0]);

        tag1[gather] = start_gather_site(F_OFFSET(link[a]), sizeof(matrix),
                                         goffset[b], EVENANDODD,
                                         local_pt[gather][1]);

        tag2[gather] = start_gather_field(src[a], sizeof(matrix), goffset[b],
                                          EVENANDODD, local_pt[gather][2]);

        tag3[gather] = start_gather_site(F_OFFSET(link[b]), sizeof(matrix),
                                         goffset[a], EVENANDODD,
                                         local_pt[gather][3]);
      }

      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      wait_gather(tag2[flip]);
      wait_gather(tag3[flip]);
      FORALLSITES(i, s) {
        // Initialize dest[index][i]
        scalar_mult_nn(&(s->link[mu]), (matrix *)(local_pt[flip][0][i]),
                       s->bc[mu], &(plaq_dest[index][i]));

        // Add or subtract the other three terms
        mult_nn_dif(&(src[nu][i]), (matrix *)(local_pt[flip][1][i]),
                    &(plaq_dest[index][i]));

        scalar_mult_nn_dif(&(s->link[nu]), (matrix *)(local_pt[flip][2][i]),
                           s->bc[nu], &(plaq_dest[index][i]));

        mult_nn_sum(&(src[mu][i]), (matrix *)(local_pt[flip][3][i]),
                    &(plaq_dest[index][i]));
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
// Term in action connecting the plaquette fermions to the link fermions
// Given src chi_ab,
//   dest_b(n) = U_a(n+b) chi_ab(n) - bc(opp_a) chi_ab(n-a) U_a(n-a)
// Initialize dest
// Use tempmat and tempmat2 for temporary storage
#ifdef VP
void Dminus(matrix *src[NPLAQ], matrix *dest[NUMLINK]) {
  register int i;
  register site *s;
  char **local_pt[2][2];
  int mu, nu, index, gather, flip = 0, a, b, next, opp_mu;
  matrix *mat[2];
  msg_tag *tag0[2], *tag1[2];

  for (mu = 0; mu < 2; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[2 + mu];
  }
  mat[0] = tempmat;
  mat[1] = tempmat2;

  // Start first set of gathers (mu = 1 and nu = 0)
  index = plaq_index[1][0];
  tag0[0] = start_gather_site(F_OFFSET(link[1]), sizeof(matrix),
                              goffset[0], EVENANDODD, local_pt[0][0]);

  FORALLSITES(i, s) {   // mu = 1 > nu = 0
    scalar_mult_nn(&(src[index][i]), &(s->link[1]), -1.0, &(mat[0][i]));
    FORALLDIR(mu)
      clear_mat(&(dest[mu][i]));        // Initialize
  }
  tag1[0] = start_gather_field(mat[0], sizeof(matrix),
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
        tag0[gather] = start_gather_site(F_OFFSET(link[a]), sizeof(matrix),
                                         goffset[b], EVENANDODD,
                                         local_pt[gather][0]);

        FORALLSITES(i, s) {
          if (a > b) {      // src is anti-symmetric under a <--> b
            scalar_mult_nn(&(src[next][i]), &(s->link[a]), -1.0,
                           &(mat[gather][i]));
          }
          else {
            mult_nn(&(src[next][i]), &(s->link[a]), &(mat[gather][i]));
          }
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix),
                                          goffset[a] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      opp_mu = OPP_LDIR(mu);
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu)      // src is anti-symmetric under mu <--> nu
          mult_nn_dif((matrix *)(local_pt[flip][0][i]), &(src[index][i]),
                      &(dest[nu][i]));
        else
          mult_nn_sum((matrix *)(local_pt[flip][0][i]), &(src[index][i]),
                      &(dest[nu][i]));

        scalar_mult_dif_matrix((matrix *)(local_pt[flip][1][i]),
                               s->bc[opp_mu], &(dest[nu][i]));
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
// Term in action connecting the plaquette fermions amongst themselves
// Given src chi_de and recalling a + b = -c - d - e,
//   dest_ab(n) = bc3[a][b][c] chi_de(n+a+b+c) Ubar_c(n+a+b)
//              - bc2[a][b] Ubar_c(n-c) chi_de(n+a+b)
// Add to dest instead of overwriting; note factor of 1/2
#ifdef QCLOSED
void DbplusPtoP(matrix *src[NPLAQ], matrix *dest[NPLAQ]) {
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
  tag0[0] = start_gather_site(F_OFFSET(link[c]), sizeof(matrix),
                              DbpP_d1[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(src[i_de], sizeof(matrix),
                               DbpP_d2[0], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(src[i_de], sizeof(matrix),
                               DbpP_d1[0], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_site(F_OFFSET(link[c]), sizeof(matrix),
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

      tag0[gather] = start_gather_site(F_OFFSET(link[c]), sizeof(matrix),
                                       DbpP_d1[next], EVENANDODD,
                                       local_pt[gather][0]);
      tag1[gather] = start_gather_field(src[i_de], sizeof(matrix),
                                        DbpP_d2[next], EVENANDODD,
                                        local_pt[gather][1]);
      tag2[gather] = start_gather_field(src[i_de], sizeof(matrix),
                                        DbpP_d1[next], EVENANDODD,
                                        local_pt[gather][2]);
      tag3[gather] = start_gather_site(F_OFFSET(link[c]), sizeof(matrix),
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
      scalar_mult_na_sum((matrix *)(local_pt[flip][1][i]),
                         (matrix *)(local_pt[flip][0][i]),
                         tr * s->bc3[a][b][c], &(dest[i_ab][i]));

      scalar_mult_an_dif((matrix *)(local_pt[flip][3][i]),
                         (matrix *)(local_pt[flip][2][i]),
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
// Term in action connecting the plaquette fermions amongst themselves
// Given src chi_ab and recalling d + e = -a - b - c,
//   dest_de(n) = bc2[opp_a][opp_b] chi_ab(n-a-b) Ubar_c(n-a-b-c)
//              - bc3[opp_a][opp_b][opp_c] Ubar_c(n-c) chi_ab(n-a-b-c)
// Add to dest instead of overwriting; note factor of 1/2
#ifdef QCLOSED
void DbminusPtoP(matrix *src[NPLAQ], matrix *dest[NPLAQ]) {
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

  tag0[0] = start_gather_site(F_OFFSET(link[c]), sizeof(matrix),
                              DbmP_d1[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(src[i_ab], sizeof(matrix),
                               DbmP_d2[0], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(src[i_ab], sizeof(matrix),
                               DbmP_d1[0], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_site(F_OFFSET(link[c]), sizeof(matrix),
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
      tag0[gather] = start_gather_site(F_OFFSET(link[c]), sizeof(matrix),
                                       DbmP_d1[next], EVENANDODD,
                                       local_pt[gather][0]);
      tag1[gather] = start_gather_field(src[i_ab], sizeof(matrix),
                                        DbmP_d2[next], EVENANDODD,
                                        local_pt[gather][1]);
      tag2[gather] = start_gather_field(src[i_ab], sizeof(matrix),
                                        DbmP_d1[next], EVENANDODD,
                                        local_pt[gather][2]);
      tag3[gather] = start_gather_site(F_OFFSET(link[c]), sizeof(matrix),
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
      scalar_mult_na_sum((matrix *)(local_pt[flip][1][i]),
                         (matrix *)(local_pt[flip][0][i]),
                         tr * s->bc2[opp_a][opp_b], &(dest[i_de][i]));

      scalar_mult_an_dif((matrix *)(local_pt[flip][3][i]),
                         (matrix *)(local_pt[flip][2][i]),
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
// Term in action connecting the site fermions to the link fermions
// Given src eta,
//   dest_a(n) = bc[a] eta(n+a) Ubar_a(n) - Ubar_a(n) eta(n)
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
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_field(src, sizeof(matrix), goffset[mu + 1],
                                       EVENANDODD, gen_pt[mu + 1]);

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
  int a, b, opp_b, next;
  Real localG = -0.5 * C2 * G;
  complex tc;
  msg_tag *tag[NUMLINK];

  // Save Tr[eta(x)] plaqdet[a][b](x)
  //   or Tr[eta(x)] ZWstar[a][b](x) in tempdet[a][b]
  FORALLDIR(a) {
    for (b = a + 1; b < NUMLINK; b++) {
      FORALLSITES(i, s) {
        CMUL(tr_eta[i], plaqdet[a][b][i], tempdet[a][b][i]);
        CMUL(tr_eta[i], plaqdet[b][a][i], tempdet[b][a][i]);
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
// Given src psi_a
//   dest(n) = psi_a(n) Ubar_a(n) - bc[opp_a] Ubar_a(n-a) psi_a(n-a)
// Initialize dest; note factor of 1/2
// Use tempmat and tempmat2 for temporary storage
#ifdef SV
void DbminusLtoS(matrix *src[NUMLINK], matrix *dest) {
  register int i, mu, nu, opp_mu;
  register site *s;
  int gather = 1, flip = 0;
  matrix *mat[2];
  msg_tag *tag[NUMLINK];

  mat[0] = tempmat;
  mat[1] = tempmat2;

  FORALLSITES(i, s) {           // Set up first gather
    clear_mat(&(dest[i]));      // Initialize
    mult_an(&(s->link[0]), &(src[0][i]), &(mat[0][i]));
  }
  tag[0] = start_gather_field(mat[0], sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);

  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {   // Start next gather
      nu = mu + 1;
      gather = (flip + 1) % 2;
      FORALLSITES(i, s)
        mult_an(&(s->link[nu]), &(src[nu][i]), &(mat[gather][i]));
      tag[nu] = start_gather_field(mat[gather], sizeof(matrix),
                                   goffset[nu] + 1, EVENANDODD, gen_pt[nu]);
    }

    opp_mu = OPP_LDIR(mu);
    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_dif_matrix((matrix *)(gen_pt[mu][i]), s->bc[opp_mu],
                             &(dest[i]));
      mult_na_sum(&(src[mu][i]), &(s->link[mu]), &(dest[i]));
    }
    cleanup_gather(tag[mu]);
    flip = gather;
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
  int a, b, next;
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
      for (mu = 0; mu < NPLAQ; mu++)
        mat_copy(&(src[i].Fplaq[mu]), &(plaq_src[mu][i]));
    }
  }
  else if (sign == -1) {
    FORALLSITES(i, s) {
      adjoint(&(src[i].Fsite), &(site_src[i]));
      FORALLDIR(mu)
        adjoint(&(src[i].Flink[mu]), &(link_src[mu][i]));
      for (mu = 0; mu < NPLAQ; mu++)
        adjoint(&(src[i].Fplaq[mu]), &(plaq_src[mu][i]));
    }
  }
  else {
    node0_printf("Error: incorrect sign in fermion_op: %d\n", sign);
    terminate(1);
  }
  FORALLSITES(i, s) {
    // Optionally rescale to make fermion operator more symmetric
#ifdef RESCALE
    scalar_mult_matrix(&(site_src[i]), 2.0, &(site_src[i]));
#endif
    tr_eta[i] = trace(&(site_src[i]));
  }

  // Assemble separate routines for each term in the fermion operator
#ifdef VP
  // Plaquette-to-link and link-to-plaquette contribution
  Dplus(link_src, plaq_dest);             // Overwrites plaq_dest
  Dminus(plaq_src, link_dest);            // Overwrites link_dest
#else
  FORALLSITES(i, s) {                     // Zero link_dest and plaq_dest
    FORALLDIR(mu)
      clear_mat(&(link_dest[mu][i]));
    for (mu = 0; mu < NPLAQ; mu++)
      clear_mat(&(plaq_dest[mu][i]));
  }
#endif

#ifdef SV
  // Site-to-link and link-to-site contribution
  DbplusStoL(site_src, link_dest);        // Adds to link_dest

  // Site-to-link plaquette determinant contribution if G is non-zero
  // Only depends on Tr[eta(x)]
  if (doG)
    detStoL(link_dest);                   // Adds to link_dest

  DbminusLtoS(link_src, site_dest);       // Overwrites site_dest

  // Link-to-site plaquette determinant contribution if G is non-zero
  if (doG)
    detLtoS(link_src, site_dest);         // Adds to site_dest
#else
  FORALLSITES(i, s)                       // Zero site_dest
    clear_mat(&(site_dest[i]));
#endif

#ifdef QCLOSED
  // Plaquette-to-plaquette contributions
  DbminusPtoP(plaq_src, plaq_dest);       // Adds to plaq_dest
  DbplusPtoP(plaq_src, plaq_dest);        // Adds to plaq_dest
#endif

  // Optionally rescale to make fermion operator more symmetric
#ifdef RESCALE
  FORALLSITES(i, s)
    scalar_mult_matrix(&(site_dest[i]), 2.0, &(site_dest[i]));
#endif

  // Copy local plaquette, link and site fermions into dest TwistFermion
  if (sign == 1) {
    FORALLSITES(i, s) {
      mat_copy(&(site_dest[i]), &(dest[i].Fsite));
      FORALLDIR(mu)
        mat_copy(&(link_dest[mu][i]), &(dest[i].Flink[mu]));
      for (mu = 0; mu < NPLAQ; mu++)
        mat_copy(&(plaq_dest[mu][i]), &(dest[i].Fplaq[mu]));
    }
  }
  else if (sign == -1) {    // Both negate and conjugate
    FORALLSITES(i, s) {
      neg_adjoint(&(site_dest[i]), &(dest[i].Fsite));
      FORALLDIR(mu)
        neg_adjoint(&(link_dest[mu][i]), &(dest[i].Flink[mu]));
      for (mu = 0; mu < NPLAQ; mu++)
        neg_adjoint(&(plaq_dest[mu][i]), &(dest[i].Fplaq[mu]));
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
