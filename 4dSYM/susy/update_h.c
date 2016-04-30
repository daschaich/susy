// -----------------------------------------------------------------
// Update the momentum matrices
// Uncomment to print out debugging messages
//#define FORCE_DEBUG
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update mom[NUMLINK] with the gauge force
// Include tunable coefficient C2 in the d^2 term of the action
// Use tr_dest, tempmat and tempdet for temporary storage
// Assume compute_plaqdet(), compute_DmuUmu()
// and compute_Fmunu() have already been run
double gauge_force(Real eps) {
  register int i, mu, nu;
  register site *s;
  char **local_pt[2][2];
  int a, b, gather, flip = 0, index, next;
  double returnit = 0.0;
  complex tc, tc2;
  msg_tag *tag[NUMLINK], *tag0[2], *tag1[2];
  matrix_f tmat, tmat2, *mat[2];

  // Three contributions from d^2 term
  // All three terms need a factor of C2
  // First we have the finite difference operator derivative times DmuUmu
  // Ubar_a(x) DmuUmu(x) - DmuUmu(x + a) Ubar_a(x)
  tag[0] = start_gather_field(DmuUmu, sizeof(matrix_f),
                              goffset[0], EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_field(DmuUmu, sizeof(matrix_f),
                                       goffset[mu + 1], EVENANDODD,
                                       gen_pt[mu + 1]);

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      mult_an_f(&(s->linkf[mu]), &(DmuUmu[i]), &tmat);
      mult_na_f((matrix_f *)(gen_pt[mu][i]), &(s->linkf[mu]), &tmat2);
      sub_matrix_f(&tmat, &tmat2, &(s->f_U[mu]));   // Initialize
    }
    cleanup_gather(tag[mu]);
  }

  // Second we have the plaquette determinant derivative contribution
  //   U_mu^{-1}(x) 2G sum_nu {D[nu][mu](x) + D[mu][nu](x-nu)}
  // In the global case D is Tr[DmuUmu] plaqdet[mu][nu]
  // In the local case D is 2Tr[DmuUmu] ZWstar[mu][nu]
  // In both cases we save D in tempdet[mu][nu]
  // Only compute if G is non-zero
  // Use tr_dest for temporary storage
  if (doG) {
    FORALLSITES(i, s) {
      tc = trace_f(&DmuUmu[i]);
      FORALLDIR(mu) {
        for (nu = mu + 1; nu < NUMLINK; nu++) {
#ifdef LINEAR_DET
          CMUL(tc, plaqdet[mu][nu][i], tempdet[mu][nu][i]);
          CMUL(tc, plaqdet[nu][mu][i], tempdet[nu][mu][i]);
#else
          CMUL(tc, ZWstar[mu][nu][i], tempdet[mu][nu][i]);
          CMUL(tc, ZWstar[nu][mu][i], tempdet[nu][mu][i]);
#endif
        }
      }
    }

    // Start first gather (mu = 0 and nu = 1), labelled by nu
    // Gather D[mu][nu] from x - nu
    tag[1] = start_gather_field(tempdet[0][1], sizeof(complex),
                                goffset[1] + 1, EVENANDODD, gen_pt[1]);

    // Main loop
    FORALLDIR(mu) {
      // Zero tr_dest to hold sum
      FORALLSITES(i, s)
        tr_dest[i] = cmplx(0.0, 0.0);

      FORALLDIR(nu) {
        if (mu == nu)
          continue;

        if (mu < NUMLINK - 1 || nu < NUMLINK - 2) { // Start next gather
          if (nu == NUMLINK - 1) {
            a = mu + 1;
            b = 0;
          }
          else if (nu == mu - 1) {
            a = mu;
            b = nu + 2;
          }
          else {
            a = mu;
            b = nu + 1;
          }
          tag[b] = start_gather_field(tempdet[a][b], sizeof(complex),
                                      goffset[b] + 1, EVENANDODD,
                                      gen_pt[b]);
        }

        // Add D[nu][mu](x) to sum while gather runs
        FORALLSITES(i, s)
          CSUM(tr_dest[i], tempdet[nu][mu][i]);

        // Add D[mu][nu](x - nu) to sum
        wait_gather(tag[nu]);
        FORALLSITES(i, s)
          CSUM(tr_dest[i], *((complex *)(gen_pt[nu][i])));
        cleanup_gather(tag[nu]);
      }

      // Now add to force
      FORALLSITES(i, s) {
        invert(&(s->linkf[mu]), &tmat);
#ifdef LINEAR_DET
        CMULREAL(tr_dest[i], G, tc);
#else
        CMULREAL(tr_dest[i], 2.0 * G, tc);
#endif
        c_scalar_mult_add_mat_f(&(s->f_U[mu]), &tmat, &tc, &(s->f_U[mu]));
      }
    }
  }

  // Third we have the scalar potential derivative contribution
  //   Udag_mu(x) 2B^2/N Tr[DmuUmu](x) Y(x)
  // where Y(x) = Tr[U_mu(x) Udag_mu(x)] / N - 1
  // Only compute if B is non-zero
  if (doB) {
    Real tr, twoBSqOvN = 2.0 * B * B / (Real)NCOL;

    FORALLSITES(i, s) {
      tc = trace_f(&DmuUmu[i]);
      for (mu = XUP; mu < NUMLINK; mu++) {
        tr = 1.0 / (Real)NCOL;
        tr *= realtrace_f(&(s->linkf[mu]), &(s->linkf[mu]));
        tr -= 1.0;
        CMULREAL(tc, twoBSqOvN * tr, tc2);

        adjoint_f(&(s->linkf[mu]), &tmat);
        c_scalar_mult_add_mat_f(&(s->f_U[mu]), &tmat, &tc2,
                                   &(s->f_U[mu]));
      }
    }
  }

  // Overall factor of C2 on all three potential d^2 contributions
  if (C2 - 1.0 > IMAG_TOL) {
    FORALLSITES(i, s) {
      FORALLDIR(mu)
        scalar_mult_matrix_f(&(s->f_U[mu]), C2, &(s->f_U[mu]));
    }
  }

  // Contribution from Fbar_{ab} F_{ab} term
  // Start first set of gathers (mu = 0 and nu = 1)
  for (mu = 0; mu < 2; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[2 + mu];
  }
  mat[0] = tempmat;
  mat[1] = tempmat2;

  tag0[0] = start_gather_site(F_OFFSET(linkf[1]), sizeof(matrix_f),
                              goffset[0], EVENANDODD, local_pt[0][0]);

  index = plaq_index[0][1];
  FORALLSITES(i, s)     // mu = 0 < nu = 1
    mult_an_f(&(Fmunu[index][i]), &(s->linkf[1]), &(mat[0][i]));
  tag1[0] = start_gather_field(mat[0], sizeof(matrix_f),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]);

  // Main loop
  FORALLDIR(mu) {
    FORALLDIR(nu) {
      if (mu == nu)
        continue;

      gather = (flip + 1) % 2;
      if (mu < NUMLINK - 1 || nu < NUMLINK - 2) { // Start next gathers
        if (nu == NUMLINK - 1) {
          a = mu + 1;
          b = 0;
        }
        else if (nu == mu - 1) {
          a = mu;
          b = nu + 2;
        }
        else {
          a = mu;
          b = nu + 1;
        }

        tag0[gather] = start_gather_site(F_OFFSET(linkf[b]),
                                         sizeof(matrix_f), goffset[a],
                                         EVENANDODD, local_pt[gather][0]);

        next = plaq_index[a][b];
        FORALLSITES(i, s) {
          if (a > b) {
            scalar_mult_matrix_f(&(Fmunu[next][i]), -1.0, &tmat);
            mult_an_f(&tmat, &(s->linkf[b]), &(mat[gather][i]));
          }
          else
            mult_an_f(&(Fmunu[next][i]), &(s->linkf[b]), &(mat[gather][i]));
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix_f),
                                          goffset[b] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu) {
          scalar_mult_matrix_f(&(Fmunu[index][i]), -1.0, &tmat);
          mult_na_f((matrix_f *)local_pt[flip][0][i], &tmat, &tmat2);
        }
        else
          mult_na_f((matrix_f *)local_pt[flip][0][i],
                        &(Fmunu[index][i]), &tmat2);

        sub_matrix_f(&tmat2, (matrix_f *)local_pt[flip][1][i], &tmat);
        scalar_mult_add_matrix_f(&(s->f_U[mu]), &tmat, 2.0, &(s->f_U[mu]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }

  // Factor of kappa = N / (2lambda) on both d^2 and F^2 terms
  FORALLSITES(i, s) {
    FORALLDIR(mu)
      scalar_mult_matrix_f(&(s->f_U[mu]), kappa, &(s->f_U[mu]));
  }

  // Only compute U(1) mass term if non-zero -- note factor of kappa
  if (bmass > IMAG_TOL) {
    Real tr, dmu = 2.0 * kappa * bmass * bmass / (Real)(NCOL * NCOL);
    FORALLSITES(i, s) {
      FORALLDIR(mu) {
        tr = realtrace_f(&(s->linkf[mu]), &(s->linkf[mu])) - (Real)NCOL;
        tr *= dmu;
        scalar_mult_add_adj_matrix_f(&(s->f_U[mu]), &(s->linkf[mu]),
                                         tr, &(s->f_U[mu]));
      }
    }
  }

  // Finally take adjoint and update the momentum
  // Subtract to reproduce -Adj(f_U)
  FORALLSITES(i, s) {
    FORALLDIR(mu)
      scalar_mult_sub_adj_matrix_f(&(s->mom[mu]), &(s->f_U[mu]), eps,
                                   &(s->mom[mu]));
  }

  // Compute average gauge force
  FORALLSITES(i, s) {
    FORALLDIR(mu)
      returnit += realtrace_f(&(s->f_U[mu]), &(s->f_U[mu]));
  }
  g_doublesum(&returnit);

  // Add in force from separate determinant term if kappa_u1 non-zero
  if (kappa_u1 > IMAG_TOL)
    returnit += det_force(eps);

  return (eps * sqrt(returnit) / volume);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Separate routines for each term in the fermion force
// All called by assemble_fermion_force below
// First Q-closed piece: chi_ab D_c chi_de epsilon_{abcde}
// Note factor of -1/2
void F1Q(vector *plaq_sol[NPLAQ], vector *plaq_psol[NPLAQ]) {
  register int i;
  register site *s;
  char **local_pt[2][4];
  int a, b, c, d, e, j, l, m, i_ab, i_de, gather, next, flip = 0;
  Real permm, tr;
  complex tc;
  msg_tag *tag0[2], *tag1[2], *tag2[2], *tag3[2];
  vector *vec0, *vec2, tvec, tvec2;
  matrix_f tmat;

  for (a = 0; a < 4; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[4 + a];
  }

  // Start first set of gathers
  // From setup_lamba.c, we see b > a and e > d
  a = FQ_lookup[0][0];
  b = FQ_lookup[0][1];
  c = FQ_lookup[0][2];
  d = FQ_lookup[0][3];
  e = FQ_lookup[0][4];
  i_ab = plaq_index[a][b];
  i_de = plaq_index[d][e];

  tag0[0] = start_gather_field(plaq_psol[i_de], sizeof(vector),
                               F1Q_d2[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(plaq_sol[i_ab], sizeof(vector),
                               goffset[c], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(plaq_psol[i_de], sizeof(vector),
                               goffset[c], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_field(plaq_sol[i_ab], sizeof(vector),
                               F1Q_d1[0], EVENANDODD, local_pt[0][3]);

  // Loop over lookup table
  for (j = 0; j < NTERMS; j++) {
    gather = (flip + 1) % 2;
    if (j < NTERMS - 1) {         // Start next set of gathers
      next = j + 1;
      a = FQ_lookup[next][0];
      b = FQ_lookup[next][1];
      c = FQ_lookup[next][2];
      d = FQ_lookup[next][3];
      e = FQ_lookup[next][4];
      i_ab = plaq_index[a][b];
      i_de = plaq_index[d][e];

      tag0[gather] = start_gather_field(plaq_psol[i_de], sizeof(vector),
                                        F1Q_d2[next], EVENANDODD,
                                        local_pt[gather][0]);
      tag1[gather] = start_gather_field(plaq_sol[i_ab], sizeof(vector),
                                        goffset[c], EVENANDODD,
                                        local_pt[gather][1]);
      tag2[gather] = start_gather_field(plaq_psol[i_de], sizeof(vector),
                                        goffset[c], EVENANDODD,
                                        local_pt[gather][2]);
      tag3[gather] = start_gather_field(plaq_sol[i_ab], sizeof(vector),
                                        F1Q_d1[next], EVENANDODD,
                                        local_pt[gather][3]);
    }

    // Do this set of computations while next set of gathers runs
    a = FQ_lookup[j][0];
    b = FQ_lookup[j][1];
    c = FQ_lookup[j][2];
    d = FQ_lookup[j][3];
    e = FQ_lookup[j][4];
    permm = perm[d][e][c][a][b];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);
    wait_gather(tag3[flip]);
    FORALLSITES(i, s) {
      clear_mat_f(&(tempmat[i]));
      vec0 = (vector *)(local_pt[flip][0][i]);
      vec2 = (vector *)(local_pt[flip][2][i]);

      tr = permm * (s->bc3[a][b][c]) * (s->bc1[c]);
      scalar_mult_vector((vector *)(local_pt[flip][1][i]), tr, &tvec);

      tr = -1.0 * permm * (s->bc2[OPP_LDIR(a)][OPP_LDIR(b)]) * (s->bc1[c]);
      scalar_mult_vector((vector *)(local_pt[flip][3][i]), tr, &tvec2);

      for (l = 0; l < DIMF; l++) {
        for (m = 0; m < DIMF; m++) {
          CMULJ_(vec0->c[l], tvec.c[m], tc);
          c_scalar_mult_add_mat_f(&(tempmat[i]), &(Lambda_prod[l][m]), &tc,
                                  &(tempmat[i]));

          CMULJ_(vec2->c[l], tvec2.c[m], tc);
          c_scalar_mult_add_mat_f(&(tempmat[i]), &(Lambda_prod[m][l]), &tc,
                                  &(tempmat[i]));
        }
      }
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    cleanup_gather(tag3[flip]);
    flip = (flip + 1) % 2;

    FORALLSITES(i, s) {
      adjoint_f(&(tempmat[i]), &tmat);
      scalar_mult_add_matrix_f(&(s->f_U[c]), &tmat, -0.5, &(s->f_U[c]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Second Q-closed piece
// Note factor of -1/2
void F2Q(vector *plaq_sol[NPLAQ], vector *plaq_psol[NPLAQ]) {
  register int i;
  register site *s;
  char **local_pt[2][4];
  int a, b, c, d, e, j, l, m, i_ab, i_de, gather, next, flip = 0;
  Real permm, tr;
  complex tc;
  msg_tag *tag0[2], *tag1[2], *tag2[2], *tag3[2];
  vector *vec0, *vec2, tvec, tvec2;
  matrix_f tmat;

  for (a = 0; a < 4; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[4 + a];
  }

  // Start first set of gathers
  // From setup_lamba.c, we see b > a and e > d
  a = FQ_lookup[0][0];
  b = FQ_lookup[0][1];
  c = FQ_lookup[0][2];
  d = FQ_lookup[0][3];
  e = FQ_lookup[0][4];
  i_ab = plaq_index[a][b];
  i_de = plaq_index[d][e];

  tag0[0] = start_gather_field(plaq_psol[i_ab], sizeof(vector),
                               F2Q_d1[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(plaq_sol[i_de], sizeof(vector),
                               goffset[c], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(plaq_psol[i_ab], sizeof(vector),
                               goffset[c], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_field(plaq_sol[i_de], sizeof(vector),
                               F2Q_d2[0], EVENANDODD, local_pt[0][3]);

  // Loop over lookup table
  for (j = 0; j < NTERMS; j++) {
    gather = (flip + 1) % 2;
    if (j < NTERMS - 1) {         // Start next set of gathers
      next = j + 1;
      a = FQ_lookup[next][0];
      b = FQ_lookup[next][1];
      c = FQ_lookup[next][2];
      d = FQ_lookup[next][3];
      e = FQ_lookup[next][4];
      i_ab = plaq_index[a][b];
      i_de = plaq_index[d][e];

      tag0[gather] = start_gather_field(plaq_psol[i_ab], sizeof(vector),
                                        F2Q_d1[next], EVENANDODD,
                                        local_pt[gather][0]);
      tag1[gather] = start_gather_field(plaq_sol[i_de], sizeof(vector),
                                        goffset[c], EVENANDODD,
                                        local_pt[gather][1]);
      tag2[gather] = start_gather_field(plaq_psol[i_ab], sizeof(vector),
                                        goffset[c], EVENANDODD,
                                        local_pt[gather][2]);
      tag3[gather] = start_gather_field(plaq_sol[i_de], sizeof(vector),
                                        F2Q_d2[next], EVENANDODD,
                                        local_pt[gather][3]);
    }

    // Do this set of computations while next set of gathers runs
    a = FQ_lookup[j][0];
    b = FQ_lookup[j][1];
    c = FQ_lookup[j][2];
    d = FQ_lookup[j][3];
    e = FQ_lookup[j][4];
    permm = perm[a][b][c][d][e];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);
    wait_gather(tag3[flip]);
    FORALLSITES(i, s) {
      clear_mat_f(&(tempmat[i]));
      vec0 = (vector *)(local_pt[flip][0][i]);
      vec2 = (vector *)(local_pt[flip][2][i]);

      tr = permm * (s->bc2[OPP_LDIR(a)][OPP_LDIR(b)]) * (s->bc1[c]);
      scalar_mult_vector((vector *)(local_pt[flip][1][i]), tr, &tvec);

      tr = -1.0 * permm * (s->bc3[a][b][c]) * (s->bc1[c]);
      scalar_mult_vector((vector *)(local_pt[flip][3][i]), tr, &tvec2);

      for (l = 0; l < DIMF; l++) {
        for (m = 0; m < DIMF; m++) {
          CMULJ_(vec0->c[l], tvec.c[m], tc);
          c_scalar_mult_add_mat_f(&(tempmat[i]), &(Lambda_prod[l][m]), &tc,
                                  &(tempmat[i]));

          CMULJ_(vec2->c[l], tvec2.c[m], tc);
          c_scalar_mult_add_mat_f(&(tempmat[i]), &(Lambda_prod[m][l]), &tc,
                                  &(tempmat[i]));
        }
      }
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    cleanup_gather(tag3[flip]);
    flip = (flip + 1) % 2;

    FORALLSITES(i, s) {
      adjoint_f(&(tempmat[i]), &tmat);
      scalar_mult_add_matrix_f(&(s->f_U[c]), &tmat, -0.5, &(s->f_U[c]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant contributions to the fermion force
// Use Uinv, Udag_inv, UpsiU, Tr_Uinv and tr_dest for temporary storage
// Also use tempdet and (if global det) tempZW for temporary storage
// The accumulator names refer to the corresponding derivatives
// Assume compute_plaqdet() has already been run
// If sign = 1 then take adjoint of eta
// If sign = -1 then take adjoint of psi
// A bit more code reuse may be possible
void detF(vector *eta, vector *psi[NUMLINK], int sign) {
  register int i;
  register site *s;
  int a, b, j;
  complex tc, tc2, tc3, Gc = cmplx(0.0, C2 * G * sqrt((Real)NCOL));
#ifdef LINEAR_DET
  CMULREAL(Gc, 0.5, Gc);                  // Since not squared
#else
  Real tr;
#endif
  msg_tag *mtag[8];
  matrix_f tmat, tmat2;

  // Check sign while giving Gc proper sign
  if (sign == 1) {    // Braces suppress compiler error
    CNEGATE(Gc, Gc);
  }
  else if (sign != -1) {
    node0_printf("Error: incorrect sign in detF: %d\n", sign);
    terminate(1);
  }

  // Set up and store some basic ingredients
  // Need all five directions for upcoming sums
  FORALLDIR(a) {
    // Save U_a(x)^{-1} in Uinv[a] and Udag_a(x)^{-1} in Udag_inv[a]
    // Save sum_j Tr[U_a(x)^{-1} Lambda^j] psi_a^j(x) in Tr_Uinv[a]
    // Save sum_j U_a(x)^{-1} Lambda^j psi_a^j(x) U_a(x)^{-1} in UpsiU[a]
    FORALLSITES(i, s) {
      invert(&(s->linkf[a]), &(Uinv[a][i]));
      adjoint_f(&(s->linkf[a]), &tmat);
      invert(&tmat, &(Udag_inv[a][i]));

      // Initialize accumulators for sum over j
      Tr_Uinv[a][i] = cmplx(0.0, 0.0);
      clear_mat_f(&(UpsiU[a][i]));
      for (j = 0; j < DIMF; j++) {
        mult_nn_f(&(Uinv[a][i]), &(Lambda[j]), &tmat);
        if (sign == 1)
          tc = psi[a][i].c[j];
        else
          CONJG(psi[a][i].c[j], tc);

        // Accumulate trace
        tc2 = trace_f(&tmat);
        CMUL(tc, tc2, tc3);
        CSUM(Tr_Uinv[a][i], tc3);

        // Accumulate product
        mult_nn_f(&tmat, &(Uinv[a][i]), &tmat2);
        c_scalar_mult_add_mat_f(&(UpsiU[a][i]), &tmat2, &tc, &(UpsiU[a][i]));
      }

      // tempdet holds either eta^D(x) plaqdet[a][b](x) (global)
      //                   or eta^D(x) |plaqdet[a][b](x)|^2 (local)
      if (sign == 1) {    // Braces suppress compiler error
        CONJG(eta[i].c[DIMF - 1], tc);
      }
      else
        tc = eta[i].c[DIMF - 1];

      for (b = a + 1; b < NUMLINK; b++) {
#ifdef LINEAR_DET
        CMUL(tc, plaqdet[a][b][i], tempdet[a][b][i]);
        CMUL(tc, plaqdet[b][a][i], tempdet[b][a][i]);
#else
        tr = cabs_sq(&(plaqdet[a][b][i]));
        CMULREAL(tc, tr, tempdet[a][b][i]);
        // Square is symmetric under a<-->b
        tempdet[b][a][i] = tempdet[a][b][i];
#endif
      }
    }
  }

  // Now we are ready to gather, accumulate and add to force
  // This is specialized in two big chunks, first global then local
#ifdef LINEAR_DET
  complex *plaq_term = malloc(sites_on_node * sizeof(*plaq_term));
  complex *inv_term = malloc(sites_on_node * sizeof(*inv_term));
  complex *adj_term = malloc(sites_on_node * sizeof(*adj_term));

  // Now we are ready to gather, accumulate and add to force
  // TODO: Could try to overlap these gathers, but that looks nasty...
  FORALLDIR(a) {
    // Initialize accumulators for sums over b
    FORALLSITES(i, s) {
      plaq_term[i] = cmplx(0.0, 0.0);
      inv_term[i] = cmplx(0.0, 0.0);
      adj_term[i] = cmplx(0.0, 0.0);
    }
    FORALLDIR(b) {
      if (a == b)
        continue;

      // Summary of gathers and shorthand:
      //   D[a][b](x) is eta^D(x) tempdet[a][b](x)
      //   T[a](x) is Tr[U_a(x)^{-1} psi_a(x)]
      // 0) T[b](x + a - b) in two steps
      // 1) D[a][b](x - b)
      // 2) D[b][a](x - b)
      // 3) T[a](x + b)
      // 4) T[a](x - b)
      // 5) T[b](x + a)
      // 6) T[b](x - b)
      mtag[0] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[a], EVENANDODD, gen_pt[0]);
      mtag[1] = start_gather_field(tempdet[a][b], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[1]);
      mtag[2] = start_gather_field(tempdet[b][a], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[2]);
      mtag[3] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                   goffset[b], EVENANDODD, gen_pt[3]);
      mtag[4] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[4]);
      mtag[5] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[a], EVENANDODD, gen_pt[5]);
      mtag[6] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[6]);

      // Step two of Tr_Uinv[b](x - b + a) gather, including BC
      // Use tr_dest for temporary storage
      wait_gather(mtag[0]);
      FORALLSITES(i, s)
        CMULREAL(*((complex *)(gen_pt[0][i])), s->bc1[a], tr_dest[i]);
      cleanup_gather(mtag[0]);
      mtag[0] = start_gather_field(tr_dest, sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[0]);

      // Now accumulate all three terms
      wait_gather(mtag[1]);         // 1) D[a][b](x - b)
      wait_gather(mtag[2]);         // 2) D[b][a](x - b)
      wait_gather(mtag[3]);         // 5) T[a](x + b)
      wait_gather(mtag[4]);         // 6) T[a](x - b)
      wait_gather(mtag[5]);         // 3) T[b](x + a)
      wait_gather(mtag[6]);         // 4) T[b](x - b)
      wait_gather(mtag[0]);         // 0) T[b](x + a - b)
      FORALLSITES(i, s) {
        // Accumulate plaq_term
        // D[b][a](x) {T[a](x) + T[b](x + a)}
        // gen_pt[5] is T[b](x + a)
        CMULREAL(*((complex *)(gen_pt[5][i])), s->bc1[a], tc);
        CADD(Tr_Uinv[a][i], tc, tc2);
        CMUL(tempdet[b][a][i], tc2, tc);
        CSUM(plaq_term[i], tc);

        // D[a][b](x - b) {T[a](x) + T[b](x - b)}
        CMULREAL(Tr_Uinv[a][i], s->bc1[OPP_LDIR(b)], tc);
        // gen_pt[6] is T[b](x - b)
        CADD(tc, *((complex *)(gen_pt[6][i])), tc2);
        // gen_pt[1] is D[a][b](x - b)
        CMUL(*((complex *)(gen_pt[1][i])), tc2, tc);
        CSUM(plaq_term[i], tc);

        // Accumulate adj_term
        // D[a][b](x) {T[b](x) + T[a](x + b)}
        // gen_pt[3] is T[a](x + b)
        CMULREAL(*((complex *)(gen_pt[3][i])), s->bc1[b], tc);
        CADD(Tr_Uinv[b][i], tc, tc2);
        CMUL(tempdet[a][b][i], tc2, tc);
        CSUM(adj_term[i], tc);

        // D[b][a](x - b) {T[a](x - b) + T[b](x + a - b) bc1[a](x - b)}
        // gen_pt[0] is T[b](x + a - b)
        // gen_pt[4] is T[a](x - b)
        CADD(*((complex *)(gen_pt[0][i])), *((complex *)(gen_pt[4][i])), tc);
        // gen_pt[2] is D[b][a](x - b)
        CMUL(*((complex *)(gen_pt[2][i])), tc, tc2);
        CSUM(adj_term[i], tc2);

        // Accumulate inv_term = sum_b D[b][a](x) + D[a][b](x - b)
        // gen_pt[1] is D[a][b](x - b)
        CMULREAL(*((complex *)(gen_pt[1][i])), s->bc1[OPP_LDIR(b)], tc);
        CSUM(inv_term[i], tc);
        CSUM(inv_term[i], tempdet[b][a][i]);
      }
      cleanup_gather(mtag[0]);
      cleanup_gather(mtag[1]);
      cleanup_gather(mtag[2]);
      cleanup_gather(mtag[3]);
      cleanup_gather(mtag[4]);
      cleanup_gather(mtag[5]);
      cleanup_gather(mtag[6]);
    }

    // Now add to force
    // Include complex coupling Gc before taking adjoint
    FORALLSITES(i, s) {
      // Start with plaq_term hitting U_a(x)^{-1}
      CMUL(plaq_term[i], Gc, tc);
      c_scalar_mult_add_mat_f(&(s->f_U[a]), &(Uinv[a][i]), &tc, &(s->f_U[a]));

      // Add adj_term hitting Udag_a(x)^{-1} followed by adjoint
      CMUL(adj_term[i], Gc, tc);
      c_scalar_mult_add_adj_mat_f(&(s->f_U[a]), &(Udag_inv[a][i]), &tc,
                                  &(s->f_U[a]));

      // Finally subtract inv_term hitting U_a(x)^{-1} psi_a(x) U_a(x)^{-1}
      CMUL(inv_term[i], Gc, tc);
      c_scalar_mult_sub_mat_f(&(s->f_U[a]), &(UpsiU[a][i]), &tc, &(s->f_U[a]));
    }
  }
  free(plaq_term);
  free(inv_term);
  free(adj_term);
#else     // Local case
  complex *dZdU = malloc(sites_on_node * sizeof(*dZdU));
  complex *dWdU = malloc(sites_on_node * sizeof(*dWdU));
  complex *dZdUdag = malloc(sites_on_node * sizeof(*dZdUdag));
  complex *dWdUdag = malloc(sites_on_node * sizeof(*dWdUdag));
  complex *dTdU = malloc(sites_on_node * sizeof(*dTdU));

  // Set up and store one more ingredient
  for (a = XUP; a < NUMLINK; a++) {
    FORALLSITES(i, s) {
      // Save eta^D(x) ZWstar[a][b](x) in tempZW[a][b](x)
      if (sign == 1) {    // Braces suppress compiler error
        CONJG(eta[i].c[DIMF - 1], tc);
      }
      else
        tc = eta[i].c[DIMF - 1];

      for (b = a + 1; b < NUMLINK; b++) {
        CMUL(tc, ZWstar[a][b][i], tempZW[a][b][i]);
        CMUL(tc, ZWstar[b][a][i], tempZW[b][a][i]);
      }
    }
  }

  // Now we are ready to gather, accumulate and add to force
  // TODO: Could try to overlap these gathers, but that looks nasty...
  FORALLDIR(a) {
  for (a = XUP; a < NUMLINK; a++) {
    // Initialize accumulators for sums over b
    FORALLSITES(i, s) {
      dZdU[i] = cmplx(0.0, 0.0);
      dWdU[i] = cmplx(0.0, 0.0);
      dZdUdag[i] = cmplx(0.0, 0.0);
      dWdUdag[i] = cmplx(0.0, 0.0);
      dTdU[i] = cmplx(0.0, 0.0);
    }
    FORALLDIR(b) {
      if (a == b)
        continue;

      // Summary of gathers and shorthand:
      //   ZSq[a][b](x) is eta^{D*}(x) |tempdet[a][b](x)|^2
      //   ZW[a][b](x) is eta^{D*}(x)tempdet[a][b](x)[tempdet[a][b](x)-1]^*
      //   T[a](x) is Tr[U_a(x)^{-1} psi_a(x)]
      // 0) T[b](x - b + a) in two steps
      // 1) ZSq[a][b](x - b)
      // 2) ZW[a][b](x - b)
      // 3) ZW[b][a](x - b)
      // 4) T[a](x + b)
      // 5) T[a](x - b)
      // 6) T[b](x + a)
      // 7) T[b](x - b)
      mtag[0] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[a], EVENANDODD, gen_pt[0]);
      mtag[1] = start_gather_field(tempdet[a][b], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[1]);
      mtag[2] = start_gather_field(tempZW[a][b], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[2]);
      mtag[3] = start_gather_field(tempZW[b][a], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[3]);
      mtag[4] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                   goffset[b], EVENANDODD, gen_pt[4]);
      mtag[5] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[5]);
      mtag[6] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[a], EVENANDODD, gen_pt[6]);
      mtag[7] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[7]);

      // Step two of Tr_Uinv[b](x - b + a) gather, including BC
      // Use tr_dest for temporary storage
      wait_gather(mtag[0]);
      FORALLSITES(i, s) {
        tc = *((complex *)(gen_pt[0][i]));
        CMULREAL(tc, s->bc1[a], tr_dest[i]);
      }
      cleanup_gather(mtag[0]);
      mtag[0] = start_gather_field(tr_dest, sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[0]);

      // Now accumulate all five terms
      wait_gather(mtag[1]);       // 1) ZSq[a][b](x - b)
      wait_gather(mtag[2]);       // 2) ZW[a][b](x - b)
      wait_gather(mtag[3]);       // 3) ZW[b][a](x - b)
      wait_gather(mtag[4]);       // 4) T[a](x + b)
      wait_gather(mtag[5]);       // 5) T[a](x - b)
      wait_gather(mtag[6]);       // 6) T[b](x + a)
      wait_gather(mtag[7]);       // 7) T[b](x - b)
      wait_gather(mtag[0]);       // 0) T[b](x - b + a)
      FORALLSITES(i, s) {
        // dZdU and dWdUdag have same sums of traces
        // hit by ZW and ZSq, respectively
        // Z(x) {T[a](x) + BC[a](x) T[b](x + a)}
        tc = *((complex *)(gen_pt[6][i]));    // T[b](x + a)
        CMULREAL(tc, s->bc1[a], tc);
        CADD(Tr_Uinv[a][i], tc, tc2);
        CMUL(tempZW[b][a][i], tc2, tc);       // dZdU
        CSUM(dZdU[i], tc);
        CMUL(tempdet[b][a][i], tc2, tc);      // dWdUdag
        CSUM(dWdUdag[i], tc);

        // Z(x - b) {T[b](x - b) + BC[-b](x) T[a](x)}
        tc = *((complex *)(gen_pt[7][i]));    // T[b](x - b)
        CMULREAL(Tr_Uinv[a][i], s->bc1[OPP_LDIR(b)], tc2);
        CADD(tc, tc2, tc3);
        tc = *((complex *)(gen_pt[2][i]));    // ZW[a][b](x - b)
        CMUL(tc, tc3, tc2);
        CSUM(dZdU[i], tc2);
        tc = *((complex *)(gen_pt[1][i]));    // ZSq[a][b](x - b)
        CMUL(tc, tc3, tc2);
        CSUM(dWdUdag[i], tc2);

        // dWdU and dZdUdag have same sums of traces
        // hit by ZSq and ZW, respectively
        // Z(x) {T[b](x) + BC[b](x) T[a](x + b)}
        tc = *((complex *)(gen_pt[4][i]));    // T[a](x + b)
        CMULREAL(tc, s->bc1[b], tc);
        CADD(Tr_Uinv[b][i], tc, tc2);
        CMUL(tempdet[a][b][i], tc2, tc);      // dWdU
        CSUM(dWdU[i], tc);
        CMUL(tempZW[a][b][i], tc2, tc);       // dZdUdag
        CSUM(dZdUdag[i], tc);

        // Z(x - b) {T[a](x - b) + BC[a](x - b) T[b](x - b + a)}
        tc = *((complex *)(gen_pt[0][i]));    // T[b](x - b + a)
        tc2 = *((complex *)(gen_pt[5][i]));   // T[a](x - b)
        CADD(tc2, tc, tc3);
        tc = *((complex *)(gen_pt[1][i]));    // ZSq[a][b](x - b)
        CMUL(tc, tc3, tc2);
        CSUM(dWdU[i], tc2);
        tc = *((complex *)(gen_pt[3][i]));    // ZW[b][a](x - b)
        CMUL(tc, tc3, tc2);
        CSUM(dZdUdag[i], tc2);

        // Finally dTdU accumulates ZW[b][a](x) + BC[-b](x) ZW[a][b](x - b)
        tc = *((complex *)(gen_pt[2][i]));    // ZW[a][b](x - b)
        CMULREAL(tc, s->bc1[OPP_LDIR(b)], tc);
        CADD(tempZW[b][a][i], tc, tc2);
        CSUM(dTdU[i], tc2);
      }
      cleanup_gather(mtag[0]);
      cleanup_gather(mtag[1]);
      cleanup_gather(mtag[2]);
      cleanup_gather(mtag[3]);
      cleanup_gather(mtag[4]);
      cleanup_gather(mtag[5]);
      cleanup_gather(mtag[6]);
      cleanup_gather(mtag[7]);
    }

    // Now add to force
    // Include complex coupling Gc before taking adjoint
    FORALLSITES(i, s) {
      // Start with dZdU and dWdU hitting U_a(x)^{-1}
      CADD(dZdU[i], dWdU[i], tc);
      CMUL(tc, Gc, tc2);
      c_scalar_mult_add_mat_f(&(s->f_U[a]), &(Uinv[a][i]), &tc2, &(s->f_U[a]));

      // Add dZdUdag and dWdUdag hitting Udag_a(x)^{-1} followed by adjoint
      CADD(dZdUdag[i], dWdUdag[i], tc);
      CMUL(tc, Gc, tc2);
      c_scalar_mult_add_adj_mat_f(&(s->f_U[a]), &(Udag_inv[a][i]), &tc2,
                                  &(s->f_U[a]));

      // Finally subtract dTdU hitting U_a(x)^{-1} psi_a(x) U_a(x)^{-1}
      CMUL(dTdU[i], Gc, tc);
      c_scalar_mult_sub_mat_f(&(s->f_U[a]), &(UpsiU[a][i]), &tc, &(s->f_U[a]));
    }
  }
  free(dZdU);
  free(dWdU);
  free(dZdUdag);
  free(dWdUdag);
  free(dTdU);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar potential contributions to the fermion force
// Use tempmat and Tr_Uinv for temporary storage
// If sign = 1 then take adjoint of eta
// If sign = -1 then take adjoint of psi
void pot_force(vector *eta, vector *psi[NUMLINK], int sign) {
  register int i;
  register site *s;
  int a, j;
  Real tr;
  complex tc, tc2, tc3, Bc = cmplx(0.0, C2 * B * B / sqrt((Real)NCOL));
  matrix_f tmat;

  // Check sign while giving Bc proper sign
  if (sign == 1) {    // Braces suppress compiler error
    CNEGATE(Bc, Bc);
  }
  else if (sign != -1) {
    node0_printf("Error: incorrect sign in detF: %d\n", sign);
    terminate(1);
  }

  for (a = XUP; a < NUMLINK; a++) {
    // Save sum_j psi_a^j(x) Tr[Lambda^j Udag_a(x)] in tr_dest
    // Save sum_j psi_a^j(x) Lambda^j in tempmat
    FORALLSITES(i, s) {
      // Initialize accumulators for sum over j
      clear_mat_f(&(tempmat[i]));
      tr_dest[i] = cmplx(0.0, 0.0);
      for (j = 0; j < DIMF; j++) {
        mult_na_f(&(Lambda[j]), &(s->linkf[a]), &tmat);
        if (sign == 1)
          tc = psi[a][i].c[j];
        else
          CONJG(psi[a][i].c[j], tc);

        // Accumulate trace
        tc2 = trace_f(&tmat);
        CMUL(tc, tc2, tc3);
        CSUM(tr_dest[i], tc3);

        // Accumulate psi itself
        c_scalar_mult_add_mat_f(&(tempmat[i]), &(Lambda[j]), &tc,
                                &(tempmat[i]));
      }

      // Hit both with eta^{D*} and divide trace by N
      if (sign == 1) {    // Braces suppress compiler error
        CONJG(eta[i].c[DIMF - 1], tc);
      }
      else
        tc = eta[i].c[DIMF - 1];

      CMUL(tc, tr_dest[i], tc2);
      CDIVREAL(tc2, (Real)NCOL, tr_dest[i]);
      mat_copy_f(&(tempmat[i]), &tmat);
      c_scalar_mult_mat_f(&tmat, &tc, &(tempmat[i]));

      // Compute Y(x) = Tr[U_a(x) Udag_a(x)] / N - 1
      tr = 1.0 / (Real)NCOL;
      tr *= realtrace_f(&(s->linkf[a]), &(s->linkf[a]));
      tr -= 1.0;

      // We're already ready to add to force
      // Start with eta Tr / N hitting Udag_a(x)
      CMUL(Bc, tr_dest[i], tc);
      c_scalar_mult_add_mat_adj_f(&(s->f_U[a]), &(s->linkf[a]), &tc,
                                     &(s->f_U[a]));

      // Add eta Tr / N hitting U_a(x) and eta Y hitting psi_a(x)
      // and take the adjoint of the sum
      c_scalar_mult_mat_f(&(s->linkf[a]), &(tr_dest[i]), &tmat);
      scalar_mult_add_matrix_f(&tmat, &(tempmat[i]), tr, &tmat);
      c_scalar_mult_add_adj_mat_f(&(s->f_U[a]), &tmat, &Bc, &(s->f_U[a]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Assemble fermion contributions to gauge link force,
//   f_U = Adj(Ms).D_U M(U, Ub).s - Adj[Adj(Ms).D_Ub M(U, Ub).s]
// "s" is sol while "Ms" is psol
// Copy these into persistent vectors for easier gathering
// Use tempmat, Uinv, Tr_Uinv, tr_dest and Ddet[012] for temporary storage
void assemble_fermion_force(Twist_Fermion *sol, Twist_Fermion *psol) {
  register int i;
  register site *s;
  char **local_pt[2][2];
  int mu, nu, a, b, c, d, gather, flip = 0, index, next;
  complex tc;
  msg_tag *mtag[NUMLINK], *tag0[2], *tag1[2];
  vector tvec;
  matrix_f *mat[2];

  for (mu = 0; mu < 2; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[2 + mu];
  }
  mat[0] = tempmat;
  mat[1] = tempmat2;

  // For gathering it is convenient to copy the input Twist_Fermions
  // into persistent site, link and plaquette vectors
  // We can use "src" and "dest" vectors for this storage,
  // though we want to call them "sol" and "psol" for clarity
  vector *site_sol = site_src, *site_psol = site_dest;
  vector *link_sol[NUMLINK], *link_psol[NUMLINK];
  vector *plaq_sol[NPLAQ], *plaq_psol[NPLAQ];
  FORALLDIR(mu) {
    link_sol[mu] = link_src[mu];
    link_psol[mu] = link_dest[mu];
  }
  for (mu = 0; mu < NPLAQ; mu++) {
    plaq_sol[mu] = plaq_src[mu];
    plaq_psol[mu] = plaq_dest[mu];
  }
  FORALLSITES(i, s) {
    vec_copy(&(sol[i].Fsite), &(site_sol[i]));
    vec_copy(&(psol[i].Fsite), &(site_psol[i]));
    FORALLDIR(mu) {
      vec_copy(&(sol[i].Flink[mu]), &(link_sol[mu][i]));
      vec_copy(&(psol[i].Flink[mu]), &(link_psol[mu][i]));
    }
    for (mu = 0; mu < NPLAQ; mu++) {
      vec_copy(&(sol[i].Fplaq[mu]), &(plaq_sol[mu][i]));
      vec_copy(&(psol[i].Fplaq[mu]), &(plaq_psol[mu][i]));
    }
  }

#ifdef SV
  // Accumulate both terms in Uinv[mu], use to initialize f_U[mu]
  // First calculate DUbar on eta Dbar_mu psi_mu (LtoS)
  mtag[0] = start_gather_field(site_psol, sizeof(vector),
                               goffset[0], EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(site_psol, sizeof(vector),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);
    }
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      clear_mat_f(&(Uinv[mu][i]));
      for (a = 0; a < DIMF; a++) {
        for (b = 0; b < DIMF; b++) {
          CMULJ_((site_psol[i]).c[a], (link_sol[mu][i]).c[b], tc);
          c_scalar_mult_sub_mat_f(&(Uinv[mu][i]), &(Lambda_prod[a][b]), &tc,
                                  &(Uinv[mu][i]));

          // gen_pt[mu] is site_psol(x + mu)
          CMULJ_(((vector *)(gen_pt[mu][i]))->c[a],
                 (link_sol[mu][i]).c[b], tc);
          CMULREAL(tc, s->bc1[mu], tc);
          c_scalar_mult_add_mat_f(&(Uinv[mu][i]), &(Lambda_prod[b][a]), &tc,
                                  &(Uinv[mu][i]));
        }
      }
    }
    cleanup_gather(mtag[mu]);
  }

  // 2nd term, DUbar on psi_mu Dbar_mu eta (StoL)
  mtag[0] = start_gather_field(site_sol, sizeof(vector),
                               goffset[0], EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(site_sol, sizeof(vector),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);
    }
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_vector((vector *)(gen_pt[mu][i]), s->bc1[mu], &tvec);
      for (a = 0; a < DIMF; a++) {
        for (b = 0; b < DIMF; b++) {
          CMULJ_((link_psol[mu][i]).c[a], tvec.c[b], tc);
          c_scalar_mult_sub_mat_f(&(Uinv[mu][i]), &(Lambda_prod[a][b]), &tc,
                                  &(Uinv[mu][i]));

          CMULJ_((link_psol[mu][i]).c[a], (site_sol[i]).c[b], tc);
          c_scalar_mult_add_mat_f(&(Uinv[mu][i]), &(Lambda_prod[b][a]), &tc,
                                  &(Uinv[mu][i]));
        }
      }
      // Initialize the force collectors
      scalar_mult_adj_matrix_f(&(Uinv[mu][i]), 0.5, &(s->f_U[mu]));
    }
    cleanup_gather(mtag[mu]);
  }
#endif
#ifdef VP
  // Now calculate DU on chi_{munu} D_mu(U) psi_nu
  // Start first set of gathers (mu = 0 and nu = 1)
  tag0[0] = start_gather_field(link_sol[1], sizeof(vector),
                               goffset[0], EVENANDODD, local_pt[0][0]);
  index = plaq_index[0][1];
  FORALLSITES(i, s) {   // mu = 0 < nu = 1
    clear_mat_f(&(mat[0][i]));
    for (a = 0; a < DIMF; a++) {
      for (b = 0; b < DIMF; b++) {
        CMULJ_((plaq_psol[index][i]).c[a], (link_sol[1][i]).c[b], tc);
        c_scalar_mult_add_mat_f(&(mat[0][i]), &(Lambda_prod[a][b]), &tc,
                                &(mat[0][i]));
      }
    }
  }
  tag1[0] = start_gather_field(mat[0], sizeof(matrix_f),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]);

  FORALLDIR(mu) {
    FORALLDIR(nu) {
      if (mu == nu)
        continue;

      gather = (flip + 1) % 2;
      if (mu < NUMLINK - 1 || nu < NUMLINK - 2) {   // Start next gathers
        if (nu == NUMLINK - 1) {
          c = mu + 1;
          d = 0;
        }
        else if (nu == mu - 1) {
          c = mu;
          d = nu + 2;
        }
        else {
          c = mu;
          d = nu + 1;
        }
        tag0[gather] = start_gather_field(link_sol[d], sizeof(vector),
                                          goffset[c], EVENANDODD,
                                          local_pt[gather][0]);
        next = plaq_index[c][d];
        FORALLSITES(i, s) {
          clear_mat_f(&(mat[gather][i]));
          if (c > d) {    // plaq_psol is anti-symmetric under c <--> d
            scalar_mult_vector(&(plaq_psol[next][i]), -1.0, &tvec);
          }                 // Suppress compiler error
          else
            vec_copy(&(plaq_psol[next][i]), &tvec);
          for (a = 0; a < DIMF; a++) {
            for (b = 0; b < DIMF; b++) {
              CMULJ_(tvec.c[a], (link_sol[d][i]).c[b], tc);
              c_scalar_mult_add_mat_f(&(mat[gather][i]), &(Lambda_prod[a][b]), &tc,
                                      &(mat[gather][i]));
            }
          }
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix_f),
                                          goffset[d] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu) {    // plaq_psol is anti-symmetric under mu <--> nu
          scalar_mult_vector((vector *)(local_pt[flip][0][i]),
                             s->bc1[mu], &tvec);
        }                 // Suppress compiler error
        else
          scalar_mult_vector((vector *)(local_pt[flip][0][i]),
                             -1.0 * s->bc1[mu], &tvec);
        for (a = 0; a < DIMF; a++) {
          for (b = 0; b < DIMF; b++) {
            CMULJ_((plaq_psol[index][i]).c[a], tvec.c[b], tc);
            c_scalar_mult_add_mat_f(&(s->f_U[mu]), &(Lambda_prod[b][a]), &tc,
                                    &(s->f_U[mu]));
          }
        }
        add_matrix_f(&(s->f_U[mu]), (matrix_f *)(local_pt[flip][1][i]),
                     &(s->f_U[mu]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }

  // 2nd term
  // Start first set of gathers (mu = 0 and nu = 1)
  tag0[0] = start_gather_field(link_psol[1], sizeof(vector),
                               goffset[0], EVENANDODD, local_pt[0][0]);

  index = plaq_index[0][1];
  FORALLSITES(i, s) {   // mu = 0 < nu = 1
    clear_mat_f(&(mat[0][i]));
    for (a = 0; a < DIMF; a++) {
      for (b = 0; b < DIMF; b++) {
        CMULJ_((link_psol[1][i]).c[a], (plaq_sol[index][i]).c[b], tc);
        c_scalar_mult_add_mat_f(&(mat[0][i]), &(Lambda_prod[b][a]), &tc,
                                &(mat[0][i]));
      }
    }
  }
  tag1[0] = start_gather_field(mat[0], sizeof(matrix_f),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]);

  FORALLDIR(mu) {
    FORALLDIR(nu) {
      if (mu == nu)
        continue;

      gather = (flip + 1) % 2;
      if (mu < NUMLINK - 1 || nu < NUMLINK - 2) {   // Start next gathers
        if (nu == NUMLINK - 1) {
          c = mu + 1;
          d = 0;
        }
        else if (nu == mu - 1) {
          c = mu;
          d = nu + 2;
        }
        else {
          c = mu;
          d = nu + 1;
        }
        tag0[gather] = start_gather_field(link_psol[d], sizeof(vector),
                                          goffset[c], EVENANDODD,
                                          local_pt[gather][0]);

        next = plaq_index[c][d];
        FORALLSITES(i, s) {
          clear_mat_f(&(mat[gather][i]));
          if (c > d) {      // plaq_sol is anti-symmetric under c <--> d
            scalar_mult_vector(&(plaq_sol[next][i]), -1.0, &tvec);
          }                 // Suppress compiler error
          else
            vec_copy(&(plaq_sol[next][i]), &tvec);
          for (a = 0; a < DIMF; a++) {
            for (b = 0; b < DIMF; b++) {
              CMULJ_((link_psol[d][i]).c[a], tvec.c[b], tc);
              c_scalar_mult_add_mat_f(&(mat[gather][i]), &(Lambda_prod[b][a]), &tc,
                                      &(mat[gather][i]));
            }
          }
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix_f),
                                          goffset[d] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu) {    // plaq_sol is anti-symmetric under mu <--> nu
          scalar_mult_vector(&(plaq_sol[index][i]), -1.0 * s->bc1[mu],
                                 &tvec);
        }                 // Suppress compiler error
        else
          scalar_mult_vector(&(plaq_sol[index][i]), s->bc1[mu], &tvec);
        for (a = 0; a < DIMF; a++) {
          for (b = 0; b < DIMF; b++) {
            CMULJ_(((vector *)(local_pt[flip][0][i]))->c[a], tvec.c[b], tc);
            c_scalar_mult_add_mat_f(&(s->f_U[mu]), &(Lambda_prod[a][b]), &tc,
                                    &(s->f_U[mu]));
          }
        }
        sub_matrix_f(&(s->f_U[mu]), (matrix_f *)(local_pt[flip][1][i]),
                     &(s->f_U[mu]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }
#endif

  // Plaquette determinant contributions if G is non-zero
  if (doG) {
    // First connect link_sol with site_psol[DIMF - 1]^dag (LtoS)
    detF(site_psol, link_sol, PLUS);

    // Second connect site_sol[DIMF - 1] with link_psol^dag (StoL)
    detF(site_sol, link_psol, MINUS);
  }

  // Scalar potential contributions if B is non-zero
  // Use tempmat and Tr_Uinv for temporary storage
  if (doB) {
    // First connect link_sol with site_psol[DIMF - 1]^dag (LtoS)
    pot_force(site_psol, link_sol, PLUS);

    // Second connect site_sol[DIMF - 1] with link_psol^dag (StoL)
    pot_force(site_sol, link_psol, MINUS);
  }

#ifdef QCLOSED
  if (NUMLINK != 5) {
    node0_printf("ERROR: NUMLINK IS %d != 5\n", NUMLINK);
    terminate(1);
  }
  F1Q(plaq_sol, plaq_psol);
  F2Q(plaq_sol, plaq_psol);
#endif
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Update the momenta with the fermion force
// Assume that the multiCG has been run (updating the adjoint links),
// with the solution vectors in sol[j]
// Accumulate f_U for each pole into fullforce, add to momenta
// Use fullforce-->Fmunu and tempTF for temporary storage
double fermion_force(Real eps, Twist_Fermion *src, Twist_Fermion **sol) {
  register int i;
  register site *s;
  int mu, n;
  double returnit = 0.0;
  matrix_f **fullforce = malloc(NUMLINK * sizeof(*fullforce));

#ifdef FORCE_DEBUG
  int kick, ii, jj, iters = 0;
  Real final_rsq;
  double individ_force, old_action, new_action = 0.0;
  matrix_f tmat, tprint, tprint2;
  clear_mat_f(&tprint);
  clear_mat_f(&tmat);
#endif

  FORALLDIR(mu)
    fullforce[mu] = Fmunu[mu];    // Use Fmunu for temporary storage

  // Initialize fullforce[mu]
  fermion_op(sol[0], tempTF, PLUS);
  FORALLSITES(i, s)
    scalar_mult_TF(&(tempTF[i]), amp4[0], &(tempTF[i]));
  assemble_fermion_force(sol[0], tempTF);
  FORALLDIR(mu) {
    FORALLSITES(i, s)
      adjoint_f(&(s->f_U[mu]), &(fullforce[mu][i]));
  }
  for (n = 1; n < Norder; n++) {
    fermion_op(sol[n], tempTF, PLUS);
    // Makes sense to multiply here by amp4[n]...
    FORALLSITES(i, s)
      scalar_mult_TF(&(tempTF[i]), amp4[n], &(tempTF[i]));

    assemble_fermion_force(sol[n], tempTF);
#ifdef FORCE_DEBUG
    individ_force = 0.0;
#endif
    FORALLDIR(mu) {
      FORALLSITES(i, s) {
        // Take adjoint but don't negate yet...
        add_adj_matrix_f(&(fullforce[mu][i]), &(s->f_U[mu]),
                             &(fullforce[mu][i]));
#ifdef FORCE_DEBUG
//      if (s->x == 0 && s->y == 0 && s->z == 0 && s->t == 0 && mu == 3) {
//        printf("Fermion force mu=%d on site (%d, %d, %d, %d)\n",
//               mu, s->x, s->y, s->z ,s->t);
//        dumpmat_f(&(s->f_U[mu]));
//      }
        // Compute average gauge force
        individ_force += realtrace_f(&(s->f_U[mu]), &(s->f_U[mu]));
#endif
      }
    }
#ifdef FORCE_DEBUG
    g_doublesum(&individ_force);
    node0_printf("Individ_force %d %.4g\n",
                 n, eps * sqrt(individ_force) / volume);

    // Check that force syncs with fermion action
    // congrad_multi_field calls fermion_rep()
    old_action = d_fermion_action(src, sol);
    iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq);
    new_action = d_fermion_action(src, sol);
    node0_printf("EXITING  %.4g\n", new_action - old_action);
    if (fabs(new_action - old_action) > 1e-3)
      terminate(1);                             // Don't go further for now

#if 0
    // Do a scan of the fermion action
    for (mu = XUP; mu < NUMLINK; mu++) {
      FORALLSITES(i, s) {
        node0_printf("mu=%d on site (%d, %d, %d, %d)\n",
                     mu, s->x, s->y, s->z, s->t);
        tmat = s->linkf[mu];
        dumpmat_f(&(s->f_U[mu]));

        for (ii = 0; ii < NCOL; ii++) {
          for (jj = 0; jj < NCOL; jj++) {
            for (kick = -1; kick <= 1; kick += 2) {
              s->linkf[mu] = tmat;
              s->linkf[mu].e[ii][jj].real += 0.001 * (Real)kick;

              iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq);
              if (kick == -1)
                new_action -= d_fermion_action(src, sol);
              if (kick == 1) {
                new_action += d_fermion_action(src, sol);
                tprint.e[ii][jj].real = -250.0 * new_action;
              }
            }

            for (kick = -1; kick <= 1; kick += 2) {
              s->linkf[mu] = tmat;
              s->linkf[mu].e[ii][jj].imag += 0.001 * (Real)kick;

              iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq);
              if (kick == -1)
                new_action -= d_fermion_action(src, sol);
              if (kick == 1) {
                new_action += d_fermion_action(src, sol);
                node0_printf("XXXG%d%dI %.4g %.4g\n",
                             ii, jj, 0.001 * (Real)kick, 500 * new_action);
                tprint.e[ii][jj].imag = -250 * new_action;
              }
            }
          }
        }
        sub_matrix_f(&tprint, &(s->f_U[mu]), &tprint2);
        node0_printf("mu=%d on site (%d, %d, %d, %d): %.4g\n",
                     mu, s->x, s->y, s->z, s->t,
                     realtrace_f(&tprint2, &tprint2));
        dumpmat_f(&tprint);
        s->linkf[mu] = tmat;

        iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq);
      }
    }   // End scan of the fermion action
#endif
#endif
  }

  // Update the momentum from the fermion force -- sum or eps
  // Opposite sign as to gauge force,
  // because dS_G / dU = 2F_g while ds_F / dU = -2F_f
  // Move negation here as well, though adjoint remains above
  FORALLSITES(i, s) {
    FORALLDIR(mu) {
      scalar_mult_sub_matrix_f(&(s->mom[mu]), &(fullforce[mu][i]), eps,
                                   &(s->mom[mu]));
      returnit += realtrace_f(&(fullforce[mu][i]), &(fullforce[mu][i]));
    }
  }
  g_doublesum(&returnit);

  free(fullforce);

  // Reset Fmunu
  compute_Fmunu();
  return (eps * sqrt(returnit) / volume);
}
// -----------------------------------------------------------------
