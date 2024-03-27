// -----------------------------------------------------------------
// Update the momentum matrices
// Uncomment to print out debugging messages
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
  double returnit = 0.0, tr;
  complex tc;
  matrix tmat, tmat2, *mat[2];
  msg_tag *tag[NUMLINK], *tag0[2], *tag1[2];

  // All contributions from d^2 term need a factor of C2
  // First we have the finite difference operator derivative times DmuUmu
  // Ubar_a(x) DmuUmu(x) - DmuUmu(x + a) Ubar_a(x)
  tag[0] = start_gather_field(DmuUmu, sizeof(matrix),
                              goffset[0], EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_field(DmuUmu, sizeof(matrix),
                                       goffset[mu + 1], EVENANDODD,
                                       gen_pt[mu + 1]);

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      mult_an(&(s->link[mu]), &(DmuUmu[i]), &(s->f_U[mu]));   // Initialize
      mult_na_dif((matrix *)(gen_pt[mu][i]), &(s->link[mu]), &(s->f_U[mu]));
    }
    cleanup_gather(tag[mu]);
  }

  // Next we have the plaquette determinant derivative contribution
  //   U_mu^{-1}(x) 2G sum_nu {D[nu][mu](x) + D[mu][nu](x-nu)}
  // D is Tr[DmuUmu] plaqdet[mu][nu], saved in tempdet[mu][nu]
  // Only compute if G is non-zero
  // Use tr_dest for temporary storage
  if (doG) {
    FORALLSITES(i, s) {
      tc = trace(&DmuUmu[i]);
      FORALLDIR(mu) {
        for (nu = mu + 1; nu < NUMLINK; nu++) {
          CMUL(tc, plaqdet[mu][nu][i], tempdet[mu][nu][i]);
          CMUL(tc, plaqdet[nu][mu][i], tempdet[nu][mu][i]);
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
        CMULREAL(tr_dest[i], G, tc);
        c_scalar_mult_sum_mat(&(Uinv[mu][i]), &tc, &(s->f_U[mu]));
      }
    }
  }

  // Overall factor of C2 on all d^2 contributions
  if (C2 - 1.0 > IMAG_TOL) {
    FORALLSITES(i, s) {
      FORALLDIR(mu)
        scalar_mult_matrix(&(s->f_U[mu]), C2, &(s->f_U[mu]));
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

  tag0[0] = start_gather_site(F_OFFSET(link[1]), sizeof(matrix),
                              goffset[0], EVENANDODD, local_pt[0][0]);

  index = plaq_index[0][1];
  FORALLSITES(i, s)     // mu = 0 < nu = 1
    mult_an(&(Fmunu[index][i]), &(s->link[1]), &(mat[0][i]));
  tag1[0] = start_gather_field(mat[0], sizeof(matrix),
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

        tag0[gather] = start_gather_site(F_OFFSET(link[b]),
                                         sizeof(matrix), goffset[a],
                                         EVENANDODD, local_pt[gather][0]);

        next = plaq_index[a][b];
        FORALLSITES(i, s) {
          if (a > b) {
            scalar_mult_matrix(&(Fmunu[next][i]), -1.0, &tmat);
            mult_an(&tmat, &(s->link[b]), &(mat[gather][i]));
          }
          else
            mult_an(&(Fmunu[next][i]), &(s->link[b]), &(mat[gather][i]));
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix),
                                          goffset[b] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu) {
          scalar_mult_matrix(&(Fmunu[index][i]), -1.0, &tmat);
          mult_na((matrix *)local_pt[flip][0][i], &tmat, &tmat2);
        }
        else
          mult_na((matrix *)local_pt[flip][0][i], &(Fmunu[index][i]), &tmat2);

        sub_matrix(&tmat2, (matrix *)local_pt[flip][1][i], &tmat);
        scalar_mult_sum_matrix(&tmat, 2.0, &(s->f_U[mu]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }

  // Only compute susy-breaking scalar potential term if bmass non-zero
  if (bmass > IMAG_TOL) {
    Real dmu;
#ifdef EIG_POT
    dmu = 2.0 * bmass * bmass;
#else
    Real tr;
    dmu = 2.0 * one_ov_N * bmass * bmass;
#endif

    FORALLSITES(i, s) {
      FORALLDIR(mu) {
#ifdef EIG_POT
        // Ubar_a(x) [U_a(x) Ubar_a(x) - I]
        mult_na(&(s->link[mu]), &(s->link[mu]), &tmat);
        scalar_add_diag(&tmat, -1.0);
        scalar_mult_an_sum(&(s->link[mu]), &tmat, dmu, &(s->f_U[mu]));
#else
        // Ubar_a(x) (Tr[U_a(x) Ubar_a(x)] / N - 1)
        tr = one_ov_N * realtrace(&(s->link[mu]), &(s->link[mu])) - 1.0;
        tr *= dmu;
        scalar_mult_sum_adj_matrix(&(s->link[mu]), tr, &(s->f_U[mu]));
#endif
      }
    }
  }

#ifdef VARPHI
  //IMPORTANT NOTE:
  //f_var is supposed to be NxF
//TODO:Use the 2step gathers
  FORALLSITES(i, s) {
    //Init
    fun_clear_mat(&(s->f_var));
    FORALLDIR(mu) {

      tag0[0] = start_gather_site(F_OFFSET(link[mu]), sizeof(matrix),
                                  allminus_offset, EVENANDODD, local_pt[0][0]);
      wait_gather(tag0[0]);
      //tmat = &((matrix *) local_pt[0][0]);

      //fun_at_mult_an_sub(&(tmat), &(Fmu5[i]), &(s->f_var));
      fun_at_mult_an_sub((matrix*) (local_pt[0][0]), &(Fmu5[i]), &(s->f_var));
      fun_o_sum(&(Fmu5[i]), &(s->f_var));
      cleanup_gather(tag0[0]);  
    }
  }
#endif

  // Finally take adjoint and update the momentum
  // Include overall factor of kappa = N / (4lambda)
  // Subtract to reproduce -Adj(f_U)
  // Compute average gauge force in same loop
  tr = kappa * eps;
  FORALLSITES(i, s) {
    FORALLDIR(mu) {
      scalar_mult_dif_adj_matrix(&(s->f_U[mu]), tr, &(s->mom[mu]));
      returnit += realtrace(&(s->f_U[mu]), &(s->f_U[mu]));
    }
#ifdef VARPHI
    fun_scaled_adjoint_sum(&(s->f_var), -tr, &(s->varmom));
    returnit += fun_nn_realtrace(&(s->f_var),&(s->f_var));
#endif
  }
  g_doublesum(&returnit);
  returnit *= kappa * kappa;

  // Add in force from separate determinant term if kappa_u1 non-zero
  if (kappa_u1 > IMAG_TOL)
    returnit += det_force(eps);

  return (eps * sqrt(returnit) / volume);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Separate routines for each term in the fermion force
// All called by assemble_fermion_force below
// First plaq--vol piece: epsilon_{cde} theta D_c^+ chi_de
// Use tempmat, tempmat2 and UpsiU[0] for temporary storage
#ifdef PV
void F_PtoV(matrix *plaq_sol[NPLAQ], matrix *vol_psol) {
  register int i;
  register site *s;
  char **local_pt[2][3];
  int a, b, c, i_ab, j, gather, next, flip = 0;
  Real tr;
  msg_tag *tag0[2], *tag1[2], *tag2[2];
  matrix tmat, tmat2;

  for (a = 0; a < 3; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[3 + a];
  }

  //Start first set of gathers
  //From setup_lambda, b>a and b!=c, a!=c

  a = FQ_lookup[0][0];
  b = FQ_lookup[0][1];
  c = FQ_lookup[0][2];
  i_ab = plaq_index[a][b];

  tag0[0] = start_gather_field(plaq_sol[i_ab], sizeof(matrix),
                            goffset[c], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(plaq_sol[i_ab], sizeof(matrix),
                            FQ_d1[0], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(vol_psol, sizeof(matrix),
                            FQ_d1[0], EVENANDODD, local_pt[0][2]);

  for (j = 0; j < NTERMS; j++) {
    gather = (flip + 1) % 2;
    if ( j < NTERMS - 1) {
      next = j + 1;
      a = FQ_lookup[next][0];
      b = FQ_lookup[next][1];
      c = FQ_lookup[next][2];
      i_ab = plaq_index[a][b];

      tag0[gather] = start_gather_field(plaq_sol[i_ab], sizeof(matrix),
                                goffset[c], EVENANDODD, local_pt[gather][0]);
      tag1[gather] = start_gather_field(plaq_sol[i_ab], sizeof(matrix),
                                FQ_d1[next], EVENANDODD, local_pt[gather][1]);
      tag2[gather] = start_gather_field(vol_psol, sizeof(matrix),
                                FQ_d1[next], EVENANDODD, local_pt[gather][2]);
    }

    a = FQ_lookup[j][0];
    b = FQ_lookup[j][1];
    c = FQ_lookup[j][2];
    i_ab = plaq_index[a][b];
    tr = perm[a][b][c];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);

    FORALLSITES(i, s) {
      mult_nn((matrix *)(local_pt[flip][1][i]),
              (matrix *)(local_pt[flip][2][i]), &(tmat2));
      scalar_mult_matrix((matrix *)(local_pt[flip][0][i]), tr * s->bc[c], &tmat);
      mult_nn(&(vol_psol[i]), &tmat, &(tempmat[i]));

      scalar_mult_sum_matrix(&(tmat2), -tr, &(tempmat[i]));
      scalar_mult_sum_adj_matrix(&(tempmat[i]), -0.5, &(s->f_U[c]));
    }

    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    flip   = (flip + 1) % 2;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// First plaq--vol piece: epsilon_{cde} theta D_c^+ chi_de
// Use tempmat, tempmat2 and UpsiU[0] for temporary storage
// TODO:fix this, doesn't pass unit test
void F_VtoP(matrix *plaq_psol[NPLAQ], matrix *vol_sol) {
  register int i;
  register site *s;
  int a, b, c, j, i_ab, next, gather, flip = 0;
  Real tr;
  msg_tag *tag0[2], *tag1[2], *tag2[2];
  matrix tmat, tmat2;
  char **local_pt[2][3];

  for (a = 0; a < 3; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[3 + a];
  }

  //Start first set of gathers
  //From setup_lambda, b>a and b!=c, a!=c

  a = FQ_lookup[0][0];
  b = FQ_lookup[0][1];
  c = FQ_lookup[0][2];
  i_ab = plaq_index[a][b];

  tag0[0] = start_gather_field(plaq_psol[i_ab], sizeof(matrix),
                            goffset[c], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(plaq_psol[i_ab], sizeof(matrix),
                            FQ_d1[0], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(vol_sol, sizeof(matrix),
                            FQ_d1[0], EVENANDODD, local_pt[0][2]);

  for (j = 0; j < NTERMS; j++) {
    gather = (flip + 1) % 2;
    if (j < NTERMS - 1) {
      next = j + 1;
      a = FQ_lookup[next][0];
      b = FQ_lookup[next][1];
      c = FQ_lookup[next][2];
      i_ab = plaq_index[a][b];

      tag0[gather] = start_gather_field(plaq_psol[i_ab], sizeof(matrix),
                                goffset[c], EVENANDODD, local_pt[gather][0]);
      tag1[gather] = start_gather_field(plaq_psol[i_ab], sizeof(matrix),
                                FQ_d1[next], EVENANDODD, local_pt[gather][1]);
      tag2[gather] = start_gather_field(vol_sol, sizeof(matrix),
                                FQ_d1[next], EVENANDODD, local_pt[gather][2]);
    }
    a = FQ_lookup[j][0];
    b = FQ_lookup[j][1];
    c = FQ_lookup[j][2];
    i_ab = plaq_index[a][b];

    tr = perm[a][b][c];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);

    FORALLSITES(i, s) {
      mult_nn((matrix *) local_pt[flip][1][i]
            , (matrix *) local_pt[flip][2][i], &(tmat2));

      scalar_mult_matrix(&(tmat2), tr, &(tempmat[i]));

      scalar_mult_matrix(&(vol_sol[i]), -tr * s->bc[c], &tmat);
      mult_nn_sum(&tmat, (matrix *)(local_pt[flip][0][i]), &(tempmat[i]));
      scalar_mult_sum_adj_matrix(&(tempmat[i]), -0.5, &(s->f_U[c]));
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    flip   = (flip + 1) % 2;
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant contributions to the fermion force
// Use Uinv, Udag_inv, UpsiU, Tr_Uinv and tr_dest for temporary storage
// Also use tempdet for temporary storage
// The accumulator names refer to the corresponding derivatives
// Assume compute_plaqdet() has already been run
// Appropriate adjoints set up in assemble_fermion_force
// A bit more code reuse may be possible
#ifdef SL
void detF(matrix *eta, matrix *psi[NUMLINK], int sign) {
  register int i;
  register site *s;
  int a, b, opp_b;
  Real localG = 0.5 * C2 * G;
  complex tc, tc2;
  msg_tag *mtag[8];
  matrix tmat;

  // Check sign while giving G proper sign
  if (sign == 1)
    localG *= -1.0;
  else if (sign != -1) {
    node0_printf("Error: incorrect sign in detF: %d\n", sign);
    terminate(1);
  }

  // Set up and store some basic ingredients
  FORALLSITES(i, s)
    tr_eta[i] = trace(&(eta[i]));
  // Need all directions for upcoming sums
  FORALLDIR(a) {
    // U_a(x)^{-1} and Udag_a(x)^{-1} are already in Uinv[a] and Udag_inv[a]
    // Save Tr[U_a(x)^{-1} psi_a(x) in Tr_Uinv[a]
    // Save U_a(x)^{-1} psi_a(x) U_a(x)^{-1} in UpsiU[a]
    FORALLSITES(i, s) {
      mult_nn(&(Uinv[a][i]), &(psi[a][i]), &tmat);
      mult_nn(&tmat, &(Uinv[a][i]), &(UpsiU[a][i]));
      Tr_Uinv[a][i] = trace(&tmat);

      // tempdet holds Tr[eta(x)] plaqdet[a][b](x)
      for (b = a + 1; b < NUMLINK; b++) {
        CMUL(tr_eta[i], plaqdet[a][b][i], tempdet[a][b][i]);
        CMUL(tr_eta[i], plaqdet[b][a][i], tempdet[b][a][i]);
      }
    }
  }

  // Now we are ready to gather, accumulate and add to force
  complex *plaq_term = malloc(sizeof *plaq_term * sites_on_node);
  complex *inv_term = malloc(sizeof *inv_term * sites_on_node);
  complex *adj_term = malloc(sizeof *adj_term * sites_on_node);

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
        CMULREAL(*((complex *)(gen_pt[0][i])), s->bc[a], tr_dest[i]);
      cleanup_gather(mtag[0]);
      mtag[0] = start_gather_field(tr_dest, sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[0]);

      // Now accumulate all three terms
      opp_b = OPP_LDIR(b);
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
        tc = *((complex *)(gen_pt[5][i]));
        tc2.real = Tr_Uinv[a][i].real + s->bc[a] * tc.real;
        tc2.imag = Tr_Uinv[a][i].imag + s->bc[a] * tc.imag;
        plaq_term[i].real += tempdet[b][a][i].real * tc2.real
                           - tempdet[b][a][i].imag * tc2.imag;
        plaq_term[i].imag += tempdet[b][a][i].imag * tc2.real
                           + tempdet[b][a][i].real * tc2.imag;

        // D[a][b](x - b) {T[a](x) + T[b](x - b)}
        // gen_pt[6] is T[b](x - b)
        tc = *((complex *)(gen_pt[6][i]));
        tc2.real = tc.real + s->bc[opp_b] * Tr_Uinv[a][i].real;
        tc2.imag = tc.imag + s->bc[opp_b] * Tr_Uinv[a][i].imag;
        // gen_pt[1] is D[a][b](x - b)
        tc = *((complex *)(gen_pt[1][i]));
        plaq_term[i].real += tc.real * tc2.real - tc.imag * tc2.imag;
        plaq_term[i].imag += tc.imag * tc2.real + tc.real * tc2.imag;

        // Accumulate adj_term
        // D[a][b](x) {T[b](x) + T[a](x + b)}
        // gen_pt[3] is T[a](x + b)
        tc = *((complex *)(gen_pt[3][i]));
        tc2.real = Tr_Uinv[b][i].real + s->bc[b] * tc.real;
        tc2.imag = Tr_Uinv[b][i].imag + s->bc[b] * tc.imag;
        adj_term[i].real += tempdet[a][b][i].real * tc2.real
                          - tempdet[a][b][i].imag * tc2.imag;
        adj_term[i].imag += tempdet[a][b][i].imag * tc2.real
                          + tempdet[a][b][i].real * tc2.imag;

        // D[b][a](x - b) {T[a](x - b) + T[b](x + a - b) bc[a](x - b)}
        // gen_pt[0] is T[b](x + a - b)
        // gen_pt[4] is T[a](x - b)
        CADD(*((complex *)(gen_pt[0][i])), *((complex *)(gen_pt[4][i])), tc);
        // gen_pt[2] is D[b][a](x - b)
        tc2 = *((complex *)(gen_pt[2][i]));
        adj_term[i].real += tc.real * tc2.real - tc.imag * tc2.imag;
        adj_term[i].imag += tc.imag * tc2.real + tc.real * tc2.imag;

        // Accumulate inv_term = sum_b D[b][a](x) + D[a][b](x - b)
        // gen_pt[1] is D[a][b](x - b)
        tc = *((complex *)(gen_pt[1][i]));
        inv_term[i].real += tempdet[b][a][i].real + s->bc[opp_b] * tc.real;
        inv_term[i].imag += tempdet[b][a][i].imag + s->bc[opp_b] * tc.imag;
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
    FORALLSITES(i, s) {
      // Start with plaq_term hitting U_a(x)^{-1}
      CMULREAL(plaq_term[i], localG, tc);
      c_scalar_mult_sum_mat(&(Uinv[a][i]), &tc, &(s->f_U[a]));

      // Add adj_term hitting Udag_a(x)^{-1} followed by adjoint
      CMULREAL(adj_term[i], localG, tc);
      c_scalar_mult_sum_adj_mat(&(Udag_inv[a][i]), &tc, &(s->f_U[a]));

      // Finally subtract inv_term hitting U_a(x)^{-1} psi_a(x) U_a(x)^{-1}
      CMULREAL(inv_term[i], localG, tc);
      c_scalar_mult_dif_mat(&(UpsiU[a][i]), &tc, &(s->f_U[a]));
    }
  }
  free(plaq_term);
  free(inv_term);
  free(adj_term);
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Assemble fermion contributions to gauge link force,
//   f_U = Adj(Ms).D_U M(U, Ub).s - Adj[Adj(Ms).D_Ub M(U, Ub).s]
// "s" is sol while "Ms" is psol
// Copy these into persistent matrices for easier gathering
// Use tempmat, tempmat2, UpsiU, Tr_Uinv,
// tr_dest and Ddet[012] for temporary storage
// (many through calls to detF)
// Almost the same as 4d. Only diff is closed term
#ifndef PUREGAUGE
void assemble_fermion_force(Twist_Fermion *sol, Twist_Fermion *psol) {
  register int i;
  register site *s;
  char **local_pt[2][2];
  int mu, nu, a, b, gather, flip = 0, index, next;
  msg_tag *mtag[NUMLINK], *tag0[2], *tag1[2];
  matrix *mat[2], tmat;

  for (mu = 0; mu < 2; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[2 + mu];
  }
  mat[0] = tempmat;
  mat[1] = tempmat2;

  // For gathering it is convenient to copy the input Twist_Fermions
  // into persistent site, link and plaq fermions
  // We can reuse "src" and "dest" for this storage,
  // corresponding to "sol" and "psol", respectively
  FORALLSITES(i, s) {
    mat_copy(&(sol[i].Fsite), &(site_src[i]));
    adjoint(&(psol[i].Fsite), &(site_dest[i]));
    mat_copy(&(sol[i].Fvolume), &(volume_src[i]));
    adjoint(&(psol[i].Fvolume), &(volume_dest[i]));
    FORALLDIR(mu) {
      mat_copy(&(sol[i].Flink[mu]), &(link_src[mu][i]));
      adjoint(&(psol[i].Flink[mu]), &(link_dest[mu][i]));
    }
    for (mu = 0; mu < NPLAQ; mu++) {
      mat_copy(&(sol[i].Fplaq[mu]), &(plaq_src[mu][i]));
      adjoint(&(psol[i].Fplaq[mu]), &(plaq_dest[mu][i]));
    }
  }

#ifdef SL
  // Accumulate both terms in UpsiU[mu], use to initialize f_U[mu]
  // First calculate DUbar on eta Dbar_mu psi_mu (LtoS)
  // [psi_mu(x) eta(x + mu) - eta(x) psi_mu(x)]^dag
  mtag[0] = start_gather_field(site_dest, sizeof(matrix),
                               goffset[0], EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(site_dest, sizeof(matrix),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);
    }
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_matrix((matrix *)(gen_pt[mu][i]), s->bc[mu], &tmat);
      mult_nn(&(link_src[mu][i]), &tmat, &(UpsiU[mu][i]));   // Initialize
      mult_nn_dif(&(site_dest[i]), &(link_src[mu][i]), &(UpsiU[mu][i]));
    }
    cleanup_gather(mtag[mu]);
  }

  // 2nd term, DUbar on psi_mu Dbar_mu eta (StoL)
  // [eta(x) psi_mu(x) - psi_mu(x) eta(x + mu)]^dag
  mtag[0] = start_gather_field(site_src, sizeof(matrix),
                               goffset[0], EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(site_src, sizeof(matrix),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);
    }
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_matrix((matrix *)(gen_pt[mu][i]), s->bc[mu], &tmat);
      mult_nn_dif(&(link_dest[mu][i]), &tmat, &(UpsiU[mu][i]));
      mult_nn_sum(&(site_src[i]), &(link_dest[mu][i]), &(UpsiU[mu][i]));

      // Initialize the force collectors---done with UpsiU[mu]
      scalar_mult_adj_matrix(&(UpsiU[mu][i]), 0.5, &(s->f_U[mu]));
    }
    cleanup_gather(mtag[mu]);
  }
#endif
#ifdef LP
  // Now calculate DU on chi_{munu} D_mu(U) psi_nu
  // Start first set of gathers (mu = 0 and nu = 1)
  tag0[0] = start_gather_field(link_src[1], sizeof(matrix),
                               goffset[0], EVENANDODD, local_pt[0][0]);

  // Prepare and gather other term in tempmat*
  index = plaq_index[0][1];
  FORALLSITES(i, s)     // mu = 0 < nu = 1
    mult_nn(&(plaq_dest[index][i]), &(link_src[1][i]), &(mat[0][i]));
  tag1[0] = start_gather_field(mat[0], sizeof(matrix),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]);

  // Main loop
  FORALLDIR(mu) {
    FORALLDIR(nu) {
      if (mu == nu)
        continue;

      gather = (flip + 1) % 2;
      if (mu < NUMLINK - 1 || nu < NUMLINK - 2) {   // Start next gathers
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
        tag0[gather] = start_gather_field(link_src[b], sizeof(matrix),
                                          goffset[a], EVENANDODD,
                                          local_pt[gather][0]);

        // Prepare and gather other term in tempmat*
        next = plaq_index[a][b];
        FORALLSITES(i, s) {
          if (a > b) {      // plaq_dest is anti-symmetric under a <--> b
            scalar_mult_matrix(&(plaq_dest[next][i]), -1.0, &tmat);
          }                 // Suppress compiler error
          else
            mat_copy(&(plaq_dest[next][i]), &tmat);
          mult_nn(&tmat, &(link_src[b][i]), &(mat[gather][i]));
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix),
                                          goffset[b] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu) {    // plaq_dest is anti-symmetric under mu <--> nu
          scalar_mult_matrix((matrix *)(local_pt[flip][0][i]),
                             s->bc[mu], &tmat);
        }                 // Suppress compiler error
        else
          scalar_mult_matrix((matrix *)(local_pt[flip][0][i]),
                             -1.0 * s->bc[mu], &tmat);

        mult_nn_sum(&tmat, &(plaq_dest[index][i]), &(s->f_U[mu]));
        sum_matrix((matrix *)(local_pt[flip][1][i]), &(s->f_U[mu]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }

  // 2nd term
  // Start first set of gathers (mu = 0 and nu = 1)
  tag0[0] = start_gather_field(link_dest[1], sizeof(matrix),
                               goffset[0], EVENANDODD, local_pt[0][0]);

  // Prepare and gather other term in tempmat*
  index = plaq_index[0][1];
  FORALLSITES(i, s)     // mu = 0 < nu = 1
    mult_nn(&(plaq_src[index][i]), &(link_dest[1][i]), &(mat[0][i]));
  tag1[0] = start_gather_field(mat[0], sizeof(matrix),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]);

  // Main loop
  FORALLDIR(mu) {
    FORALLDIR(nu) {
      if (mu == nu)
        continue;

      gather = (flip + 1) % 2;
      if (mu < NUMLINK - 1 || nu < NUMLINK - 2) {   // Start next gathers
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
        tag0[gather] = start_gather_field(link_dest[b], sizeof(matrix),
                                          goffset[a], EVENANDODD,
                                          local_pt[gather][0]);

        // Prepare and gather other term in tempmat*
        next = plaq_index[a][b];
        FORALLSITES(i, s) {
          if (a > b) {      // plaq_src is anti-symmetric under a <--> b
            scalar_mult_matrix(&(plaq_src[next][i]), -1.0, &tmat);
          }                 // Suppress compiler error
          else
            mat_copy(&(plaq_src[next][i]), &tmat);
          mult_nn(&tmat, &(link_dest[b][i]), &(mat[gather][i]));
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix),
                                          goffset[b] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu) {    // plaq_src is anti-symmetric under mu <--> nu
          scalar_mult_matrix(&(plaq_src[index][i]), -1.0 * s->bc[mu], &tmat);
        }                 // Suppress compiler error
        else
          scalar_mult_matrix(&(plaq_src[index][i]), s->bc[mu], &tmat);

        mult_nn_sum((matrix *)(local_pt[flip][0][i]), &tmat, &(s->f_U[mu]));
        dif_matrix((matrix *)(local_pt[flip][1][i]), &(s->f_U[mu]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }
#endif
#ifdef SL
  // Plaquette determinant contributions if G is non-zero
  if (doG) {
    // First connect link_src with site_dest[DIMF - 1]^dag (LtoS)
    detF(site_dest, link_src, PLUS);

    // Second connect site_src[DIMF - 1] with link_dest^dag (StoL)
    detF(site_src, link_dest, MINUS);
  }
#endif

#ifdef PV
  if (NUMLINK != 3) {
    node0_printf("ERROR: NUMLINK IS %d != 3\n", NUMLINK);
    terminate(1);
  }
  F_PtoV(plaq_src, volume_dest);
  F_VtoP(plaq_dest, volume_src);
#endif
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Update the momenta with the fermion force
// Assume that the multiCG has been run, with the solution in sol[j]
// Accumulate f_U for each pole into fullforce, add to momenta
// Use fullforce-->Fmunu and tempTF for temporary storage
// (Calls assemble_fermion_force, which uses many more temporaries)
double fermion_force(Real eps, Twist_Fermion *src, Twist_Fermion **sol) {
  register int i;
  register site *s;
  int mu, n;
  double returnit = 0.0;
  matrix **fullforce = malloc(sizeof(matrix*) * NUMLINK);

  FORALLDIR(mu)
    fullforce[mu] = Fmunu[mu];    // Use Fmunu for temporary storage

  // Initialize fullforce[mu]
  fermion_op(sol[0], tempTF, PLUS);
  FORALLSITES(i, s)
    scalar_mult_TF(&(tempTF[i]), amp4[0], &(tempTF[i]));
  assemble_fermion_force(sol[0], tempTF);
  FORALLDIR(mu) {
    FORALLSITES(i, s)
      adjoint(&(s->f_U[mu]), &(fullforce[mu][i]));
  }
  for (n = 1; n < Norder; n++) {
    fermion_op(sol[n], tempTF, PLUS);
    // Makes sense to multiply here by amp4[n]...
    FORALLSITES(i, s)
      scalar_mult_TF(&(tempTF[i]), amp4[n], &(tempTF[i]));

    assemble_fermion_force(sol[n], tempTF);
    // Take adjoint but don't negate yet...
    FORALLDIR(mu) {
      FORALLSITES(i, s)
        sum_adj_matrix(&(s->f_U[mu]), &(fullforce[mu][i]));
    }
  }

  // Update the momentum from the fermion force -- sum or eps
  // Opposite sign as to gauge force,
  // because dS_G / dU = 2F_g while ds_F / dU = -2F_f
  // Move negation here as well, though adjoint remains above
  FORALLSITES(i, s) {
    FORALLDIR(mu) {
      scalar_mult_dif_matrix(&(fullforce[mu][i]), eps, &(s->mom[mu]));
      returnit += realtrace(&(fullforce[mu][i]), &(fullforce[mu][i]));
    }
  }
  g_doublesum(&returnit);

  free(fullforce);

  // Reset Fmunu
  compute_Fmunu();
  return (eps * sqrt(returnit) / volume);
}
#endif
// -----------------------------------------------------------------
