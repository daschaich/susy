// -----------------------------------------------------------------
// Update the momentum matrices
// Uncomment to print out debugging messages
//#define FORCE_DEBUG
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update mom[NUMLINK] with the gauge force
// Include tunable coefficient C2 in the d^2 term of the action
double gauge_force(Real eps) {
  register int i, mu, nu;
  register site *s;
  register Real eb3, eb3UdU, normK, normS, trace_sub;
  double returnit = 0;
  complex trUdU;
  msg_tag *tag[NUMLINK], *tag0, *tag1;
  su3_matrix_f tmat1, tmat2, tmat3, tmatK, UdU;

  // Contribution from (D_a U_a)^2 / 2 term
  compute_DmuUmu();
  // To reproduce:
  //   f_U.set(x, mu, f_U.get(x, mu) + Udag.get(x, mu) * DmuUmu.get(x));
  //   f_U.set(x, mu, f_U.get(x, mu) - DmuUmu.get(x + e_mu) * Udag.get(x, mu));
  tag[0] = start_gather_site(F_OFFSET(DmuUmu), sizeof(su3_matrix_f),
                             goffset[0], EVENANDODD, gen_pt[0]);
  for (mu = 0; mu < NUMLINK; mu++) {
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_site(F_OFFSET(DmuUmu), sizeof(su3_matrix_f),
                                      goffset[mu + 1], EVENANDODD,
                                      gen_pt[mu + 1]);

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      mult_su3_an_f(&(s->linkf[mu]), &(s->DmuUmu), &tmat1);
      mult_su3_na_f((su3_matrix_f *)gen_pt[mu][i], &(s->linkf[mu]), &tmat2);
      sub_su3_matrix_f(&tmat1, &tmat2, &tmat3);
      scalar_mult_su3_matrix_f(&tmat3, C2, &(s->f_U[mu]));    // Initialize
    }
    cleanup_gather(tag[mu]);
  }

  // Contribution from Fbar_{ab} F_{ab} term
  compute_Fmunu();
  // To reproduce:
  //   f_U.set(x, mu, f_U.get(x, mu) + 2 * U.get(x + e_mu, nu) * Adj(Fmunu.get(x, mu, nu)));
  //   f_U.set(x, mu, f_U.get(x, mu) - 2 * Adj(Fmunu.get(x - e_nu, mu, nu)) * U.get(x - e_nu, nu));
  for (mu = 0; mu < NUMLINK; mu++) {
    for (nu = 0; nu < NUMLINK; nu++) {
      if (mu == nu)
        continue;
      FORALLSITES(i, s)
        mult_su3_an_f(&(s->Fmunu[mu][nu]), &(s->linkf[nu]), &(s->tempmat2));

      tag0 = start_gather_site(F_OFFSET(linkf[nu]), sizeof(su3_matrix_f),
                               goffset[mu], EVENANDODD, gen_pt[0]);
      tag1 = start_gather_site(F_OFFSET(tempmat2), sizeof(su3_matrix_f),
                               goffset[nu] + 1, EVENANDODD, gen_pt[1]);
      wait_gather(tag0);
      wait_gather(tag1);

      FORALLSITES(i, s) {
        mult_su3_na_f((su3_matrix_f *)gen_pt[0][i], &(s->Fmunu[mu][nu]),
                      &tmat1);
        sub_su3_matrix_f(&tmat1, (su3_matrix_f *)gen_pt[1][i], &tmat3);
        scalar_mult_add_su3_matrix_f(&(s->f_U[mu]), &tmat3, 2.0, &(s->f_U[mu]));
      }
      cleanup_gather(tag0);
      cleanup_gather(tag1);
    }
  }

  // Factor of kapp = Nc / (2lambda) on both (D_a U_a)^2 and F^2 terms
  for (mu = 0; mu < NUMLINK; mu++) {
    FORALLSITES(i, s)
      scalar_mult_su3_matrix_f(&(s->f_U[mu]), kappa, &(s->f_U[mu]));
  }

  // Only compute U(1) mass term if non-zero -- note factor of kappa
  eb3 = 2.0 * kappa * bmass * bmass / (Real)(NCOL * NCOL);
  if (eb3 > 1.e-8) {
    for (mu = 0; mu < NUMLINK; mu++) {
      FORALLSITES(i, s) {
        mult_su3_an_f(&(s->linkf[mu]), &(s->linkf[mu]), &UdU);
        trUdU = trace_su3_f(&UdU);
        eb3UdU = eb3 * (trUdU.real -(Real)NCOL);

        su3_adjoint_f(&(s->linkf[mu]), &tmat1);
        scalar_mult_add_su3_matrix_f(&(s->f_U[mu]), &tmat1, eb3UdU,
                                     &(s->f_U[mu]));
      }
    }
  }

  // Only compute Konishi and SUGRA contributions if non-zero
  // No factor of kappa
  if (CK > 1.0e-8 || CS > 1.0e-8) {
    normK = CK;
    normS = CS / (Real)(NUMLINK * NUMLINK);
    trace_sub = -1.0 / (Real)NUMLINK;
    compute_Bmu();

    // The SUGRA trace subtraction makes this a little awkward
    // Store Konishi and SUGRA contributions in the site
    // (tempmat1 and Fmunu[mu][nu], respectively)
    // and only add SUGRA to force in a second loop over mu and nu
    // As always, average over all 25 SUGRA components
    FORALLSITES(i, s) {
      clear_su3mat_f(&(s->tempmat1));   // Konishi contribution
      for (mu = 0; mu < NUMLINK; mu++) {
        for (nu = 0; nu < NUMLINK; nu++)
          clear_su3mat_f(&(s->Fmunu[mu][nu]));   // SUGRA contributions
      }
    }

    // Accumulate and store Konishi and SUGRA contributions
    // We're after f_U[mu], which involves Ubar_mu hit with all B_nu
    // So we need to gather every B_nu from site x + mu
    for (mu = 0; mu < NUMLINK; mu++) {
      for (nu = 0; nu < NUMLINK; nu++) {
        tag0 = start_gather_site(F_OFFSET(B[nu]), sizeof(su3_matrix_f),
                                 goffset[mu], EVENANDODD, gen_pt[0]);
        wait_gather(tag0);
        FORALLSITES(i, s) {
          mult_su3_an_f(&(s->linkf[mu]), &(s->B[nu]), &tmat1);
          mult_su3_na_f((su3_matrix_f *)gen_pt[nu][i], &(s->linkf[mu]),
                        &tmat2);
          add_su3_matrix_f(&tmat1, &tmat2, &tmat3);

          // Accumulate Konishi and SUGRA
          if (mu == nu)
            add_su3_matrix_f(&(s->tempmat1), &tmat3, &(s->tempmat1));
//          add_su3_matrix_f(&(s->tempmat2), &tmat3, &(s->tempmat2));
        }
      }
    }
    // SUGRA trace subtraction (already summed diagonal contributions)
//    FORALLSITES(i, s)
//      sub_su3_matrix_f(&(s->tempmat2), &(s->tempmat1), &(s->tempmat2));

    // Now add Konishi and SUGRA contributions to f_U[mu]
    FORALLSITES(i, s) {
        
      //
      //
      // How to do trace subtraction ...
      //
      //
      //
      for (mu = 0; mu < NUMLINK; mu++) {
        scalar_mult_add_su3_matrix_f(&tmat3, &tmatK, trace_sub, &tmat3);
        // Update force with averaged SUGRA contribution
        scalar_mult_add_su3_matrix_f(&(s->f_U[mu]), &tmat3, normS,
                                     &(s->f_U[mu]));
      }
    }
  }

  // Finally take adjoint and update the momentum
  // Subtract to reproduce -Adj(f_U)
  for (mu = 0; mu < NUMLINK; mu++) {
    cleanup_gather(tag[mu]);
    FORALLSITES(i, s) {
      su3_adjoint_f(&(s->f_U[mu]), &tmat2);
      scalar_mult_sub_su3_matrix_f(&(s->mom[mu]), &tmat2, eps, &(s->mom[mu]));
    }
  }

  // Compute average gauge force
  for (mu = 0; mu < NUMLINK; mu++) {
    FORALLSITES(i, s)
      returnit += realtrace_su3_f(&(s->f_U[mu]), &(s->f_U[mu]));
  }
  g_doublesum(&returnit);

  // Add in force from determinant term
  returnit += det_force(eps);
  return (eps * sqrt(returnit) / volume);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Create the fermion force
void assemble_fermion_force(Twist_Fermion *sol, Twist_Fermion *psol) {
  register int i;
  register site *s;
  char **local_pt[2][4];
  int mu, nu, a, b, c, d, e, l, m, counter, gather, next, flip = 0;
  Real permm, BC;
  complex cterm;
  msg_tag *mtag[NUMLINK], *tag0[2], *tag1[2], *tag2[2], *tag3[2];
  msg_tag *mtag0 = NULL, *mtag1 = NULL;
  su3_vector *vec, *vec0, *vec1, *vec2, *vec3;
  su3_matrix_f tmpmat, tmpmat2;
  su3_matrix_f *tempmat = malloc(sites_on_node * sizeof(*tempmat));

  // Copy Twist_Fermions into persistent site, link and plaquette fermions
  FORALLSITES(i, s) {
    for (mu = 0; mu < NUMLINK; mu++) {
      for (nu = 0; nu < NUMLINK; nu++) {
        su3vec_copy(&(sol[i].Fplaq[mu][nu]), &(s->plaq_sol[mu][nu]));
        su3vec_copy(&(psol[i].Fplaq[mu][nu]), &(s->plaq_psol[mu][nu]));
      }
    }
    for (mu = 0; mu < NUMLINK; mu++) {
      su3vec_copy(&(sol[i].Flink[mu]), &(s->link_sol[mu]));
      su3vec_copy(&(psol[i].Flink[mu]), &(s->link_psol[mu]));
    }
    su3vec_copy(&(sol[i].Fsite), &(s->site_sol));
    su3vec_copy(&(psol[i].Fsite), &(s->site_psol));
  }

  // Clear the force collectors
  for (mu = 0; mu < NUMLINK; mu++) {
    FORALLSITES(i, s)
      clear_su3mat_f(&(s->f_U[mu]));
  }

  // Compute fermion contribution to gauge link force,
  // f_U = Adj(Ms).D_U M(U, Ub) s - Adj[Adj(Ms). D_Ub M(U, Ub) s]
  // First calculate DU on chi_{munu} D_mu(U) psi_nu
#ifdef VP
  for (mu = 0; mu < NUMLINK; mu++) {
    for (nu = 0; nu < NUMLINK; nu++) {
      if (mu == nu)
        continue;

      mtag0 = start_gather_site(F_OFFSET(link_sol[nu]), sizeof(su3_vector),
                                 goffset[mu], EVENANDODD, gen_pt[0]);
      wait_gather(mtag0);
      FORALLSITES(i, s) {
        vec = (su3_vector *)(gen_pt[0][i]);
        for (a = 0; a < NUMGEN; a++) {
          for (b = 0; b < NUMGEN; b++) {
            CMULJ_((s->plaq_psol[mu][nu]).c[a], vec->c[b], cterm);
            CNEGATE(cterm, cterm);
            CMULREAL(cterm, s->bc1[mu], cterm);
            c_scalar_mult_add_su3mat_f(&(s->f_U[mu]), &(Lambda_prod[b][a]),
                                       &cterm, &(s->f_U[mu]));
          }
        }
      }
      cleanup_gather(mtag0);

      FORALLSITES(i, s) {
        clear_su3mat_f(&(tempmat[i]));
        for (a = 0; a < NUMGEN; a++) {
          for (b = 0; b < NUMGEN; b++) {
            CMULJ_((s->plaq_psol[mu][nu]).c[a],
                   (s->link_sol[nu]).c[b], cterm);
            c_scalar_mult_add_su3mat_f(&(tempmat[i]), &(Lambda_prod[a][b]),
                                       &cterm, &(tempmat[i]));
          }
        }
      }
      mtag1 = start_gather_field(tempmat, sizeof(su3_matrix_f),
                                 goffset[nu] + 1, EVENANDODD, gen_pt[1]);
      wait_gather(mtag1);
      FORALLSITES(i, s)
        add_su3_matrix_f(&(s->f_U[mu]), (su3_matrix_f *)(gen_pt[1][i]),
                         &(s->f_U[mu]));

      cleanup_gather(mtag1);
    }
  }

  // 2nd term
  for (mu = 0; mu < NUMLINK; mu++) {
    for (nu = 0; nu < NUMLINK; nu++) {
      if (mu == nu)
        continue;

      mtag0 = start_gather_site(F_OFFSET(link_psol[nu]), sizeof(su3_vector),
                                 goffset[mu], EVENANDODD, gen_pt[0]);

      FORALLSITES(i, s) {
        clear_su3mat_f(&(tempmat[i]));
        for (a = 0; a < NUMGEN; a++) {
          for (b = 0; b < NUMGEN; b++) {
            CMULJ_((s->link_psol[nu]).c[a],
                   (s->plaq_sol[mu][nu]).c[b], cterm);
            c_scalar_mult_add_su3mat_f(&(tempmat[i]), &(Lambda_prod[b][a]),
                                       &cterm, &(tempmat[i]));
          }
        }
      }
      mtag1 = start_gather_field(tempmat, sizeof(su3_matrix_f),
                                 goffset[nu] + 1, EVENANDODD, gen_pt[1]);

      wait_gather(mtag1);
      FORALLSITES(i, s)
        sub_su3_matrix_f(&(s->f_U[mu]), (su3_matrix_f *)(gen_pt[1][i]),
                         &(s->f_U[mu]));

      cleanup_gather(mtag1);

      wait_gather(mtag0);
      FORALLSITES(i, s) {
        vec = (su3_vector *)(gen_pt[0][i]);
        for (a = 0; a < NUMGEN; a++) {
          for (b = 0; b < NUMGEN; b++) {
            CMULJ_((vec->c[a]), (s->plaq_sol[mu][nu]).c[b], cterm);
            CMULREAL(cterm, s->bc1[mu], cterm);
            c_scalar_mult_add_su3mat_f(&(s->f_U[mu]), &(Lambda_prod[a][b]),
                                       &cterm, &(s->f_U[mu]));
          }
        }
      }
      cleanup_gather(mtag0);
    }
  }
#endif
#ifdef SV
  // 3rd term, DUb on psi_mu Db_mu(U) eta
  mtag[0] = start_gather_site(F_OFFSET(site_psol), sizeof(su3_vector),
                              goffset[0], EVENANDODD, gen_pt[0]);
  for (mu = 0; mu < NUMLINK; mu++) {
    if (mu < NUMLINK - 1)
      mtag[mu + 1] = start_gather_site(F_OFFSET(site_psol),
                           sizeof(su3_vector), goffset[mu + 1],
                           EVENANDODD, gen_pt[mu + 1]);

    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      clear_su3mat_f(&tmpmat);
      vec = (su3_vector *)(gen_pt[mu][i]);
      for (a = 0; a < NUMGEN; a++) {
        for (b = 0; b < NUMGEN; b++) {
          CMULJ_((s->site_psol).c[a],
                 (s->link_sol[mu]).c[b], cterm);
          CNEGATE(cterm, cterm);
          c_scalar_mult_add_su3mat_f(&tmpmat, &(Lambda_prod[a][b]),
                                     &cterm, &tmpmat);

          CMULJ_((vec->c[a]), (s->link_sol[mu]).c[b], cterm);
          CMULREAL(cterm, s->bc1[mu], cterm);
          c_scalar_mult_add_su3mat_f(&tmpmat, &(Lambda_prod[b][a]),
                                     &cterm, &tmpmat);
        }
      }
      su3_adjoint_f(&tmpmat, &tmpmat2);
      scalar_mult_add_su3_matrix_f(&(s->f_U[mu]), &tmpmat2, 0.5,
                                   &(s->f_U[mu]));
    }
    cleanup_gather(mtag[mu]);
  }

  // 4th term
  mtag[0] = start_gather_site(F_OFFSET(site_sol), sizeof(su3_vector),
                              goffset[0], EVENANDODD, gen_pt[0]);
  for (mu = 0; mu < NUMLINK; mu++) {
    if (mu < NUMLINK - 1)
      mtag[mu + 1] = start_gather_site(F_OFFSET(site_sol),
                           sizeof(su3_vector), goffset[mu + 1],
                           EVENANDODD, gen_pt[mu + 1]);

    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      clear_su3mat_f(&tmpmat);
      vec = (su3_vector *)(gen_pt[mu][i]);
      for (a = 0; a < NUMGEN; a++) {
        for (b = 0; b < NUMGEN; b++) {
          CMULJ_((s->link_psol[mu]).c[a], vec->c[b], cterm);
          CMULREAL(cterm, s->bc1[mu], cterm);
          CNEGATE(cterm, cterm);
          c_scalar_mult_add_su3mat_f(&tmpmat, &(Lambda_prod[a][b]),
                                     &cterm, &tmpmat);

          CMULJ_((s->link_psol[mu]).c[a],
                 (s->site_sol).c[b], cterm);
          c_scalar_mult_add_su3mat_f(&tmpmat, &(Lambda_prod[b][a]),
                                     &cterm, &tmpmat);
        }
      }
      su3_adjoint_f(&tmpmat, &tmpmat2);
      scalar_mult_add_su3_matrix_f(&(s->f_U[mu]), &tmpmat2, 0.5,
                                   &(s->f_U[mu]));
    }
    cleanup_gather(mtag[mu]);
  }
#endif
#ifdef QCLOSED
  if (NUMLINK != 5) {
    node0_printf("ERROR: NUMLINK IS %d != 5\n", NUMLINK);
    terminate(1);
  }

  for (a = 0; a < 4; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[4 + a];
  }

  // First Q-closed piece: chi_ab D_c chi_de epsilon_{abcde}
  // Start first set of gathers
  a = FQ_lookup[0][0];
  b = FQ_lookup[0][1];
  c = FQ_lookup[0][2];
  d = FQ_lookup[0][3];
  e = FQ_lookup[0][4];
  tag0[0] = start_gather_site(F_OFFSET(plaq_psol[d][e]), sizeof(su3_vector),
                              F1Q_d2[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_site(F_OFFSET(plaq_sol[a][b]), sizeof(su3_vector),
                              goffset[c], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_site(F_OFFSET(plaq_psol[d][e]), sizeof(su3_vector),
                              goffset[c], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_site(F_OFFSET(plaq_sol[a][b]), sizeof(su3_vector),
                              F1Q_d1[0], EVENANDODD, local_pt[0][3]);

  // Loop over lookup table
  for (counter = 0; counter < NTERMS; counter++) {
    gather = (flip + 1) % 2;
    if (counter < NTERMS - 1) {     // Start next set of gathers
      next = counter + 1;
      a = FQ_lookup[next][0];
      b = FQ_lookup[next][1];
      c = FQ_lookup[next][2];
      d = FQ_lookup[next][3];
      e = FQ_lookup[next][4];

      tag0[gather] = start_gather_site(F_OFFSET(plaq_psol[d][e]),
                      sizeof(su3_vector), F1Q_d2[next], EVENANDODD,
                      local_pt[gather][0]);
      tag1[gather] = start_gather_site(F_OFFSET(plaq_sol[a][b]),
                      sizeof(su3_vector), goffset[c], EVENANDODD,
                      local_pt[gather][1]);
      tag2[gather] = start_gather_site(F_OFFSET(plaq_psol[d][e]),
                      sizeof(su3_vector), goffset[c], EVENANDODD,
                      local_pt[gather][2]);
      tag3[gather] = start_gather_site(F_OFFSET(plaq_sol[a][b]),
                      sizeof(su3_vector), F1Q_d1[next], EVENANDODD,
                      local_pt[gather][3]);
    }

    // Do this set of computations while next set of gathers runs
    a = FQ_lookup[counter][0];
    b = FQ_lookup[counter][1];
    c = FQ_lookup[counter][2];
    d = FQ_lookup[counter][3];
    e = FQ_lookup[counter][4];
    permm = perm[d][e][c][a][b];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);
    wait_gather(tag3[flip]);
    FORALLSITES(i, s) {
      clear_su3mat_f(&(tempmat[i]));
      vec0 = (su3_vector *)(local_pt[flip][0][i]);
      vec1 = (su3_vector *)(local_pt[flip][1][i]);
      vec2 = (su3_vector *)(local_pt[flip][2][i]);
      vec3 = (su3_vector *)(local_pt[flip][3][i]);
      for (l = 0; l < NUMGEN; l++) {
        for (m = 0; m < NUMGEN; m++) {
          CMULJ_(vec0->c[l], vec1->c[m], cterm);
          CMULREAL(cterm, permm, cterm);
          BC = (s->bc3[a][b][c]) * (s->bc1[c]);
          CMULREAL(cterm, BC, cterm);
          c_scalar_mult_add_su3mat_f(&(tempmat[i]), &(Lambda_prod[l][m]),
                                     &cterm, &(tempmat[i]));

          CMULJ_(vec2->c[l], vec3->c[m], cterm);
          CMULREAL(cterm, -permm, cterm);
          BC = (s->bc2[OPP_LDIR(a)][OPP_LDIR(b)]) * (s->bc1[c]);
          CMULREAL(cterm, BC, cterm);
          c_scalar_mult_add_su3mat_f(&(tempmat[i]), &(Lambda_prod[m][l]),
                                     &cterm, &(tempmat[i]));
        }
      }
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    cleanup_gather(tag3[flip]);
    flip = (flip + 1) % 2;

    // Reproduces f_U.set(x, c, f_U.get(x, c) - 0.5 * Adj(tmp));
    FORALLSITES(i, s) {
      su3_adjoint_f(&(tempmat[i]), &tmpmat);
      scalar_mult_add_su3_matrix_f(&(s->f_U[c]), &tmpmat, -0.5, &(s->f_U[c]));
    }
  }

  // Second Q-closed piece
  flip = 0;   // Reset
  // Start first set of gathers
  a = FQ_lookup[0][0];
  b = FQ_lookup[0][1];
  c = FQ_lookup[0][2];
  d = FQ_lookup[0][3];
  e = FQ_lookup[0][4];
  tag0[0] = start_gather_site(F_OFFSET(plaq_psol[a][b]), sizeof(su3_vector),
                              F2Q_d1[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_site(F_OFFSET(plaq_sol[d][e]), sizeof(su3_vector),
                              goffset[c], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_site(F_OFFSET(plaq_psol[a][b]), sizeof(su3_vector),
                              goffset[c], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_site(F_OFFSET(plaq_sol[d][e]), sizeof(su3_vector),
                              F2Q_d2[0], EVENANDODD, local_pt[0][3]);

  // Loop over lookup table
  for (counter = 0; counter < NTERMS; counter++) {
    gather = (flip + 1) % 2;
    if (counter < NTERMS - 1) {     // Start next set of gathers
      next = counter + 1;
      a = FQ_lookup[next][0];
      b = FQ_lookup[next][1];
      c = FQ_lookup[next][2];
      d = FQ_lookup[next][3];
      e = FQ_lookup[next][4];

      tag0[gather] = start_gather_site(F_OFFSET(plaq_psol[a][b]),
                      sizeof(su3_vector), F2Q_d1[next], EVENANDODD,
                      local_pt[gather][0]);
      tag1[gather] = start_gather_site(F_OFFSET(plaq_sol[d][e]),
                      sizeof(su3_vector), goffset[c], EVENANDODD,
                      local_pt[gather][1]);
      tag2[gather] = start_gather_site(F_OFFSET(plaq_psol[a][b]),
                      sizeof(su3_vector), goffset[c], EVENANDODD,
                      local_pt[gather][2]);
      tag3[gather] = start_gather_site(F_OFFSET(plaq_sol[d][e]),
                      sizeof(su3_vector), F2Q_d2[next], EVENANDODD,
                      local_pt[gather][3]);
    }

    // Do this set of computations while next set of gathers runs
    a = FQ_lookup[counter][0];
    b = FQ_lookup[counter][1];
    c = FQ_lookup[counter][2];
    d = FQ_lookup[counter][3];
    e = FQ_lookup[counter][4];
    permm = perm[a][b][c][d][e];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);
    wait_gather(tag3[flip]);
    FORALLSITES(i, s) {
      clear_su3mat_f(&(tempmat[i]));
      vec0 = (su3_vector *)(local_pt[flip][0][i]);
      vec1 = (su3_vector *)(local_pt[flip][1][i]);
      vec2 = (su3_vector *)(local_pt[flip][2][i]);
      vec3 = (su3_vector *)(local_pt[flip][3][i]);
      for (l = 0; l < NUMGEN; l++) {
        for (m = 0; m < NUMGEN; m++) {
          CMULJ_(vec0->c[l], vec1->c[m], cterm);
          CMULREAL(cterm, permm, cterm);
          BC = (s->bc2[OPP_LDIR(a)][OPP_LDIR(b)]) * (s->bc1[c]);
          CMULREAL(cterm, BC, cterm);
          c_scalar_mult_add_su3mat_f(&(tempmat[i]), &(Lambda_prod[l][m]),
                                     &cterm, &(tempmat[i]));

          CMULJ_(vec2->c[l], vec3->c[m], cterm);
          CMULREAL(cterm, -permm, cterm);
          BC = (s->bc3[a][b][c]) * (s->bc1[c]);
          CMULREAL(cterm, BC, cterm);
          c_scalar_mult_add_su3mat_f(&(tempmat[i]), &(Lambda_prod[m][l]),
                                     &cterm, &(tempmat[i]));
        }
      }
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    cleanup_gather(tag3[flip]);
    flip = (flip + 1) % 2;

    // Reproduces f_U.set(x, c, f_U.get(x, c) - 0.5 * Adj(tmp));
    FORALLSITES(i, s) {
      su3_adjoint_f(&(tempmat[i]), &tmpmat);
      scalar_mult_add_su3_matrix_f(&(s->f_U[c]), &tmpmat, -0.5, &(s->f_U[c]));
    }
  }
#endif

  // Final adjoint and minus sign
  FORALLSITES(i, s) {
    for (mu = 0; mu < NUMLINK; mu++) {
//      tmpmat = s->f_U[mu];
//      scalar_mult_su3_matrix_f(&tmpmat, -1.0, &tmpmat);
      scalar_mult_su3_matrix_f(&(s->f_U[mu]), -1.0, &tmpmat);
      su3_adjoint_f(&tmpmat, &(s->f_U[mu]));
    }
  }

  // Copy persistent site, link and plaquette fermions into Twist_Fermions
  FORALLSITES(i, s) {
    for (mu = 0; mu < NUMLINK; mu++) {
      for (nu = 0; nu < NUMLINK; nu++) {
        su3vec_copy(&(s->plaq_sol[mu][nu]), &(sol[i].Fplaq[mu][nu]));
        su3vec_copy(&(s->plaq_psol[mu][nu]), &(psol[i].Fplaq[mu][nu]));
      }
    }
    for (mu = 0; mu < NUMLINK; mu++) {
      su3vec_copy(&(s->link_sol[mu]), &(sol[i].Flink[mu]));
      su3vec_copy(&(s->link_psol[mu]), &(psol[i].Flink[mu]));
    }
    su3vec_copy(&(s->site_sol), &(sol[i].Fsite));
    su3vec_copy(&(s->site_psol), &(psol[i].Fsite));
  }
  free(tempmat);
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Update the momenta with the fermion force
// Assume that the multiCG has been run, with the answers in sol[j]
// Accumulate force into f_U and temporary fullforce
double fermion_force(Real eps, Twist_Fermion *src, Twist_Fermion **sol) {
  register int i;
  register site *s;
  int mu, n;
  double returnit = 0.0;
  su3_matrix_f **fullforce;
  Twist_Fermion *psol = malloc(sites_on_node * sizeof(*psol));

#ifdef FORCE_DEBUG
  int kick, ii, jj, iters = 0;
  Real final_rsq;
  double individ_force, old_action, new_action, diff_action;
  su3_matrix_f tmat1, tprint, tprint2;
  clear_su3mat_f(&tprint);
  clear_su3mat_f(&tmat1);
#endif

  fullforce = malloc(NUMLINK * sizeof(*fullforce));
  for (mu = 0; mu < NUMLINK; mu++) {
    fullforce[mu] = malloc(sites_on_node * sizeof(su3_matrix_f));
      FORALLSITES(i, s)
        clear_su3mat_f(&(fullforce[mu][i]));
  }

  for (n = 0; n < Norder; n++) {
    fermion_op(sol[n], psol, PLUS);
    // Makes sense to multiply here by amp4[n]...
    FORALLSITES(i, s)
      scalar_mult_TF(&(psol[i]), amp4[n], &(psol[i]));

    assemble_fermion_force(sol[n], psol);
#ifdef FORCE_DEBUG
    individ_force = 0.0;
#endif
    for (mu = 0; mu < NUMLINK; mu++) {
      FORALLSITES(i, s) {
//        scalar_mult_add_su3_matrix_f(&(s->mom[mu]), &(s->f_U[mu]), eps, &(s->mom[mu]));
        add_su3_matrix_f(&(fullforce[mu][i]), &(s->f_U[mu]), &(fullforce[mu][i]));
#ifdef FORCE_DEBUG
      if (s->x == 0 && s->y == 0 && s->z == 0 && s->t == 0 && mu == 3) {
        printf("Fermion force mu=%d on site (%d, %d, %d, %d)\n",
               mu, s->x, s->y, s->z ,s->t);
        dumpmat_f(&(s->f_U[mu]));
      }
      // Compute average gauge force
      individ_force += realtrace_su3_f(&(s->f_U[mu]), &(s->f_U[mu]));
#endif
      }
    }
#ifdef FORCE_DEBUG
    g_doublesum(&individ_force);
    node0_printf("Individ_force %d %e \n",n,eps*sqrt(individ_force)/volume);

    // Do a scan of the fermion action
    for (mu = 0; mu < NUMLINK; mu++) {
      FORALLSITES(i, s) {
        node0_printf("mu=%d on site (%d, %d, %d, %d)\n",
                     mu, s->x, s->y, s->z, s->t);
        tmat1 = s->linkf[mu];
        dumpmat_f(&(s->f_U[mu]));

        old_action = d_fermion_action(src, sol);
        node_printf("Old fermion action %.4g\n",old_action);
        for (ii = 0; ii < NCOL; ii++) {
          for (jj = 0; jj < NCOL; jj++) {
            for (kick = -1; kick <= 1; kick += 2) {
              s->linkf[mu] = tmat1;
              s->linkf[mu].e[ii][jj].real += 0.001 * (Real)kick;

              fermion_rep();
              iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq);
              if (kick== -1)
                new_action = -d_fermion_action(src, sol);
              if (kick== 1) {
                new_action += d_fermion_action(src, sol);
                tprint.e[ii][jj].real = -250.0 * new_action;
              }
            }

            for (kick = -1; kick <= 1; kick += 2) {
              s->linkf[mu]=tmat1;
              s->linkf[mu].e[ii][jj].imag += 0.001 * (Real)kick;

              fermion_rep();
              iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq, 0);
              if (kick == -1)
                new_action = -d_fermion_action(src, sol);
              if (kick == 1) {
                new_action += d_fermion_action(src, sol);
                node0_printf("XXXG%d%dI %.4g %.4g\n",
                             ii, jj, 0.001 * (Real)kick, 500 * new_action);
                tprint.e[ii][jj].imag = -250 * new_action;
              }
            }
          }
        }
        sub_su3_matrix_f(&(tprint), &(s->f_U[mu]), &tprint2);
        node0_printf("mu=%d on site (%d, %d, %d, %d): %.4g\n",
                     mu, s->x, s->y, s->z, s->t,
                     realtrace_su3_f(&(tprint2), &(tprint2)));
        dumpmat_f(&(tprint));
        s->linkf[mu] = tmat1;

        fermion_rep();
        iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq);
        new_action = d_fermion_action(src, sol);
        node0_printf("EXITING  %.4g\n", new_action - old_action);
      }
    }
#endif  // End of scan of the fermion action
  }

  // Update the momentum from the fermion force -- sum or eps
  // Opposite sign as to gauge force,
  // because dS_G / dU = 2F_g while ds_F / dU = -2F_f
  FORALLSITES(i, s) {
    for (mu = 0; mu < NUMLINK; mu++) {
      scalar_mult_add_su3_matrix_f(&(s->mom[mu]), &(fullforce[mu][i]), eps,
                                   &(s->mom[mu]));
      returnit += realtrace_su3_f(&(fullforce[mu][i]), &(fullforce[mu][i]));
    }
  }
  g_doublesum(&returnit);

  for (mu = 0; mu < NUMLINK; mu++)
    free(fullforce[mu]);
  free(fullforce);
  free(psol);
  return (eps * sqrt(returnit) / volume);
}
// -----------------------------------------------------------------
