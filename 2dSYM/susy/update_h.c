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
  register Real eb3, eb3UdU;
  double returnit = 0;
  complex trUdU;
  msg_tag *tag[NUMLINK], *tag0, *tag1;
  su3_matrix_f tmat1, tmat2, tmat3, UdU;

  // Contribution from (D_a U_a)^2 / 2 term
  compute_DmuUmu();
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

  // Factor of kappa = N Nt^2 / (2lambda) on both (D_a U_a)^2 and F^2 terms
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

  // Finally take adjoint and update the momentum
  // Subtract to reproduce -Adj(f_U)
  for (mu = 0; mu < NUMLINK; mu++) {
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
  int mu, nu, a, b;
  complex cterm;
  msg_tag *mtag[NUMLINK];
  msg_tag *mtag0 = NULL, *mtag1 = NULL;
  su3_vector *vec, tvec;
  su3_matrix_f tmpmat, tmpmat2;
  su3_matrix_f *tempmat = malloc(sites_on_node * sizeof(*tempmat));

  // Copy Twist_Fermions into persistent site, link and plaquette fermions
  FORALLSITES(i, s) {
    su3vec_copy(&(sol[i].Fplaq), &(s->plaq_sol));
    su3vec_copy(&(psol[i].Fplaq), &(s->plaq_psol));
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
        for (a = 0; a < DIMF; a++) {
          if (mu > nu) {    // plaq_psol is anti-symmetric under mu <--> nu
            CNEGATE((s->plaq_psol).c[a], tvec.c[a]);
          }                 // Suppress compiler error
          else
            tvec.c[a] = (s->plaq_psol).c[a];
          for (b = 0; b < DIMF; b++) {
            CMULJ_(tvec.c[a], vec->c[b], cterm);
            CNEGATE(cterm, cterm);
            CMULREAL(cterm, s->bc[mu], cterm);
            c_scalar_mult_add_su3mat_f(&(s->f_U[mu]), &(Lambda_prod[b][a]),
                                       &cterm, &(s->f_U[mu]));
          }
        }
      }
      cleanup_gather(mtag0);

      FORALLSITES(i, s) {
        clear_su3mat_f(&(tempmat[i]));
        for (a = 0; a < DIMF; a++) {
          if (mu > nu) {    // plaq_psol is anti-symmetric under mu <--> nu
            CNEGATE((s->plaq_psol).c[a], tvec.c[a]);
          }                 // Suppress compiler error
          else
            tvec.c[a] = (s->plaq_psol).c[a];
          for (b = 0; b < DIMF; b++) {
            CMULJ_(tvec.c[a], (s->link_sol[nu]).c[b], cterm);
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
        for (b = 0; b < DIMF; b++) {
          if (mu > nu) {    // plaq_sol is anti-symmetric under mu <--> nu
            CNEGATE((s->plaq_sol).c[b], tvec.c[b]);
          }                 // Suppress compiler error
          else
            tvec.c[b] = (s->plaq_sol).c[b];
        }
        for (a = 0; a < DIMF; a++) {
          for (b = 0; b < DIMF; b++) {
            CMULJ_((s->link_psol[nu]).c[a], tvec.c[b], cterm);
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
        for (b = 0; b < DIMF; b++) {
          if (mu > nu) {    // plaq_sol is anti-symmetric under mu <--> nu
            CNEGATE((s->plaq_sol).c[b], tvec.c[b]);
          }                 // Suppress compiler error
          else
            tvec.c[b] = (s->plaq_sol).c[b];
        }
        for (a = 0; a < DIMF; a++) {
          for (b = 0; b < DIMF; b++) {
            CMULJ_(vec->c[a], tvec.c[b], cterm);
            CMULREAL(cterm, s->bc[mu], cterm);
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
      for (a = 0; a < DIMF; a++) {
        for (b = 0; b < DIMF; b++) {
          CMULJ_((s->site_psol).c[a],
                 (s->link_sol[mu]).c[b], cterm);
          CNEGATE(cterm, cterm);
          c_scalar_mult_add_su3mat_f(&tmpmat, &(Lambda_prod[a][b]),
                                     &cterm, &tmpmat);

          CMULJ_((vec->c[a]), (s->link_sol[mu]).c[b], cterm);
          CMULREAL(cterm, s->bc[mu], cterm);
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
      for (a = 0; a < DIMF; a++) {
        for (b = 0; b < DIMF; b++) {
          CMULJ_((s->link_psol[mu]).c[a], vec->c[b], cterm);
          CMULREAL(cterm, s->bc[mu], cterm);
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

  // Final adjoint and minus sign
  FORALLSITES(i, s) {
    for (mu = 0; mu < NUMLINK; mu++) {
      scalar_mult_su3_matrix_f(&(s->f_U[mu]), -1.0, &tmpmat);
      su3_adjoint_f(&tmpmat, &(s->f_U[mu]));
    }
  }

  // Copy persistent site, link and plaquette fermions into Twist_Fermions
  FORALLSITES(i, s) {
    su3vec_copy(&(s->plaq_sol), &(sol[i].Fplaq));
    su3vec_copy(&(s->plaq_psol), &(psol[i].Fplaq));
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
    node0_printf("Individ_force %d %.4g\n",
                 n, eps * sqrt(individ_force) / volume);

    // Do a scan of the fermion action
    for (mu = 0; mu < NUMLINK; mu++) {
      FORALLSITES(i, s) {
        node0_printf("mu=%d on site (%d, %d)\n", mu, s->x, s->t);
        tmat1 = s->linkf[mu];
        dumpmat_f(&(s->f_U[mu]));

        old_action = d_fermion_action(src, sol);
        node_printf("Old fermion action %.4g\n", old_action);
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
        sub_su3_matrix_f(&tprint, &(s->f_U[mu]), &tprint2);
        node0_printf("mu=%d on site (%d, %d): %.4g\n", mu, s->x, s->t,
                     realtrace_su3_f(&tprint2, &tprint2));
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
