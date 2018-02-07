// -----------------------------------------------------------------
// Update the momentum matrices
// Uncomment to print out debugging messages
//#define FORCE_DEBUG
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update mom with the bosonic force
// Uses some pointer arithmetic that might be susceptible to bugs
double bosonic_force(Real eps) {
  register int i, j, k, l;
  register site *s;
  Real tr;
  double returnit = 0.0, tmp_so3 = 2.0 * mass_so3, tmp_so6 = 2.0 * mass_so6;
  matrix tmat;
  msg_tag *tag, *tag2, *tag3;

  // First we have the finite difference operator gauge derivative
  // 2 * X_i(n + 1) * U^dagger(n) * X_i(n)
  tag = start_gather_site(F_OFFSET(X[0]), sizeof(matrix) * NSCALAR,
                          TUP, EVENANDODD, gen_pt[0]);

  tag2 = start_gather_site(F_OFFSET(X[0]), sizeof(matrix) * NSCALAR,
                           TDOWN, EVENANDODD, gen_pt[1]);

  tag3 = start_gather_site(F_OFFSET(link), sizeof(matrix),
                           TDOWN, EVENANDODD, gen_pt[2]);

  wait_gather(tag);
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      mult_an(&(s->link), &(s->X[j]), &tmat);
      mult_nn((matrix *)(gen_pt[0][i] + j), &tmat, &(s->f_U));
    }
  }

  tr = 2.0 * kappa * eps;
  FORALLSITES(i, s) {
    scalar_mult_dif_adj_matrix(&(s->f_U), tr, &(s->mom));
    returnit += realtrace(&(s->f_U), &(s->f_U));
  }

  // This is the finite difference operator scalar derivative
  wait_gather(tag2);
  wait_gather(tag3);
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      mult_na((matrix *)(gen_pt[0][i] + j), &(s->link), &tmat);
      mult_nn(&(s->link), &tmat, &(s->f_X[j]));
      dif_matrix(&(s->X[j]), &(s->f_X[j]));
      mult_nn((matrix *)(gen_pt[1][i] + j), (matrix *)(gen_pt[2][i]), &tmat);
      mult_an_sum((matrix *)(gen_pt[2][i]), &tmat, &(s->f_X[j]));
      scalar_mult_matrix(&(s->f_X[j]), 2.0, &(s->f_X[j]));
    }
  }

  // The pure scalar stuff
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
#ifdef BMN
      for (k = 0; k < NSCALAR; k++) {
        if (j == k)
          continue;
        for (l = 0; l < NSCALAR; l++) {
          if((j == l) || (k == l))
            continue;

          if (epsilon[j][k][l] > 0) {
            scalar_mult_nn_dif(&(s->X[k]), &(s->X[l]),
                               mass_Myers, &(s->f_X[j]));
          }
          else if (epsilon[j][k][l] < 0) {
            scalar_mult_nn_sum(&(s->X[k]), &(s->X[l]),
                               mass_Myers, &(s->f_X[j]));
          }
        }
      }
#endif
      for (k = j + 1; k < NSCALAR; k++) {
        mult_nn(&(s->X[j]), &(s->X[k]), &tmat);
        mult_nn_dif(&(s->X[k]), &(s->X[j]), &tmat);
        mult_nn_dif(&(s->X[k]), &tmat, &(s->f_X[j]));
        mult_nn_sum(&tmat, &(s->X[k]), &(s->f_X[j]));
      }

      mult_na((matrix *)(gen_pt[0][i] + j), &(s->link), &tmat);
      mult_nn(&(s->link), &tmat, &(s->f_X[j]));
      dif_matrix(&(s->X[j]), &(s->f_X[j]));
      mult_nn((matrix *)(gen_pt[1][i] + j), (matrix *)(gen_pt[2][i]), &tmat);
      mult_an_sum((matrix *)(gen_pt[2][i]), &tmat, &(s->f_X[j]));
    }
    for (j = 0; j < 3; j++)
      scalar_mult_dif_matrix(&(s->X[j]), tmp_so3, &(s->f_X[j]));

    for (j = 3; j < NSCALAR; j++)
      scalar_mult_dif_matrix(&(s->X[j]), tmp_so6, &(s->f_X[j]));
  }
  cleanup_gather(tag);
  cleanup_gather(tag2);
  cleanup_gather(tag3);

  // Finally take adjoint and update the momentum
  // Include overall factor of kappa = N / (2lambda)
  // Subtract to reproduce -Adj(f_X)
  // Compute average scalar force in same loop
  tr = kappa * eps;
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      scalar_mult_dif_adj_matrix(&(s->f_X[j]), tr, &(s->mom_X[j]));
      returnit += realtrace(&(s->f_X[j]), &(s->f_X[j]));
    }
  }
  g_doublesum(&returnit);
  returnit *= kappa * kappa;

  return (eps * sqrt(returnit) / nt);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Assemble fermion contributions to gauge link force,
//   f_U = Adj(Ms).D_U M(U, Ub).s - Adj[Adj(Ms).D_Ub M(U, Ub).s]
// "s" is sol while "Ms" is psol
// Copy these into persistent matrices for easier gathering
// Use tempmat, tempmat2, Tr_Uinv,
// tr_dest and Ddet[012] for temporary storage
#ifndef PUREGAUGE
void assemble_fermion_force(Twist_Fermion *sol, Twist_Fermion *psol) {
  register int i;
  register site *s;
  char **local_pt[2][2];
  int mu, nu, a, b, gather, flip = 0;
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
    FORALLDIR(mu) {
      mat_copy(&(sol[i].Flink[mu]), &(link_src[mu][i]));
      adjoint(&(psol[i].Flink[mu]), &(link_dest[mu][i]));
    }
    mat_copy(&(sol[i].Fplaq), &(plaq_src[i]));
    adjoint(&(psol[i].Fplaq), &(plaq_dest[i]));
  }

#ifdef SV
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
}
#endif
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Update the momenta with the fermion force
// Assume that the multiCG has been run, with the solution in sol[j]
// Accumulate f_U for each pole into fullforce, add to momenta
// Allocate fullforce while using temp_ferm for temporary storage
// (Calls assemble_fermion_force, which uses many more temporaries)
#ifndef PUREGAUGE
double fermion_force(Real eps, Twist_Fermion *src, Twist_Fermion **sol) {
  register int i;
  register site *s;
  int mu, n;
  double returnit = 0.0;
  matrix *fullforce = malloc(sizeof *fullforce * sites_on_node);

#ifdef FORCE_DEBUG
  int kick, ii, jj, iters = 0;
  Real final_rsq;
  double individ_force, old_action, new_action = 0.0;
  matrix tmat, tprint, tprint2;
  clear_mat(&tprint);
  clear_mat(&tmat);
#endif

  // Initialize fullforce
  fermion_op(sol[0], tempTF, PLUS);
  FORALLSITES(i, s)
    scalar_mult_TF(&(tempTF[i]), amp4[0], &(tempTF[i]));
  assemble_fermion_force(sol[0], tempTF);
  FORALLSITES(i, s)
    adjoint(&(s->f_U), &(fullforce[i]));
  for (n = 1; n < Norder; n++) {
    fermion_op(sol[n], tempTF, PLUS);
    // Makes sense to multiply here by amp4[n]...
    FORALLSITES(i, s)
      scalar_mult_TF(&(tempTF[i]), amp4[n], &(tempTF[i]));

    assemble_fermion_force(sol[n], tempTF);
#ifdef FORCE_DEBUG
    individ_force = 0.0;
#endif
    FORALLSITES(i, s) {
      // Take adjoint but don't negate yet...
      sum_adj_matrix(&(s->f_U), &(fullforce[i]));
#ifdef FORCE_DEBUG
//    if (s->t == 0 && mu == 3) {
//      printf("Fermion force mu=%d on site (%d, %d)\n", mu, s->x, s->t);
//      dumpmat(&(s->f_U[mu]));
//    }
      // Compute average gauge force
      individ_force += realtrace(&(s->f_U), &(s->f_U));
#endif
    }
#ifdef FORCE_DEBUG
    g_doublesum(&individ_force);
    node0_printf("Individ_force %d %.4g\n",
                 n, eps * sqrt(individ_force) / volume);

    // Check that force syncs with fermion action
    old_action = fermion_action(src, sol);
    iters += congrad_multi(src, sol, niter, rsqmin, &final_rsq);
    new_action = fermion_action(src, sol);
    node0_printf("EXITING  %.4g\n", new_action - old_action);
    if (fabs(new_action - old_action) > 1e-3)
      terminate(1);                             // Don't go further for now

#if 0
    // Do a scan of the fermion action
    FORALLSITES(i, s) {
      node0_printf("site %d\n", s->t);
      tmat = s->link;
      dumpmat(&(s->f_U));

      for (ii = 0; ii < NCOL; ii++) {
        for (jj = 0; jj < NCOL; jj++) {
          for (kick = -1; kick <= 1; kick += 2) {
            s->link = tmat;
            s->link.e[ii][jj].real += 0.001 * (Real)kick;

            iters += congrad_multi(src, sol, niter, rsqmin, &final_rsq);
            if (kick == -1)
              new_action -= fermion_action(src, sol);
            if (kick == 1) {
              new_action += fermion_action(src, sol);
              tprint.e[ii][jj].real = -250.0 * new_action;
            }
          }

          for (kick = -1; kick <= 1; kick += 2) {
            s->link = tmat;
            s->link.e[ii][jj].imag += 0.001 * (Real)kick;

            iters += congrad_multi(src, sol, niter, rsqmin, &final_rsq);
            if (kick == -1)
              new_action -= fermion_action(src, sol);
            if (kick == 1) {
              new_action += fermion_action(src, sol);
              node0_printf("XXXG%d%dI %.4g %.4g\n",
                           ii, jj, 0.001 * (Real)kick, 500 * new_action);
              tprint.e[ii][jj].imag = -250 * new_action;
            }
          }
        }
      }
      sub_matrix(&tprint, &(s->f_U[mu]), &tprint2);
      node0_printf("site %d: %.4g\n", s->t,
                   realtrace(&tprint2, &tprint2));
      dumpmat(&tprint);
      s->link[mu] = tmat;

      iters += congrad_multi(src, sol, niter, rsqmin, &final_rsq);
    }   // End scan of the fermion action
#endif
#endif
  }

  // Update the momentum from the fermion force -- sum or eps
  // Opposite sign as to gauge force,
  // because dS_G / dU = 2F_g while ds_F / dU = -2F_f
  // Move negation here as well, though adjoint remains above
  FORALLSITES(i, s) {
    scalar_mult_dif_matrix(&(fullforce[i]), eps, &(s->mom));
    returnit += realtrace(&(fullforce[i]), &(fullforce[i]));
  }
  g_doublesum(&returnit);

  free(fullforce);
  return (eps * sqrt(returnit) / volume);
}
#endif
// -----------------------------------------------------------------
