// -----------------------------------------------------------------
// Update the momentum matrices
// Uncomment to print out debugging messages
//#define FORCE_DEBUG
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update mom with the bosonic force
double bosonic_force(Real eps) {
  register int i, j, k, l, pt, pt2 = 2 * NSCALAR + 1;
  register site *s;
  Real tr;
  double returnit = 0.0, tmp_so3 = 2.0 * mass_so3, tmp_so6 = 2.0 * mass_so6;
  matrix tmat, tmat2;
  msg_tag *tag[NSCALAR], *tag2[NSCALAR], *tag3;
#ifdef DEBUG_CHECK
  anti_hermitmat tah;
#endif

  // Clear the force collectors
  FORALLSITES(i, s) {
    clear_mat(&(s->f_U));
    for (j = 0; j < NSCALAR; j++)
      clear_mat(&(s->f_X[j]));
  }

#ifdef KINETIC
  // First we have the finite difference operator gauge derivative
  // Must transform as site variable so momenta can be exponentiated
  //   U(n) d/dU(n) Tr[2 U(t) X(t+1) Udag(t) X(t)  - X(t+1) X(t+1) - X(t) X(t)]
  //     = 2 delta_{nt} U(n) X(t+1) Udag(t) X(t) = 2 U(n) X(n+1) Udag(n) X(n)
  for (j = 0; j < NSCALAR; j++) {
    tag[j] = start_gather_site(F_OFFSET(X[j]), sizeof(matrix),
                               TUP, EVENANDODD, gen_pt[j]);

    tag2[j] = start_gather_site(F_OFFSET(X[j]), sizeof(matrix),
                                TDOWN, EVENANDODD, gen_pt[NSCALAR + j]);
  }
  tag3 = start_gather_site(F_OFFSET(link), sizeof(matrix),
                           TDOWN, EVENANDODD, gen_pt[pt2]);

  for (j = 0; j < NSCALAR; j++)
   wait_gather(tag[j]);

  FORALLSITES(i, s) {
    clear_mat(&(s->f_U));   // Clear the force collectors
    for (j = 0; j < NSCALAR; j++) {   // X(n+1) = gen_pt[j]
      mult_na((matrix *)(gen_pt[j][i]), &(s->link), &tmat);
      mult_nn(&(s->link), &tmat, &tmat2);
      mult_nn_sum(&tmat2, &(s->X[j]), &(s->f_U));
    }
  }
#endif

  // Take adjoint and update the gauge momenta
  // Make them anti-hermitian following non-susy code
  // Include overall factor of kappa = N / (4lambda), and factor of 2
  // Subtract to reproduce -Adj(f_U)
  // Compute average gauge force in same loop
  tr = 2.0 * kappa * eps;
  FORALLSITES(i, s) {
    uncompress_anti_hermitian(&(s->mom), &tmat);
    scalar_mult_sum_matrix(&(s->f_U), tr, &tmat);
    make_anti_hermitian(&tmat, &(s->mom));
    returnit += realtrace(&(s->f_U), &(s->f_U));
  }

#ifdef KINETIC
  // This is the finite difference operator scalar derivative
  //   d/dX(n) Tr[X(t) U(t) X(t+1) Udag(t) + X(t+1) Udag(t) X(t) U(t)
  //              - X(t+1) X(t+1) - X(t) X(t)]
  //     = 2 delta_{nt} U(t) X(t+1) Udag(t)
  //       + 2 delta_{n(t+1)} Udag(t) X(t) U(t)
  //       - 2 delta_{n(t+1)} X(t+1) - 2 delta_{nt} X(t)
  //     = 2 [U(n) X(n+1) Udag(n) + Udag(n-1) X(n-1) U(n-1) - 2 X(n)]
  for (j = 0; j < NSCALAR; j++)
    wait_gather(tag2[j]);
  wait_gather(tag3);
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      pt = NSCALAR + j;             // X(n-1) is gen_pt[pt]

      // Initialize force with on-site -2X(n)
      scalar_mult_matrix(&(s->X[j]), -2.0, &(s->f_X[j]));

      // Add forward hopping term using X(n+1) = gen_pt[j]
      mult_na((matrix *)(gen_pt[j][i]), &(s->link), &tmat);
      mult_nn_sum(&(s->link), &tmat, &(s->f_X[j]));

      // Add backward hopping term using X(n-1) = gen_pt[pt]
      //                             and U(n-1) = gen_pt[pt2]
      mult_nn((matrix *)(gen_pt[pt][i]), (matrix *)(gen_pt[pt2][i]), &tmat);
      mult_an_sum((matrix *)(gen_pt[pt2][i]), &tmat, &(s->f_X[j]));

      scalar_mult_matrix(&(s->f_X[j]), 2.0, &(s->f_X[j]));
    }
  }
#endif

  // The pure scalar stuff
  tr = 3.0 * mass_Myers;
  FORALLSITES(i, s) {
#ifdef BMN
#ifdef MYERS
    // Myers term
    //   d/dX_i(n) -Tr[eps_{jkl} X_j(t) X_k(t) X_l(t)]]
    //     = -eps_{jkl} [delta_{ij} X_k(n) X_l(n) + delta_{ik} X_l(n) X_j(n)
    //                                            + delta_{il} X_j(n) X_k(n)]
    //     = -eps_{ikl} X_k(n) X_l(n) - eps_{ilj} X_l(n) X_j(n)
    //                                - eps_{ijk} X_j(n) X_k(n)
    //     = -3eps_{ikl} X_k(n) X_l(n)
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        if (j == k)
          continue;
        for (l = 0; l < 3; l++) {
          if ((j == l) || (k == l))
            continue;

          if (epsilon[j][k][l] > 0)     // Overall negative sign absorbed
            scalar_mult_nn_dif(&(s->X[k]), &(s->X[l]), tr, &(s->f_X[j]));
          else if (epsilon[j][k][l] < 0)
            scalar_mult_nn_sum(&(s->X[k]), &(s->X[l]), tr, &(s->f_X[j]));
        }
      }
    }
#endif
#endif

#ifdef COMM_TERM
    // Commutator term (usual cyclic product rule trick...)
    //   d/dX_k(n) sum_{i != j} -Tr[  X_i(t) X_j(t) X_i(t) X_j(t)
    //                              + X_j(t) X_i(t) X_j(t) X_i(t)
    //                              + X_i(t) X_j(t) X_i(t) X_j(t)
    //                              + X_j(t) X_i(t) X_j(t) X_i(t)
    //                              - X_i(t) X_j(t) X_j(t) X_i(t)
    //                              - X_j(t) X_j(t) X_i(t) X_i(t)
    //                              - X_j(t) X_i(t) X_i(t) X_j(t)
    //                              - X_i(t) X_i(t) X_j(t) X_j(t)]
    //     = sum_{j != k} -[4X_j(n) X_k(n) X_j(n) - 2X_j(n) X_j(n) X_k(n)
    //                                            - 2X_k(n) X_j(n) X_j(n)]
    //     = sum_{j != k} -2[X_j(n) X_k(n) X_j(n) - X_j(n) X_j(n) X_k(n)
    //                     + X_j(n) X_k(n) X_j(n) - X_k(n) X_j(n) X_j(n)]
    //     = sum_{j != k} -2[X_j(n) [X_k(n), X_j(n)] + [X_j(n), X_k(n)] X_j(n)]
    //     = sum_{j != k} -2[X_j(n) [X_k(n), X_j(n)] - [X_k(n), X_j(n)] X_j(n)]
    //     = sum_{j != k} -2[X_j(n), [X_k(n), X_j(n)]]
    //   --> sum_{k != j} -2[X_k(n), [X_j(n), X_k(n)]], j fixed
    // Following serial code, keep k>j with no factor 2
    for (j = 0; j < NSCALAR; j++) {
      for (k = j + 1; k < NSCALAR; k++) {
        if (k == j)
          continue;

        // tmat = 2[X_j, X_k]
        mult_nn(&(s->X[j]), &(s->X[k]), &tmat);
        mult_nn_dif(&(s->X[k]), &(s->X[j]), &tmat);
//        scalar_mult_matrix(&tmat, 2.0, &tmat);
        // Overall negative sign absorbed below
        mult_nn_dif(&(s->X[k]), &tmat, &(s->f_X[j]));
        mult_nn_sum(&tmat, &(s->X[k]), &(s->f_X[j]));
      }
    }
#endif

#ifdef SCALAR_POT
    // Simple d/dX_i(n) -X_j(t)^2 = -2 X_i(n)
    // Absorb factor of two into tmp_so# = 2 * mass_so#
    // Coefficients depend on BMN vs. BFSS, set in setup.c
    for (j = 0; j < 3; j++)
      scalar_mult_dif_matrix(&(s->X[j]), tmp_so3, &(s->f_X[j]));

    for (j = 3; j < NSCALAR; j++)
      scalar_mult_dif_matrix(&(s->X[j]), tmp_so6, &(s->f_X[j]));
#endif
  }
  for (j = 0; j < NSCALAR; j++) {
    cleanup_gather(tag[j]);
    cleanup_gather(tag2[j]);
  }
  cleanup_gather(tag3);

  // Take adjoint and update the scalar momenta
  // Include overall factor of kappa = N / (4lambda)
  // Subtract to reproduce -Adj(f_X)
  // Compute average scalar force in same loop (combine with gauge from above)
  tr = kappa * eps;
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
#ifdef DEBUG_CHECK
      // Make f_X traceless anti-hermitian, which it should be already
      make_anti_hermitian(&(s->f_X[j]), &tah);
      uncompress_anti_hermitian(&tah, &(s->f_X[j]));
#endif
      scalar_mult_sum_matrix(&(s->f_X[j]), tr, &(s->mom_X[j]));
      returnit += realtrace(&(s->f_X[j]), &(s->f_X[j]));
    }
  }
  g_doublesum(&returnit);
  returnit *= kappa * kappa;

  return (eps * sqrt(returnit) / (double)nt);
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
