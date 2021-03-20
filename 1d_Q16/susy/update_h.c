// -----------------------------------------------------------------
// Update the momentum matrices
// Uncomment to print out debugging messages
//#define FORCE_DEBUG
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update mom with the bosonic force
double bosonic_force(Real eps) {
  register int i, j, k, l;
  register site *s;
  Real tr;
  double returnit = 0.0, tmp_so3 = 2.0 * mass_so3, tmp_so6 = 2.0 * mass_so6;
  matrix tmat;
#ifndef UNGAUGED
  matrix tmat2;
#endif
  msg_tag *tag[NSCALAR], *tag2[NSCALAR];
#ifdef DEBUG_CHECK
  anti_hermitmat tah;
#endif

  // Clear the force collectors
  FORALLSITES(i, s) {
#ifndef UNGAUGED
    clear_mat(&(s->f_U));
#endif
    for (j = 0; j < NSCALAR; j++)
      clear_mat(&(s->f_X[j]));
  }

  // First we have the finite difference operator gauge derivative
  // Must transform as site variable so momenta can be exponentiated
  //   U(n) d/dU(n) Tr[2 U(t) X(t+1) Udag(t) X(t) - X(t+1) X(t+1) - X(t) X(t)]
  //     = 2 delta_{nt} U(n) X(t+1) Udag(t) X(t) = 2 U(n) X(n+1) Udag(n) X(n)
  // Also used in scalar force, so gather even if ungauged
  for (j = 0; j < NSCALAR; j++) {
    tag[j] = start_gather_site(F_OFFSET(X[j]), sizeof(matrix),
                               TUP, EVENANDODD, gen_pt[j]);
  }

  for (j = 0; j < NSCALAR; j++) {
    // For scalar force term, compute and gather Udag(n-1) X(n-1) U(n-1)
    FORALLSITES(i, s) {
#ifndef UNGAUGED
      mult_nn(&(s->X[j]), &(s->link), &tmat);
      mult_an(&(s->link), &tmat, &(temp_X[j][i]));
#else
      mat_copy(&(s->X[j]), &(temp_X[j][i]));
#endif
    }
    tag2[j] = start_gather_field(temp_X[j], sizeof(matrix),
                                 TDOWN, EVENANDODD, gen_pt[NSCALAR + j]);
  }

  for (j = 0; j < NSCALAR; j++) {   // X(n+1) = gen_pt[j]
    wait_gather(tag[j]);
#ifndef UNGAUGED
    FORALLSITES(i, s) {
      mult_na((matrix *)(gen_pt[j][i]), &(s->link), &tmat);
      mult_nn(&(s->link), &tmat, &tmat2);
      mult_nn_sum(&(s->X[j]), &tmat2, &(s->f_U));
    }
#endif
  }

#ifndef UNGAUGED
  // Take adjoint and update the gauge momenta
  // Make them anti-hermitian following non-susy code
  // Include overall factor of kappa = N / (4lambda), and factor of 2
  // !!! Another factor of 2 needed for conservation (real vs. complex?)
  // Compute average gauge force in same loop
  tr = 4.0 * kappa * eps;
  FORALLSITES(i, s) {
    uncompress_anti_hermitian(&(s->mom), &tmat);
    scalar_mult_dif_matrix(&(s->f_U), tr, &tmat);
    make_anti_hermitian(&tmat, &(s->mom));
    returnit += realtrace(&(s->f_U), &(s->f_U));
  }
#endif

  // This is the finite difference operator scalar derivative
  //   d/dX(n) Tr[X(t) U(t) X(t+1) Udag(t) + X(t+1) Udag(t) X(t) U(t)
  //              - X(t+1) X(t+1) - X(t) X(t)]
  //     = 2 delta_{nt} U(t) X(t+1) Udag(t)
  //       + 2 delta_{n(t+1)} Udag(t) X(t) U(t)
  //       - 2 delta_{n(t+1)} X(t+1) - 2 delta_{nt} X(t)
  //     = 2 [U(n) X(n+1) Udag(n) + Udag(n-1) X(n-1) U(n-1) - 2 X(n)]
  for (j = 0; j < NSCALAR; j++) {
    wait_gather(tag2[j]);
    FORALLSITES(i, s) {
      // Initialize force with on-site -2X(n)
      scalar_mult_matrix(&(s->X[j]), -2.0, &(s->f_X[j]));

      // Add forward hopping term using X(n+1) = gen_pt[j]
#ifndef UNGAUGED
      mult_na((matrix *)(gen_pt[j][i]), &(s->link), &tmat);
      mult_nn_sum(&(s->link), &tmat, &(s->f_X[j]));
#else
      sum_matrix((matrix *)(gen_pt[j][i]), &(s->f_X[j]));
#endif

      // Add backward hopping term
      //   Udag(n-1) X(n-1) U(n-1) = gen_pt[NSCALAR + j]
      sum_matrix((matrix *)(gen_pt[NSCALAR + j][i]), &(s->f_X[j]));

      scalar_mult_matrix(&(s->f_X[j]), 2.0, &(s->f_X[j]));
    }
    cleanup_gather(tag[j]);
    cleanup_gather(tag2[j]);
  }

  // The pure scalar stuff
  tr = 3.0 * mass_Myers;
  FORALLSITES(i, s) {
#ifdef BMN
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

          if (epsilon[j][k][l] > 0)     // Following 02 Jun 2019 notes
            scalar_mult_nn_sum(&(s->X[k]), &(s->X[l]), tr, &(s->f_X[j]));
          else if (epsilon[j][k][l] < 0)
            scalar_mult_nn_dif(&(s->X[k]), &(s->X[l]), tr, &(s->f_X[j]));
        }
      }
    }
#endif

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
    // Overall factor of 2 cancels with 1/2 of Eq. 14 of 02 Jun 2019 notes
    for (j = 0; j < NSCALAR; j++) {
      for (k = 0; k < NSCALAR; k++) {
        if (k == j)
          continue;

        // tmat = [X_j, X_k]
        mult_nn(&(s->X[j]), &(s->X[k]), &tmat);
        mult_nn_dif(&(s->X[k]), &(s->X[j]), &tmat);
        // Overall negative sign absorbed below
        mult_nn_dif(&(s->X[k]), &tmat, &(s->f_X[j]));
        mult_nn_sum(&tmat, &(s->X[k]), &(s->f_X[j]));
      }
    }

    // Simple d/dX_i(n) -X_j(t)^2 = -2 X_i(n)
    // Absorb factor of two into tmp_so# = 2 * mass_so#
    // Coefficients depend on BMN vs. BFSS, set in setup.c
    for (j = 0; j < 3; j++)
      scalar_mult_dif_matrix(&(s->X[j]), tmp_so3, &(s->f_X[j]));

    for (j = 3; j < NSCALAR; j++)
      scalar_mult_dif_matrix(&(s->X[j]), tmp_so6, &(s->f_X[j]));
  }

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
// Take Adj(Ms) for easier gathering
// Accumulate contribution into f_U and f_X[k], no overwriting
// Serial code has a factor of 2 * 0.5 = 1.0 which we ignore
#ifndef PUREGAUGE
void assemble_fermion_force(matrix **sol, matrix *psol[NFERMION]) {
  register int i, j, k, l, m, n;
  register site *s;
  Real tr;
  // sqrt(2) factor from Eq. 14 of 2 Jun 2019 notes
  complex tc = cmplx(0.0, -1.0 / sqrt(2.0));
  matrix tmat;
#ifndef UNGAUGED
  matrix tmat2;
  msg_tag *tag[NCHIRAL_FERMION], *tag2[NCHIRAL_FERMION];
#endif

  // First gauge force s->f_U
  // For gathering it is convenient to overwrite psol by its adjoint
  for (j = 0; j < NFERMION; j++) {
    FORALLSITES(i, s) {
      adjoint(&(psol[j][i]), &tmat);
      mat_copy(&tmat, &(psol[j][i]));
    }
  }
#ifndef UNGAUGED
  for (j = 0; j < NCHIRAL_FERMION; j++) {
    k = j + NCHIRAL_FERMION;
    tag[j] = start_gather_field(sol[k], sizeof(matrix),
                                TUP, EVENANDODD, gen_pt[j]);
    tag2[j] = start_gather_field(psol[k], sizeof(matrix),
                                 TUP, EVENANDODD, gen_pt[k]);
  }

  for (j = 0; j < NCHIRAL_FERMION; j++) {
    k = j + NCHIRAL_FERMION;
    wait_gather(tag[j]);
    FORALLSITES(i, s) {
      mult_nn(&(s->link), (matrix *)(gen_pt[j][i]), &tmat);
      scalar_mult_na(&tmat, &(s->link), s->bc[0], &tmat2);
      mult_nn_dif(&tmat2, &(psol[j][i]), &(s->f_U));
      mult_nn_sum(&(psol[j][i]), &tmat2, &(s->f_U));
    }
    cleanup_gather(tag[j]);

    wait_gather(tag2[j]);
    FORALLSITES(i, s) {
      mult_nn(&(s->link), (matrix *)(gen_pt[k][i]), &tmat);
      scalar_mult_na(&tmat, &(s->link), s->bc[0], &tmat2);
      mult_nn_sum(&tmat2, &(sol[j][i]), &(s->f_U));
      mult_nn_dif(&(sol[j][i]), &tmat2, &(s->f_U));
    }
    cleanup_gather(tag2[j]);
  }
#endif

  // Now scalar forces s->f_X[k] from Yukawa terms
  // Sqrt factor from Eq. 14 in 2 Jun 2019 notes
  // Gamma[0--6] embedded into -i*sigma_2 = ( 0 -1 \\ 1 0 ) in 8x8 blocks
  for (j = 0; j < NCHIRAL_FERMION; j++) {
    m = j + NCHIRAL_FERMION;
    for (k = 0; k < NCHIRAL_FERMION; k++) {
      n = k + NCHIRAL_FERMION;
      for (l = 0; l < NSCALAR - 2; l++) {
        if (Gamma[l][j][k] != 0) {
          tr = (Real)Gamma[l][j][k] / sqrt(2.0);
          FORALLSITES(i, s) {
            mult_nn(&(psol[j][i]), &(sol[n][i]), &tmat);
            mult_nn_dif(&(sol[n][i]), &(psol[j][i]), &tmat);
            scalar_mult_dif_matrix(&tmat, tr, &(s->f_X[l]));

            mult_nn(&(psol[m][i]), &(sol[k][i]), &tmat);
            mult_nn_dif(&(sol[k][i]), &(psol[m][i]), &tmat);
            scalar_mult_sum_matrix(&tmat, tr, &(s->f_X[l]));
          }
        }
      }
    }

    // Last 2 gammas are diagonal
    k = NSCALAR - 2;
    n = NSCALAR - 1;
    tr = 1.0 / sqrt(2.0);                 // Sqrt from 2 Jun 2019 notes
    FORALLSITES(i, s) {
      // Gamma[7] = ( -1 0 \\ 0 1 ) in 8x8 blocks
      scalar_mult_nn_dif(&(psol[j][i]), &(sol[j][i]), tr, &(s->f_X[k]));
      scalar_mult_nn_sum(&(sol[j][i]), &(psol[j][i]), tr, &(s->f_X[k]));

      scalar_mult_nn_sum(&(psol[m][i]), &(sol[m][i]), tr, &(s->f_X[k]));
      scalar_mult_nn_dif(&(sol[m][i]), &(psol[m][i]), tr, &(s->f_X[k]));

      // Gamma[8] = ( -i 0 \\ 0 -i ) in 8x8 blocks
      // Negative absorbed into tc above
      mult_nn(&(psol[j][i]), &(sol[j][i]), &tmat);
      mult_nn_dif(&(sol[j][i]), &(psol[j][i]), &tmat);

      mult_nn_sum(&(psol[m][i]), &(sol[m][i]), &tmat);
      mult_nn_dif(&(sol[m][i]), &(psol[m][i]), &tmat);

      c_scalar_mult_sum_mat(&tmat, &tc, &(s->f_X[n]));
    }
  }
}
#endif
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Update the momenta with the fermion force
// Assume that the multiCG has been run, with the solution in sol[j]
// Use temp_ferm for temporary storage
// assemble_fermion_force above accumulates into f_U and f_X
// Calls fermion_op, which uses tempmat
#ifndef PUREGAUGE
double fermion_force(Real eps, matrix **src, matrix ***sol) {
  register int i, k, n;
  register site *s;
  double returnit = 0.0;
#ifndef UNGAUGED
  matrix tmat;
#endif
  anti_hermitmat tah;

  // Clear the force collectors
  FORALLSITES(i, s) {
#ifndef UNGAUGED
    clear_mat(&(s->f_U));
#endif
    for (k = 0; k < NSCALAR; k++)
      clear_mat(&(s->f_X[k]));
  }

  // Accumulate forces from each pole in the force collectors
  for (n = 0; n < Norder; n++) {
    fermion_op(sol[n], temp_ferm, PLUS);
    // Makes sense to multiply here by amp4[n]...
    for (k = 0; k < NFERMION; k++) {
      FORALLSITES(i, s)
        scalar_mult_matrix(&(temp_ferm[k][i]), amp4[n], &(temp_ferm[k][i]));
    }
    assemble_fermion_force(sol[n], temp_ferm);
  }

  // Make sure forces are traceless anti-hermitian
  FORALLSITES(i, s) {
#ifndef UNGAUGED
    make_anti_hermitian(&(s->f_U), &tah);
    uncompress_anti_hermitian(&tah, &(s->f_U));
#endif
    for (k = 0; k < NSCALAR; k++) {
      make_anti_hermitian(&(s->f_X[k]), &tah);
      uncompress_anti_hermitian(&tah, &(s->f_X[k]));
    }
  }

  // Update the momentum from the fermion force -- sum or eps
  // Opposite sign as to gauge force,
  // because dS_G / dU = 2F_g while ds_F / dU = -2F_f
  // Move negation here as well, though adjoint remains above
  FORALLSITES(i, s) {
#ifndef UNGAUGED
    uncompress_anti_hermitian(&(s->mom), &tmat);
    scalar_mult_sum_matrix(&(s->f_U), eps, &tmat);
    make_anti_hermitian(&tmat, &(s->mom));
    returnit += realtrace(&(s->f_U), &(s->f_U));
#endif
    for (k = 0; k < NSCALAR; k++) {
      scalar_mult_sum_matrix(&(s->f_X[k]), eps, &(s->mom_X[k]));
      returnit += realtrace(&(s->f_X[k]), &(s->f_X[k]));
    }
  }
  g_doublesum(&returnit);
  return (eps * sqrt(returnit) / nt);
}
#endif
// -----------------------------------------------------------------
