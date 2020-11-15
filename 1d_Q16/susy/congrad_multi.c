// -----------------------------------------------------------------
// Multi-mass conjugate gradient algorithm a la B. Jegerlehner
// The number of masses is runtime input Norder
// shift[Norder] is the array of mass values, set up in setup_rhmc.c
// BEWARE: The temporary TF pm[0] is never malloced or used

// At least for now we hard-code a zero initial guess
// We check all psi for convergence and quit doing the converged ones

// In this version of the code, all the scalars are real because M = Ddag D
//#define CG_DEBUG
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Keep all Norder-dependent mallocs here so that we can change Norder
// Return number of iterations
// src is where the source is created
// psim[Norder] are working TFs for the conjugate gradient
// MaxCG is the maximum number of iterations per solve
// errormin is the target |r|^2, scaled below by source_norm = |src|^2
// size_r is the final obtained |r|^2, hopefully < errormin * source_norm
#ifndef PUREGAUGE
int congrad_multi(matrix **src, matrix ***psim,
                  int MaxCG, Real errormin, Real *size_r) {

  register int i, j, k;
  register site *s;
  int N_iter, iteration = 0;
  int *converged = malloc(sizeof *converged * Norder);
  double rsq, rsqnew, source_norm = 0.0, rsqstop, c1, c2, cd;
  double *zeta_i   = malloc(sizeof *zeta_i * Norder);
  double *zeta_im1 = malloc(sizeof *zeta_im1 * Norder);
  double *zeta_ip1 = malloc(sizeof *zeta_ip1 * Norder);
  double *beta_i   = malloc(sizeof *beta_i * Norder);
  double *beta_im1 = malloc(sizeof *beta_im1 * Norder);
  double *alpha    = malloc(sizeof *alpha * Norder);
  double rsqj;
  matrix ***pm = malloc(sizeof(matrix**) * Norder);

  for (j = 1; j < Norder; j++) { // !!!
    pm[j] = malloc(sizeof(matrix*) * NFERMION);
    for (i = 0; i < NFERMION; i++)
      pm[j][i] = malloc(sizeof(matrix) * sites_on_node);
  }

  // Initialize zero initial guess, etc.
  // dest = 0, r = source, pm[j] = r
  for (i = 0; i < Norder; i++)
    converged[i] = 0;
  for (k = 0; k < NFERMION; k++) {
    FORALLSITES(i, s) {
      mat_copy(&(src[k][i]), &(rm[k][i]));
      mat_copy(&(rm[k][i]), &(pm0[k][i]));
      clear_mat(&(psim[0][k][i]));
      for (j = 1; j < Norder; j++) {
        clear_mat(&(psim[j][k][i]));
        mat_copy(&(rm[k][i]), &(pm[j][k][i]));
      }
    }
  }

  for (k = 0; k < NFERMION; k++) {
    FORALLSITES(i, s)
      source_norm += (double)realtrace(&(src[k][i]), &(src[k][i]));
  }
  g_doublesum(&source_norm);
  rsq = source_norm;
  rsqstop = errormin * source_norm;
#ifdef CG_DEBUG
  node0_printf("congrad: source_norm = %.4g\n", source_norm);
  node0_printf("stopping when residue is %.4g\n", rsqstop);
#endif

  for (j = 0; j < Norder; j++) {
    zeta_im1[j] = 1;
    zeta_i[j] = 1;
    alpha[j] = 0;
    beta_im1[j] = 1;
  }

  for (N_iter = 0; N_iter < MaxCG && rsq > rsqstop; N_iter++) {
    // mp = (M(u) + shift[0]) pm
    DSq(pm0, mpm);
    iteration++;
    total_iters++;
    for (k = 0; k < NFERMION; k++) {
      FORALLSITES(i, s)
        scalar_mult_sum_matrix(&(pm0[k][i]), shift[0], &(mpm[k][i]));
    }
    // beta_i[0] = -(r, r) / (pm, Mpm)
    cd = 0;
    for (k = 0; k < NFERMION; k++) {
      FORALLSITES(i, s)
        cd += realtrace(&(pm0[k][i]), &(mpm[k][i]));
    }
    g_doublesum(&cd);

    beta_i[0] = -rsq / cd;
#ifdef CG_DEBUG
    node0_printf("beta_i %.4g rsq %.4g cd %.4g\n", beta_i[0], rsq, cd);
#endif

    // beta_i(sigma)
    // zeta_ip1(sigma)
    zeta_ip1[0] = 1;
    for (j = 1; j < Norder; j++) {
      if (converged[j] == 0) {
        zeta_ip1[j] = zeta_i[j] * zeta_im1[j] * beta_im1[0];
        c1 = beta_i[0] * alpha[0] * (zeta_im1[j] - zeta_i[j]);
        c2 = zeta_im1[j] * beta_im1[0] * (1 - (shift[j] - shift[0]) * beta_i[0]);
        zeta_ip1[j] /= c1 + c2;
        beta_i[j] = beta_i[0] * zeta_ip1[j] / zeta_i[j];
      }
    }

    // psim[j] = psim[j] - beta[j] * pm[j]
    for (k = 0; k < NFERMION; k++) {
      FORALLSITES(i, s) {
        scalar_mult_dif_matrix(&(pm0[k][i]), (Real)beta_i[0],
                               &(psim[0][k][i]));
        for (j = 1; j < Norder; j++) {
          if (converged[j] == 0) {
            scalar_mult_dif_matrix(&(pm[j][k][i]), (Real)beta_i[j],
                                   &(psim[j][k][i]));
          }
        }
      }
    }

    // r = r + beta[0] * mp
    for (k = 0; k < NFERMION; k++) {
      FORALLSITES(i, s)
        scalar_mult_sum_matrix(&(mpm[k][i]), (Real)beta_i[0], &(rm[k][i]));
    }
    // alpha_ip1[j]
    rsqnew = 0;
    for (k = 0; k < NFERMION; k++) {
      FORALLSITES(i, s)
        rsqnew += (double)realtrace(&(rm[k][i]), &(rm[k][i]));
    }
    g_doublesum(&rsqnew);
    alpha[0] = rsqnew / rsq;
#ifdef CG_DEBUG
    node0_printf("alpha %.4g rsqnew %.4g rsq %.4g\n", alpha[0], rsqnew, rsq);
#endif

    // alpha_ip11 -- note shifted indices with respect to Eq. 2.43!
    for (j = 1; j < Norder; j++) {
      if (converged[j] == 0)
        alpha[j] = alpha[0] * zeta_ip1[j] * beta_i[j] / (zeta_i[j] * beta_i[0]);
    }

    // pm[j] = zeta_ip1[j] * r + alpha[j] * pm[j]
    for (k = 0; k < NFERMION; k++) {
      FORALLSITES(i, s) {
        scalar_mult_matrix(&(rm[k][i]), (Real)zeta_ip1[0], &(mpm[k][i]));
        scalar_mult_matrix(&(pm0[k][i]), (Real)alpha[0], &(pm0[k][i]));
        sum_matrix(&(mpm[k][i]), &(pm0[k][i]));
        for (j = 1; j < Norder; j++) {
          if (converged[j] == 0) {
            scalar_mult_matrix(&(rm[k][i]), (Real)zeta_ip1[j], &(mpm[k][i]));
            scalar_mult_matrix(&(pm[j][k][i]), (Real)alpha[j], &(pm[j][k][i]));
            sum_matrix(&(mpm[k][i]), &(pm[j][k][i]));
          }
        }
      }
    }

    // Test for convergence
    rsq = rsqnew;
    for (j = 1; j < Norder; j++) {
      if (converged[j] == 0) {
        rsqj = rsq * zeta_ip1[j] * zeta_ip1[j];
        if (rsqj <= rsqstop) {
          converged[j] = 1;
#ifdef CG_DEBUG
          node0_printf(" psi%d converged in %d steps, rsq = %.4g\n",
                       j, N_iter, rsqj);
#endif
        }
      }
    }
#ifdef CG_DEBUG
    if ((N_iter / 10) * 10 == N_iter) {
      node0_printf("iter %d residue %.4g\n", N_iter, (double)(rsq));
      fflush(stdout);
    }
#endif

    // Scroll scalars
    for (j = 0; j < Norder; j++) {
      if (converged[j] == 0) {
        beta_im1[j] = beta_i[j];
        zeta_im1[j] = zeta_i[j];
        zeta_i[j] = zeta_ip1[j];
      }
    }
  }
  if (rsq > rsqstop)
    node0_printf(" multi CONGRAD not converged\n rsq = %.4g\n", rsq);

  *size_r = rsq;

  // Test inversion
#ifdef CG_DEBUG
  for (j = 0; j < Norder; j++) {
    source_norm = 0;
    for (k = 0; k < NFERMION; k++) {
      FORALLSITES(i, s)
        source_norm += (double)realtrace(&(psim[j][k][i]), &(psim[j][k][i]));
    }
    g_doublesum(&source_norm);
    node0_printf("Norm of psim %d shift %.4g is %.4g\n",
                 j, shift[j], source_norm);

    DSq(psim[j], mpm);              // mpm = (D^2).psim[j]
    source_norm = 0;                // Re-using for convenience
    for (k = 0; k < NFERMION; k++) {
      FORALLSITES(i, s) {           // Add shift.psi and subtract src
        scalar_mult_sum_matrix(&(psim[j][k][i]), shift[j],  &(mpm[k][i]));
        dif_matrix(&(src[k][i]), &(mpm[k][i]));
        source_norm += (double)realtrace(&(mpm[k][i]), &(mpm[k][i]));
      }
    }
    g_doublesum(&source_norm);
    node0_printf("%d Test of (D^2 + shift)psi - src: %.4g\n", j, source_norm);
  }
#endif

  for (i = 1; i < Norder; i++) {
    for (k = 0; k < NFERMION; k++)
      free(pm[i][k]);
    free(pm[i]);
  }

  free(pm);
  free(zeta_i);
  free(zeta_ip1);
  free(zeta_im1);
  free(beta_im1);
  free(beta_i);
  free(alpha);
  free(converged);
  return iteration;
}
#endif
// -----------------------------------------------------------------
