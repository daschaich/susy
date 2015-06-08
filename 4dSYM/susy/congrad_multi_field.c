// -----------------------------------------------------------------
// Multi-mass conjugate gradient algorithm a la B. Jegerlehner
// The number of masses is runtime input Norder
// shift[Norder] is the array of mass values, set up in setup_rhmc.c
// BEWARE: The temporary vector pm[0] is never malloced or used

// At least for now we hard-code a zero initial guess
// We check all vectors for convergence and quit doing the converged ones

// In this version of the code, all the scalars are real because M = Ddag D
//#define CG_DEBUG
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return number of iterations
// src is where the source is created
// psim[Norder] are working vectors for the conjugate gradient
// MaxCG is the maximum number of iterations per restart
// RsdCG is the target residual, normalized as sqrt(r * r) / sqrt(src * src)
// size_r is the final obtained residual, rsq < rsqmin * source_norm
int congrad_multi_field(Twist_Fermion *src, Twist_Fermion **psim,
                        int MaxCG, Real RsdCG, Real *size_r) {

  register int i, j;
  register site *s;
  int N_iter, iteration = 0;
  int *converged = malloc(Norder * sizeof(*converged));
  Real floatvar, floatvar2;     // SSE kluge
  Real *floatvarj = malloc(Norder * sizeof(*floatvarj));
  Real *floatvark = malloc(Norder * sizeof(*floatvark));
  double rsq = 0, rsqnew, source_norm = 0, errormin, rsqstop, c1, c2, cd;
  double *zeta_i   = malloc(Norder * sizeof(*zeta_i));
  double *zeta_im1 = malloc(Norder * sizeof(*zeta_im1));
  double *zeta_ip1 = malloc(Norder * sizeof(*zeta_ip1));
  double *beta_i   = malloc(Norder * sizeof(*beta_i));
  double *beta_im1 = malloc(Norder * sizeof(*beta_im1));
  double *alpha    = malloc(Norder * sizeof(*alpha));
  double rsqj;
  complex ctmp;
  Twist_Fermion *mpm = malloc(sites_on_node * sizeof(*mpm));
  Twist_Fermion *pm0 = malloc(sites_on_node * sizeof(*pm0));
  Twist_Fermion *rm  = malloc(sites_on_node * sizeof(*rm));
  Twist_Fermion **pm  = malloc(Norder * sizeof(**pm));
  for (i = 1; i < Norder; i++)    // !!!
    pm[i] = malloc(sites_on_node * sizeof(Twist_Fermion));

  // Make sure irrep links used by the fermion operator are up to date
  fermion_rep();

  // Initialize zero initial guess, etc.
  // dest = 0, r = source, pm[j] = r
  errormin = RsdCG * RsdCG;
  for (i = 0; i < Norder; i++)
    converged[i] = 0;
  FORALLSITES(i, s) {
    copy_TF(&(src[i]), &(rm[i]));
    copy_TF(&(rm[i]), &(pm0[i]));
    clear_TF(&(psim[0][i]));
    for (j = 1; j < Norder; j++) {
      clear_TF(&(psim[j][i]));
      copy_TF(&(rm[i]), &(pm[j][i]));
    }
  }

  FORALLSITES(i, s)
    source_norm += (double)magsq_TF(&(src[i]));

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
    hdelta0_field(pm0, mpm);
    iteration++;
    total_iters++;
    FORALLSITES(i,s)
      scalar_mult_add_TF(&(mpm[i]), &(pm0[i]), shift[0], &(mpm[i]));

    // beta_i[0] = -(r, r) / (pm, Mpm)
    cd = 0;
    FORALLSITES(i, s) {
      ctmp = TF_dot(&(pm0[i]), &(mpm[i]));
      cd += ctmp.real;
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
    floatvar = -(Real)beta_i[0];
    for (j = 1; j < Norder; j++) {
      if (converged[j] == 0)
        floatvarj[j] = -(Real)beta_i[j];
    }

    FORALLSITES(i, s) {
      scalar_mult_add_TF(&(psim[0][i]), &(pm0[i]), floatvar, &(psim[0][i]));
      for (j = 1; j < Norder; j++) {
        if (converged[j] == 0)
          scalar_mult_add_TF(&(psim[j][i]), &(pm[j][i]), floatvarj[j],
                             &(psim[j][i]));
      }
    }

    // r = r + beta[0] * mp
    floatvar = (Real)beta_i[0];
    FORALLSITES(i, s)
      scalar_mult_add_TF(&(rm[i]), &(mpm[i]),floatvar, &(rm[i]));

    // alpha_ip1[j]
    rsqnew = 0;
    FORALLSITES(i, s)
      rsqnew += (double)magsq_TF(&(rm[i]));

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
    floatvar  = (Real)zeta_ip1[0];
    floatvar2 = (Real)alpha[0];
    for (j = 1; j < Norder; j++) {
      floatvarj[j] = (Real)zeta_ip1[j];
      floatvark[j] = (Real)alpha[j];
    }
    FORALLSITES(i, s) {
      scalar_mult_TF(&(rm[i]),floatvar, &(mpm[i]));
      scalar_mult_add_TF(&(mpm[i]), &(pm0[i]), floatvar2, &(pm0[i]));
      for (j = 1; j < Norder; j++) {
        if (converged[j] == 0) {
          scalar_mult_TF(&(rm[i]), floatvarj[j], &(mpm[i]));
          scalar_mult_add_TF(&(mpm[i]), &(pm[j][i]),floatvark[j], &(pm[j][i]));
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
          node0_printf(" vector %d converged in %d steps, rsq = %.4g\n",
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
    FORALLSITES(i, s)
      source_norm += (double)magsq_TF(&(psim[j][i]));

    g_doublesum(&source_norm);
    node0_printf("Norm of psim %d shift %.4g is %.4g\n", j, shift[j], source_norm);

    hdelta0_field(psim[j], mpm);    // mpm = (D^2 + fmass^2).psim[j]
    source_norm = 0;                // Re-using for convenience
    FORALLSITES(i, s) {             // Add shift.psi and subtract src
      scalar_mult_add_TF(&(mpm[i]), &(psim[j][i]), shift[j],  &(mpm[i]));
      scalar_mult_add_TF(&(mpm[i]), &(src[i]), -1.0, &(mpm[i]));
      source_norm += (double)magsq_TF(&(mpm[i]));
    }
    g_doublesum(&source_norm);
    node0_printf("%d Test of (D^2 + shift)psi - src: %.4g\n", j, source_norm);
  }
#endif

  for (i = 1; i < Norder; i++)
    free(pm[i]);
  free(pm);
  free(rm);
  free(mpm);
  free(pm0);
  free(zeta_i);
  free(zeta_ip1);
  free(zeta_im1);
  free(beta_im1);
  free(beta_i);
  free(alpha);
  free(converged);
  free(floatvarj);
  free(floatvark);
  return iteration;
}
// -----------------------------------------------------------------
