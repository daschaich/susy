// -----------------------------------------------------------------
// Compute coefficients in Chebyshev approximation to spectral density
// Adapted from adjoint SU(N) code by Georg Bergner
#include "susy_includes.h"

// Uses tempTF, mpm, pm0 and rm for temporary storage
void chebyshev_coeff() {
  register int i, j, k;
  register site *s;
  Real diff = lambda_max - lambda_min, tr;
  Real two_ov_diff = 2.0 / diff, four_ov_diff = 4.0 / diff;
  Real msum_ov_diff = -(lambda_max + lambda_min) / diff;
  Real m2sum_ov_diff = 2.0 * msum_ov_diff;
  Real norm = 1.0 / (Real)(16 * DIMF * volume);
  // Since we don't call the CG, we can use its vectors
  // Set up pointers to them with more sensible names
  Twist_Fermion *vek = mpm;
  Twist_Fermion *vec_next = pm0;
  Twist_Fermion *vec_prev = rm;

  // Initialize results and cheb_err, then average over Nstoch
  for (j = 0; j < cheb_order; j++) {
    cheb_coeff[j] = 0.0;
    cheb_err[j] = 0.0;
  }
  for (k = 0; k < Nstoch; k++) {
    Z2source();

    // First special case: c[0]
    tr = 0.0;
    FORALLSITES(i, s)
      tr += magsq_TF(&(z_rand[i]));
    g_doublesum(&tr);
    cheb_coeff[0] += tr;
    cheb_err[0] += tr * tr;

    if (cheb_order < 1) {
      free(vek);
      free(vec_next);
      free(vec_prev);
      return;
    }

    // Second special case: c[1]
    DSq(z_rand, tempTF);
    tr = 0.0;
    FORALLSITES(i, s) {
      scalar_mult_TF(&(tempTF[i]), two_ov_diff, &(vek[i]));
      scalar_mult_sum_TF(&(z_rand[i]), msum_ov_diff, &(vek[i]));
      copy_TF(&(z_rand[i]), &(vec_prev[i]));

      TF_rdot_sum(&(z_rand[i]), &(vek[i]), &tr);
    }
    g_doublesum(&tr);
    cheb_coeff[1] += tr;
    cheb_err[1] += tr * tr;

    // General case: c[j]
    for (j = 2; j < cheb_order; j++) {
      DSq(vek, tempTF);
      tr = 0.0;
      FORALLSITES(i, s) {
        scalar_mult_TF(&(tempTF[i]), four_ov_diff, &(vec_next[i]));
        scalar_mult_sum_TF(&(vek[i]), m2sum_ov_diff, &(vec_next[i]));
        dif_TF(&(vec_prev[i]), &(vec_next[i]));
        TF_rdot_sum(&(z_rand[i]), &(vec_next[i]), &tr);

        // Set up for next coefficient
        copy_TF(&(vek[i]), &(vec_prev[i]));
        copy_TF(&(vec_next[i]), &(vek[i]));
      }
      g_doublesum(&tr);
      cheb_coeff[j] += tr;
      cheb_err[j] += tr * tr;
    }
#ifdef DEBUG_CHECK
    node0_printf("Stochastic estimator %d of %d:\n", k, Nstoch);
    for (j = 0; j < cheb_order; j++)
      node0_printf("%d %.4g:\n", j, cheb_coeff[j]);
#endif
  }

  // Average over (global) volume and stochastic estimators,
  // and estimate standard deviations
  for (j = 0; j < cheb_order; j++) {
    cheb_coeff[j] *= norm / (Real)Nstoch;
    tr = cheb_err[j] * norm * norm / (Real)Nstoch;
    cheb_err[j] = sqrt(fabs(tr - cheb_coeff[j] * cheb_coeff[j]));
    cheb_err[j] *= sqrt1_ov_Nm1;
  }
}
// -----------------------------------------------------------------
