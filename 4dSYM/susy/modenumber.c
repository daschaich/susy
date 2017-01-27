// -----------------------------------------------------------------
// Calculation of the mode number with the Giusti--Luescher method
// Adapted from adjoint SU(N) code by Georg Bergner
#include "susy_includes.h"

//#define DEBUG_TEST_CLENSH
//#define DEBUG_TEST_PROJ
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// X = 1 - 2 OmStar^2 (DDdag + OmStar^2)^(-1)
// X^2 = 1 - 4 OmStar^2 Inv + 4 OmStar^4 Inv Inv
// z(X)= 2X / (max - min) - (max + min) / (max - min)
// This is z(X^2)
void X(const Real OmStar,
    Twist_Fermion *src, Twist_Fermion *dest,
    Twist_Fermion** workv) {

  register int i;
  register site *s;
  Real OmSq = OmStar * OmStar, rescale = 2.0 / (1.0 - step_eps);
  Real OmFour = rescale * 4.0 * OmSq * OmSq;
  Real m4OmSq = -4.0 * rescale * OmSq, size_r;
  Twist_Fermion *tmpv = workv[3], **psim;

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(Twist_Fermion));
  shift[0] = OmSq;

  congrad_multi(src, psim, niter, rsqmin, &size_r);
  FORALLSITES(i, s)
    copy_TF(&(psim[0][i]), &(tmpv[i]));
  congrad_multi(tmpv, psim, niter, rsqmin, &size_r);
  FORALLSITES(i, s)
    copy_TF(&(psim[0][i]), &(dest[i]));
  FORALLSITES(i, s) {
    scalar_mult_TF(&(dest[i]), OmFour, &(dest[i]));
    scalar_mult_sum_TF(&(tmpv[i]), m4OmSq, &(dest[i]));
    sum_TF(&(src[i]), &(dest[i]));
  }
  /*Out = (rescale * 4.0 * OmStar^4) * Out - (rescale * 4.0 * OmStar^2) * tmp_
   + In*/
}

// dest = src - (2Om_*^2) (DDdag + Om_*^2)^(-1) src
void simpleX(const Real OmStar, Twist_Fermion *src, Twist_Fermion *dest) {
  register int i;
  register site *s;
  Real OmSq = OmStar * OmStar, size_r;
  Twist_Fermion **psim;

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(Twist_Fermion));
  shift[0] = OmSq;

  congrad_multi(src, psim, niter, rsqmin, &size_r);
  FORALLSITES(i, s)
    scalar_mult_add_TF(&(src[i]), &(psim[0][i]), -(2.0 * OmSq), &(dest[i]));
  /*      Out = In - (2.0 * Om_*^2) * tmp_;*/
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// P(z(x^2))
// This uses the Clenshaw-Algorithm: p_n(D) v = sum_n a_n D^n v =(a_0+D(a_1+D(a_2+\ldots(a_{n-2}+D(a_{n-1} +a_{n}D)))))v
// or: v_0 = a_n v;\; v_1 = Dv_0+a_{n-1}v;\; v_2 = Dv_1+a_{n-2}v;\;\ldots\; v_n = p_n(D)v
// Return number of CG calls
int clensh(const Real OmStar, Twist_Fermion *src,
           Twist_Fermion *dest, Twist_Fermion **workv) {

  register unsigned int i;
  register site* s;
  unsigned int k, CGcalls = 0;
  Real tr;
  Twist_Fermion *bp2, *bp1, *bn, *tmp;

  if (step_order == 0) {
    FORALLSITES(i, s)
      clear_TF(&(dest[i]));
    return CGcalls;
  }

  if (step_order == 1) {
    FORALLSITES(i, s)
      scalar_mult_TF(&(src[i]), step_coeff[0], &(dest[i]));
    return CGcalls;
  }
  X(OmStar, src, dest, workv);
  CGcalls += 2;

  if (step_order == 2) {
    FORALLSITES(i, s) {
      scalar_mult_TF(&(src[i]), step_coeff[0], &(dest[i]));
      scalar_mult_sum_TF(&(dest[i]), step_coeff[1], &(dest[i]));
    }
    return CGcalls;
  }
  bp2 = dest;
  bp1 = workv[4];
  bn = workv[5];
  tr = 2.0 * step_coeff[step_order - 1];    // Remove from site loop
  FORALLSITES(i, s) {
    scalar_mult_TF(&(bp2[i]), tr, &(bp2[i]));
    scalar_mult_sum_TF(&(src[i]), step_coeff[step_order - 2], &(bp2[i]));
  }
  X(OmStar, bp2, bp1, workv);
  CGcalls += 2;

  if (step_order == 3) {
    tr = step_coeff[step_order - 3] - step_coeff[step_order - 1];
    FORALLSITES(i, s) {
      scalar_mult_TF(&(src[i]), tr, &(dest[i]));
      sum_TF(&(bp1[i]), &(dest[i]));
    }
    return CGcalls;
  }

  tr = step_coeff[step_order - 3] - step_coeff[step_order - 1];
  FORALLSITES(i, s) {
    scalar_mult_TF(&(bp1[i]), 2.0, &(bp1[i]));
    scalar_mult_sum_TF(&(src[i]), tr, &(bp1[i]));
  }

  if (step_order == 4) {
    X(OmStar, bp1, bn, workv);
    CGcalls += 2;
    FORALLSITES(i, s) {
      scalar_mult_sum_TF(&(src[i]), step_coeff[0], &(bn[i]));
      sub_TF(&(bn[i]), &(bp2[i]), &(dest[i]));
    }
    return CGcalls;
  }

  for (k = step_order - 4; k--;) {   // TODO: Relies on unsigned int?
    X(OmStar, bp1, bn, workv);
    CGcalls += 2;
    FORALLSITES(i, s) {
      scalar_mult_TF(&(bn[i]), 2.0, &(bn[i]));
      scalar_mult_sum_TF(&(src[i]), step_coeff[k + 1], &(bn[i]));
      dif_TF(&(bp2[i]), &(bn[i]));
    }
    tmp = bp2;
    bp2 = bp1;
    bp1 = bn;
    bn = tmp;
  }
  X(OmStar, bp1, bn, workv);
  CGcalls += 2;

  FORALLSITES(i, s) {
    scalar_mult_sum_TF(&(src[i]), step_coeff[0], &(bn[i]));
    sub_TF(&(bn[i]), &(bp2[i]), &(dest[i]));
  }

  return CGcalls;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void compute_mode() {
  register int i;
  register site *s;
  int k, l, CGcalls, sav_iters;
  Real OmStar, tr, norm = 1.0 / (Real)(16 * DIMF * sites_on_node);
  double dtime;
  Twist_Fermion *v0 = malloc(sites_on_node * sizeof(Twist_Fermion));
  Twist_Fermion *v1 = malloc(sites_on_node * sizeof(Twist_Fermion));
  Twist_Fermion *v2 = malloc(sites_on_node * sizeof(Twist_Fermion));
  Twist_Fermion *workv[6];
  workv[0] = malloc(sites_on_node * sizeof(Twist_Fermion));
  workv[1] = malloc(sites_on_node * sizeof(Twist_Fermion));
  workv[2] = malloc(sites_on_node * sizeof(Twist_Fermion));
  workv[3] = malloc(sites_on_node * sizeof(Twist_Fermion));
  workv[4] = malloc(sites_on_node * sizeof(Twist_Fermion));
  workv[5] = malloc(sites_on_node * sizeof(Twist_Fermion));
#ifdef DEBUG_TEST_CLENSH
  Z2source();
  clensh(OmStar, v0, v1, workv);
  tr = 0.0;
  FORALLSITES(i, s)
    tr += (double)magsq_TF(&(v1[i]));
  tr *= 1.0 / (double)(16 * DIMF * sites_on_node);
  g_doublesum(&tr);
  node0_printf("TESTMODE Clensh: P^2(z)= %lg\n", tr);

  OmStar = 0.2;
  simpleX(OmStar, v0, v1);
  tr = 0.0;
  FORALLSITES(i, s)
    tr += (double)magsq_TF(&(v1[i]));
  tr *= 1.0 / (double)(16 * DIMF * sites_on_node);
  g_doublesum(&tr);
  node0_printf("TESTMODE Clensh: xOperator^2(0.2, 0.4)= %lg\n", tr);
  X(OmStar, v0, v1, workv);
  tr = 0.0;
  FORALLSITES(i, s)
    tr += (double)magsq_TF(&(v1[i]));

  tr *= 1.0 / (double)(16 * DIMF * sites_on_node);
  g_doublesum(&tr);
  node0_printf("TESTMODE Clensh: F(z(xOperator(0.2, 0.4)))^2 = %lg\n", tr);
  free(v0);
  free(v1);
  free(v2);
  free(workv[0]);
  free(workv[1]);
  free(workv[2]);
  free(workv[3]);
  free(workv[4]);
  free(workv[5]);
  return 0;
#endif

  // Set up stochastic Z2 random sources
  for (l = 0; l < Nstoch; l++) {
    Z2source();
    FORALLSITES(i, s)
      copy_TF(&(z_rand[i]), &(source[l][i]));
  }

  for (k = 0; k < numOmega; k++) {
    // Scale by star
    OmStar = Omega[k] / star;

    // Initialize results and err, then average over Nstoch
    mode[k] = 0.0;
    err[k] = 0.0;
    for (l = 0; l < Nstoch; l++) {
      // Setup timing and iteration counts
      dtime = -dclock();
      sav_iters = total_iters;    // total_iters incremented by CG
      CGcalls = 0;
      FORALLSITES(i, s)
        copy_TF(&(source[l][i]), &(v0[i]));   // TODO: Can replace v0

      // h(X) = 0.5 * (1 - X P(Z(X^2)))
      // Z(X) = 2X / (max - min) - (max + min) / (max - min);
      // X = 1 - 2Om_*^2 (DDdag + Om_*^2)^(-1)
      // First step function application
      CGcalls += clensh(OmStar, v0, v1, workv);
      simpleX(OmStar, v1, v2);
      CGcalls++;
      FORALLSITES(i, s) {
        sub_TF(&(v0[i]), &(v2[i]), &(v1[i]));
        scalar_mult_TF(&(v1[i]), 0.5, &(v1[i]));
      }

      // Second step function application -- mode number is now just norm
      clensh(OmStar, v1, v0, workv);
      simpleX(OmStar, v0, v2);
      CGcalls++;
      tr = 0.0;
      FORALLSITES(i, s) {
        sub_TF(&(v1[i]), &(v2[i]), &(v0[i]));
        scalar_mult_TF(&(v0[i]), 0.5, &(v0[i]));
        tr += magsq_TF(&(v0[i]));
      }
      g_doublesum(&tr);
      mode[k] += tr;
      err[k] += tr * tr;

      // Monitor iterations and timing
      dtime += dclock();
      node0_printf("Stochastic estimator %d of %d for Omega_* = %.4g : ",
                   l, Nstoch, OmStar);
      node0_printf("%d iter from %d CG calls in %.4g seconds\n",
                   total_iters - sav_iters, CGcalls, dtime);
    }

    // Average over volume and stochastic estimators,
    // and estimate standard deviations
    mode[k] *= norm / (Real)Nstoch;
    tr = err[k] * norm * norm / (Real)Nstoch;
    err[k] = sqrtN_ov_Nm1 * sqrt((tr - mode[k] * mode[k]));
  }
  free(v0);
  free(v1);
  free(v2);
  free(workv[0]);
  free(workv[1]);
  free(workv[2]);
  free(workv[3]);
  free(workv[4]);
  free(workv[5]);
}
// -----------------------------------------------------------------
