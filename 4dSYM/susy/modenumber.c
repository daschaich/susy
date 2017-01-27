// -----------------------------------------------------------------
// Calculation of the mode number with the Giusti--Luescher method
// Adapted from adjoint SU(N) code by Georg Bergner
#include "susy_includes.h"

//#define DEBUG_CG
//#define DEBUG_TEST_CLENSH
//#define DEBUG_TEST_PROJ
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Simple CG version including shifts
// TODO: Should be able to replace with fmass in single-pole congrad_multi...
void invertShiftedCG(const Real shift, const Real residgoal,
    const unsigned long maxit, Twist_Fermion *src, Twist_Fermion *dest,
    Real *resid, int *iterations, Twist_Fermion **workv) {

  int i, iter;
  register site* s;
  Real aval, beta;
  double source_norm = 0.0, rsq, rsqstop, oldrsq, dvecz, dvecdvec;
  complex tc;
  Twist_Fermion *tmp = workv[0], *dvec = workv[1], *rvec = workv[2];

  DSq(dest, tmp);
  FORALLSITES(i, s) {
    source_norm += (double)magsq_TF(&(src[i]));

    scalar_mult_mult_add_TF(&(tmp[i]), -1.0, &(dest[i]), -shift, &(rvec[i]));
    sum_TF(&(src[i]), &(rvec[i])); // dvec = rvec = src-Dsq*dest-shift*dest
    copy_TF(&(rvec[i]), &(dvec[i]));
    rsq += (double)magsq_TF(&(rvec[i]));
  }
  g_doublesum(&source_norm);
  g_doublesum(&rsq);
  if (source_norm < 0.01)     // Avoid problems for too-small initial vectors
    rsqstop = residgoal;
  else
  rsqstop = residgoal * source_norm;

#ifdef DEBUG_CG
  node0_printf("Initialized shifted CG ");
  node0_printf("with res %.4g, src_norm %.4g rsqstop %.4g\n",
               rsq, source_norm, rsqstop);
#endif

  if (rsq < rsqstop) { // already converged with initial vector
    *iterations = 0;
    *resid = rsq;
    return;
  }

  // Start iterations
  dvecdvec = rsq;
  for (iter = 0; iter < maxit; ++iter) {
    oldrsq = rsq;
    DSq(dvec, tmp);
    //aval =|rvec_i|^2/(<dvec_i z> + shift <dvec_i dvec_i>)
    dvecz = 0.0;
    FORALLSITES(i, s) {
      tc = TF_dot(&(dvec[i]), &(tmp[i]));
      dvecz += tc.real;
    }
    g_doublesum(&dvecz);
    aval = rsq / (dvecz + shift * dvecdvec);

    rsq = 0.0;
    FORALLSITES(i, s) {
      scalar_mult_sum_TF(&(dvec[i]), aval, &(dest[i])); // dest+= aval*dvec
      scalar_mult_sum_TF(&(tmp[i]), -aval, &(rvec[i])); // rvec-= aval*(tmp+shift*dvec)
      scalar_mult_sum_TF(&(dvec[i]), -aval * shift, &(rvec[i]));
      rsq += (double)magsq_TF(&(rvec[i]));
    }
    g_doublesum(&rsq);

#ifdef DEBUG_CG
    node0_printf("Mode CG it %d resid %.4g res_goal %.4g\n",
                 iter + 1, rsq, rsqstop);
#endif

    if (rsq < rsqstop) {
      // This is to ensure the correctness of the inversion.
      // If inversion is not correct this forces a complete restart of the iterations.
      DSq(dest, tmp);
      FORALLSITES(i, s) {
        scalar_mult_mult_add_TF(&(tmp[i]), -1.0, &(dest[i]), -shift,
            &(rvec[i]));
        sum_TF(&(src[i]), &(rvec[i]));
        copy_TF(&(rvec[i]), &(dvec[i]));

        rsq += magsq_TF(&(rvec[i]));
      }
      g_doublesum(&rsq);

#ifdef DEBUG_CG
      node0_printf("Checking Final Mode CG res =%g\n", rsq);
#endif

      if (rsq < rsqstop) {
#ifdef DEBUG_CG
        node0_printf("Mode CG it =%d, resid =%g, res_goal =%g\n", iter + 1,
            rsq, rsqstop);
#endif
        *iterations = iter + 1;
        *resid = rsq;
        return;
      }

      // complete re-initialization
      dvecdvec = rsq;
    }
    else {
      beta = rsq / oldrsq;

      dvecdvec = 0.0;
      FORALLSITES(i, s) {
        scalar_mult_add_TF(&(rvec[i]), &(dvec[i]), beta, &(dvec[i]));
        dvecdvec += (double)magsq_TF(&(dvec[i]));
      }
      g_doublesum(&dvecdvec);
    }
  }

#ifdef DEBUG_CG
  node0_printf("Mode CG Failed convergence it =%d, resid =%g, res_goal =%g\n",
      ite + 1, rsq, rsqstop);
#endif

  *iterations = iter + 1;
  *resid = rsq;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// X = 1 - 2 OmStar^2 (DDdag + OmStar^2)^(-1)
// X^2 = 1 - 4 OmStar^2 Inv + 4 OmStar^4 Inv Inv
// z(X)= 2X / (max - min) - (max + min) / (max - min)
// This is z(X^2)
void applyXoperator(const Real OmStar,
    Twist_Fermion *src, Twist_Fermion *dest,
    Twist_Fermion** workv) {

  Real OmSq = OmStar * OmStar, rescale = 2.0 / (1.0 - step_eps);
  Real OmFour = OmSq * OmSq;
  Real subtract = (1.0 + step_eps) / (1.0 - step_eps);
  Real idfact = rescale - subtract;
  Twist_Fermion* tmpv = workv[3];
  Real resid;
  int iter;
  register site* s;
  int i;
  invertShiftedCG(OmSq, rsqmin, niter, src, tmpv, &resid, &iter, workv);
  invertShiftedCG(OmSq, rsqmin, niter, tmpv, dest, &resid, &iter, workv);
  FORALLSITES(i, s) {
    scalar_mult_mult_add_TF(&(dest[i]), (rescale * 4.0 * OmFour), &(tmpv[i]),
                                       -(rescale * 4.0 * OmSq), &(dest[i]));
    scalar_mult_sum_TF(&(src[i]), idfact, &(dest[i]));
  }
  /*Out = (rescale * 4.0 * OmStar^4) * Out - (rescale * 4.0 * OmStar^2) * tmp_
   + idfact * In;*/
}

// dest = src - (2Om_*^2) (DDdag + Om_*^2)^(-1) src
void applyXoperatorSimple(const Real OmStar,
    Twist_Fermion *src,
    Twist_Fermion *dest, Twist_Fermion **workv) {

  Real OmSq = OmStar * OmStar;
  Twist_Fermion *tmpv = workv[3];
  Real resid;
  int iter;
  register site* s;
  int i;
  invertShiftedCG(OmSq, rsqmin, niter, src, tmpv, &resid, &iter, workv);
  FORALLSITES(i, s)
    scalar_mult_add_TF(&(src[i]), &(tmpv[i]), -(2.0 * OmSq), &(dest[i]));
  /*      Out = In - (2.0 * Om_*^2) * tmp_;*/
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// P(z(x^2))
// This uses the Clenshaw-Algorithm: p_n(D) v = sum_n a_n D^n v =(a_0+D(a_1+D(a_2+\ldots(a_{n-2}+D(a_{n-1} +a_{n}D)))))v
// or: v_0 = a_n v;\; v_1 = Dv_0+a_{n-1}v;\; v_2 = Dv_1+a_{n-2}v;\;\ldots\; v_n = p_n(D)v
void clensh(const Real OmStar, Twist_Fermion *src,
            Twist_Fermion *dest, Twist_Fermion **workv) {

  register unsigned int i, k;
  register site* s;
  Twist_Fermion *bp2, *bp1, *bn, *tmp;

  if (step_order == 0) {
    FORALLSITES(i, s)
      clear_TF(&(dest[i]));
    return;
  }

  if (step_order == 1) {
    FORALLSITES(i, s)
      scalar_mult_TF(&(src[i]), step_coeff[0], &(dest[i]));
    return;
  }
  applyXoperator(OmStar, src, dest, workv);

  if (step_order == 2) {
    FORALLSITES(i, s) {
      scalar_mult_mult_add_TF(&(src[i]), step_coeff[0], &(dest[i]),
                                         step_coeff[1], &(dest[i]));
    }
    return;
  }
  bp2 = dest;
  bp1 = workv[4];
  bn = workv[5];
  FORALLSITES(i, s) {
    scalar_mult_mult_add_TF(&(src[i]), step_coeff[step_order - 2], &(bp2[i]),
        2.0 * step_coeff[step_order - 1], &(bp2[i]));
  }
  applyXoperator(OmStar, bp2, bp1, workv);

  if (step_order == 3) {
    FORALLSITES(i, s) {
      scalar_mult_mult_add_TF(&(src[i]),
          step_coeff[step_order - 3] - step_coeff[step_order - 1], &(bp1[i]), 1.0,
          &(dest[i]));
    }
    return;
  }

  FORALLSITES(i, s) {
    scalar_mult_mult_add_TF(&(src[i]),
        step_coeff[step_order - 3] - step_coeff[step_order - 1],
        &(bp1[i]), 2.0, &(bp1[i]));
  }

  if (step_order == 4) {
    applyXoperator(OmStar, bp1, bn, workv);
    FORALLSITES(i, s) {
      scalar_mult_sum_TF(&(src[i]), step_coeff[0], &(bn[i]));
      scalar_mult_add_TF(&(bn[i]), &(bp2[i]), -1.0, &(dest[i]));
    }
    return;
  }

  for (k = step_order - 4; k--;) {   // TODO: Relies on unsigned int?
    applyXoperator(OmStar, bp1, bn, workv);
    FORALLSITES(i, s) {
      scalar_mult_mult_add_TF(&(src[i]), step_coeff[k + 1], &(bn[i]), 2.0,
          &(bn[i]));
      dif_TF(&(bp2[i]), &(bn[i]));
    }
    tmp = bp2;
    bp2 = bp1;
    bp1 = bn;
    bn = tmp;
  }
  applyXoperator(OmStar, bp1, bn, workv);

  FORALLSITES(i, s) {
    scalar_mult_sum_TF(&(src[i]), step_coeff[0], &(bn[i]));
    sub_TF(&(bn[i]), &(bp2[i]), &(dest[i]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void compute_mode() {
  register int i, k, l;
  register site *s;
  Real OmStar, tr, norm = 1.0 / (Real)(16 * DIMF * sites_on_node);
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
  applyXoperatorSimple(OmStar, v0, v1, workv);
  tr = 0.0;
  FORALLSITES(i, s)
    tr += (double)magsq_TF(&(v1[i]));
  tr *= 1.0 / (double)(16 * DIMF * sites_on_node);
  g_doublesum(&tr);
  node0_printf("TESTMODE Clensh: xOperator^2(0.2, 0.4)= %lg\n", tr);
  applyXoperator(OmStar, v0, v1, workv);
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
      FORALLSITES(i, s)
        copy_TF(&(source[l][i]), &(v0[i]));   // TODO: Can replace v0

      // h(X) = 0.5 * (1 - X P(Z(X^2)))
      // Z(X) = 2X / (max - min) - (max + min) / (max - min);
      // X = 1 - 2Om_*^2 (DDdag + Om_*^2)^(-1)
      // First step function application
      clensh(OmStar, v0, v1, workv);
      applyXoperatorSimple(OmStar, v1, v2, workv);
      FORALLSITES(i, s) {
        scalar_mult_mult_add_TF(&(v0[i]), 0.5, &(v2[i]), -0.5, &(v1[i]));
      }

      // Second step function application -- mode number is now just norm
      clensh(OmStar, v1, v0, workv);
      applyXoperatorSimple(OmStar, v0, v2, workv);
      tr = 0.0;
      FORALLSITES(i, s) {
        scalar_mult_mult_add_TF(&(v1[i]), 0.5, &(v2[i]), -0.5, &(v0[i]));
        tr += magsq_TF(&(v0[i]));
      }
      g_doublesum(&tr);
      mode[k] += tr;
      err[k] += tr * tr;

#ifdef DEBUG_CHECK
      node0_printf("Stochastic estimator %d of %d: %d %.4g\n",
                   l, Nstoch, k, mode[k]);
#endif
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
