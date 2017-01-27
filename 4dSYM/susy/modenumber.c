// -----------------------------------------------------------------
// Calculation of the mode number with the Giusti--Luescher method
// Adapted from adjoint SU(N) code by Georg Bergner
#include "susy_includes.h"

// TODO: Is there any way we can use the multi-shift solver
//       to do the inversions for many Omega simultaneously?
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dest = src - 2Om_*^2 (DDdag + Om_*^2)^(-1) src
void X(const Real OmStar, Twist_Fermion *src, Twist_Fermion *dest) {
  register int i;
  register site *s;
  Real m2OmSq = -2.0 * OmStar * OmStar, size_r;
  Twist_Fermion **psim;

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(Twist_Fermion));
  shift[0] = OmStar * OmStar;

  congrad_multi(src, psim, niter, rsqmin, &size_r);
  FORALLSITES(i, s)
    scalar_mult_add_TF(&(src[i]), &(psim[0][i]), m2OmSq, &(dest[i]));
}

// dest = (2X^2 - 1 - step_eps) src / (1 - step_eps)
// With = X^2 = 1 - 4 OmStar^2 (DDdag + OmStar^2)^(-1)
//                + 4 OmStar^4 (DDdag + OmStar^2)^(-2)
// Uses z_rand for temporary storage (tempTF used by CG!)
void Z(const Real OmStar, Twist_Fermion *src, Twist_Fermion *dest) {
  register int i;
  register site *s;
  Real OmSq = OmStar * OmStar, scale = 2.0 / (1.0 - step_eps);
  Real m4OmSq = -4.0 * scale * OmSq, size_r;
  Real OmFour = scale * 4.0 * OmSq * OmSq;
  Twist_Fermion **psim;

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(Twist_Fermion));
  shift[0] = OmSq;

  // This is more compact, but the subtraction in X(src)
  // seems to leave it relatively poorly conditioned,
  // increasing CG iterations by almost 10% even for a small 4nt4 test
//  X(OmStar, src, z_rand);
//  X(OmStar, z_rand, dest);
//  Real toAdd = (-1.0 - step_eps) / (1.0 - step_eps);
//  FORALLSITES(i, s) {
//    scalar_mult_TF(&(dest[i]), scale, &(dest[i]));
//    scalar_mult_sum_TF(&(src[i]), toAdd, &(dest[i]));
//  }

  congrad_multi(src, psim, niter, rsqmin, &size_r);
  FORALLSITES(i, s)
    copy_TF(&(psim[0][i]), &(z_rand[i]));
  congrad_multi(z_rand, psim, niter, rsqmin, &size_r);
  FORALLSITES(i, s) {
    scalar_mult_TF(&(psim[0][i]), OmFour, &(dest[i]));
    scalar_mult_sum_TF(&(z_rand[i]), m4OmSq, &(dest[i]));
    sum_TF(&(src[i]), &(dest[i]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Clenshaw algorithm:
// P(X) src = sum_i^n c[i] T[i] src = (b[0] - X b[1]) src,
// where b[i] = c[i] + 2zb[i + 1] - b[i + 2], b[n] = b[n + 1] = 0
// Use TODO... for temporary storage (z_rand and tempTF used by Z!)
// Return number of CG calls
int clenshaw(const Real OmStar, Twist_Fermion *src,
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
  Z(OmStar, src, dest);
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
  Z(OmStar, bp2, bp1);
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
    Z(OmStar, bp1, bn);
    CGcalls += 2;
    FORALLSITES(i, s) {
      scalar_mult_sum_TF(&(src[i]), step_coeff[0], &(bn[i]));
      sub_TF(&(bn[i]), &(bp2[i]), &(dest[i]));
    }
    return CGcalls;
  }

  for (k = step_order - 4; k--;) {   // TODO: Relies on unsigned int?
    Z(OmStar, bp1, bn);
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
  Z(OmStar, bp1, bn);
  CGcalls += 2;

  FORALLSITES(i, s) {
    scalar_mult_sum_TF(&(src[i]), step_coeff[0], &(bn[i]));
    sub_TF(&(bn[i]), &(bp2[i]), &(dest[i]));
  }
  return CGcalls;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wrapper for step function approximated by h(X) = [1 - X p(X)^2] / 2
// Return number of CG calls
int step(const Real OmStar, Twist_Fermion *src, Twist_Fermion *dest, Twist_Fermion *XPXSq,
   Twist_Fermion **workv) {
  register int i;
  register site *s;
  int CGcalls = 0;

  // dest = P(X^2) src temporarily
  CGcalls = clenshaw(OmStar, src, dest, workv);

  // dest = (src - X P(X^2) src) / 2
  X(OmStar, dest, XPXSq);
  CGcalls++;
  FOREVENSITES(i, s) {
    sub_TF(&(src[i]), &(XPXSq[i]), &(dest[i]));
    scalar_mult_TF(&(dest[i]), 0.5, &(dest[i]));
  }
  return CGcalls;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use TODO... for temporary storage (z_rand and tempTF used by clenshaw!)
void compute_mode() {
  register int i;
  register site *s;
  int k, l, CGcalls, sav_iters;
  Real OmStar, dtime, norm = 1.0 / (Real)Nstoch, tr;
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
      FORALLSITES(i, s)
        copy_TF(&(source[l][i]), &(v0[i]));   // TODO: Can replace v0

//      CGcalls = step(OmStar, source[l], v1, z_rand, workv);   // !!! z_rand tmp stor
//      CGcalls += step(OmStar, v1, v0, z_rand, workv);   // !!! z_rand tmp stor

      // h(X) = 0.5 * (1 - X P(Z(X^2)))
      // Z(X) = 2X / (max - min) - (max + min) / (max - min);
      // X = 1 - 2Om_*^2 (DDdag + Om_*^2)^(-1)
      // First step function application
      CGcalls = clenshaw(OmStar, v0, v1, workv);
      X(OmStar, v1, v2);
      CGcalls++;
      FORALLSITES(i, s) {
        sub_TF(&(v0[i]), &(v2[i]), &(v1[i]));
        scalar_mult_TF(&(v1[i]), 0.5, &(v1[i]));
      }

      // Second step function application
      CGcalls += clenshaw(OmStar, v1, v0, workv);
      X(OmStar, v0, v2);
      CGcalls++;
      FORALLSITES(i, s) {
        sub_TF(&(v1[i]), &(v2[i]), &(v0[i]));
        scalar_mult_TF(&(v0[i]), 0.5, &(v0[i]));
      }

      // Mode number is now just magnitude of v0
      tr = 0.0;
      FORALLSITES(i, s)
        tr += magsq_TF(&(v0[i]));
      g_doublesum(&tr);
      mode[k] += tr;
      err[k] += tr * tr;

      // Monitor iterations and timing
      dtime += dclock();
      node0_printf("Stoch est %d of %d for Omega_* = %.4g : ",
                   l, Nstoch, OmStar);
      node0_printf("%d iter from %d CG calls in %.4g seconds\n",
                   total_iters - sav_iters, CGcalls, dtime);
    }

    // Average over volume and stochastic estimators,
    // and estimate standard deviations
    mode[k] *= norm;
    tr = err[k] * norm;
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
