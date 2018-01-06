// -----------------------------------------------------------------
// Calculation of the mode number with the Giusti--Luescher method
// Adapted from adjoint SU(N) code by Georg Bergner
#include "susy_includes.h"

// TODO: Is there any way we can use the multi-shift solver
//       to do the inversions for many Omega simultaneously?
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dest = src - 2Om_*^2 (DDdag + Om_*^2)^(-1) src
void X(Twist_Fermion *src, Twist_Fermion *dest) {
  register int i;
  register site *s;
  Real m2OmSq = -2.0 * OmStar * OmStar, size_r;
  Twist_Fermion **psim;

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof **psim);
  psim[0] = malloc(sizeof(Twist_Fermion) * sites_on_node);
  shift[0] = OmStar * OmStar;

  congrad_multi(src, psim, niter, rsqmin, &size_r);
  FORALLSITES(i, s)
    scalar_mult_add_TF(&(src[i]), &(psim[0][i]), m2OmSq, &(dest[i]));
}

// dest = (2X^2 - 1 - step_eps) src / (1 - step_eps)
// With = X^2 = 1 - 4 OmStar^2 (DDdag + OmStar^2)^(-1)
//                + 4 OmStar^4 (DDdag + OmStar^2)^(-2)
// Uses z_rand for temporary storage (tempTF used by CG!)
void Z(Twist_Fermion *src, Twist_Fermion *dest) {
  register int i;
  register site *s;
  Real OmSq = OmStar * OmStar, scale = 2.0 / (1.0 - step_eps);
  Real m4OmSq = -4.0 * scale * OmSq, size_r;
  Real OmFour = scale * 4.0 * OmSq * OmSq;
  Twist_Fermion **psim;

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof **psim);
  psim[0] = malloc(sizeof(Twist_Fermion) * sites_on_node);
  shift[0] = OmSq;

  // This is more compact, but the subtraction in X(src)
  // seems to leave it relatively poorly conditioned,
  // increasing CG iterations by almost 10% even for a small 4nt4 test
//  X(src, z_rand);
//  X(z_rand, dest);
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
// P(X) src = sum_i^n c[i] T[i] src = (b[0] - Z(b[1])) src,
// where b[i] = c[i] + 2Z(b[i + 1]) - b[i + 2], b[n] = b[n + 1] = 0
// Use bj and bjp1 for temporary storage (z_rand and tempTF in use!)
// We also store bjp2 in dest for intermediate steps
// Return number of CG calls
int clenshaw(Twist_Fermion *src, Twist_Fermion *dest) {
  register unsigned int i;
  register site* s;
  int j, CGcalls = 0;

  for (j = step_order; j >= 0; j--) {
    // Construct bj src = (cj + 2Z(bjp1) - bjp2) src
    // bjp1 and bjp2 = dest come from previous iterations (initially zero)
    // We can overwrite dest with Z.bjp1
    FORALLSITES(i, s) {                       // Initialize
      scalar_mult_TF(&(src[i]), step_coeff[j], &(bj[i]));
      if (j < step_order - 1)                 // Subtract bjp2 src
        dif_TF(&(dest[i]), &(bj[i]));
    }

    // Add 2Z(bjp1) src
    if (j < step_order) {
      Z(bjp1, dest);
      CGcalls += 2;
      FORALLSITES(i, s)
        scalar_mult_sum_TF(&(dest[i]), 2.0, &(bj[i]));
    }

    // Now shift dest = bjp2 <-- bjp1 and bjp1 <-- bj for next iteration
    if (j > 0) {
      FORALLSITES(i, s) {
        copy_TF(&(bjp1[i]), &(dest[i]));
        copy_TF(&(bj[i]), &(bjp1[i]));
      }
    }
  }

  // We now have bj = b[0] src and dest = bjp2 = Z(b[1]) src
  // Complete (b[0] - Z(b[1])) src
  FORALLSITES(i, s) {
    scalar_mult_TF(&(dest[i]), -1.0, &(dest[i]));
    sum_TF(&(bj[i]), &(dest[i]));
  }
  return CGcalls;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wrapper for step function approximated by h(X) = [1 - X p(X)^2] / 2
// Use XPXSq for temporary storage (z_rand and tempTF in use!)
// Return number of CG calls
int step(Twist_Fermion *src, Twist_Fermion *dest) {
  register int i;
  register site *s;
  int CGcalls = 0;

  // dest = P(X^2) src temporarily
  CGcalls = clenshaw(src, dest);

  // dest = (src - X P(X^2) src) / 2
  X(dest, XPXSq);
  CGcalls++;
  FORALLSITES(i, s) {
    sub_TF(&(src[i]), &(XPXSq[i]), &(dest[i]));
    scalar_mult_TF(&(dest[i]), 0.5, &(dest[i]));
  }
  return CGcalls;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use hX for temporary storage (XPXSq, z_rand and tempTF in use!)
void compute_mode() {
  register int i;
  register site *s;
  int k, l, CGcalls, sav_iters;
  Real dtime, norm = 1.0 / (Real)Nstoch, tr;

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

      // Hit gaussian random vector twice with step function
      CGcalls = step(source[l], hX);
      CGcalls += step(hX, dest);

      // Mode number is now just magnitude of dest
      tr = 0.0;
      FORALLSITES(i, s)
        tr += magsq_TF(&(dest[i]));
      g_doublesum(&tr);
      mode[k] += tr;
      err[k] += tr * tr;

      // Monitor iterations and timing
      dtime += dclock();
      node0_printf("Stoch est %d of %d for Omega_* = %.4g : ",
                   l, Nstoch, OmStar);
      node0_printf("%.4g from %d iter %d inverts %.4g seconds\n",
                   tr, total_iters - sav_iters, CGcalls, dtime);
    }

    // Average over volume and stochastic estimators,
    // and estimate standard deviations
    mode[k] *= norm;
    tr = err[k] * norm;
    err[k] = sqrt1_ov_Nm1 * sqrt(fabs(tr - mode[k] * mode[k]));
  }
}
// -----------------------------------------------------------------
