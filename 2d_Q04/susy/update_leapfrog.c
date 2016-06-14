// -----------------------------------------------------------------
// Update lattice
// Leapfrog integrator
// Begin at "integral" time, with H and U evaluated at the same time

// Uncomment to print out debugging messages
//#define UPDATE_DEBUG
#include "susy_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>         // For "finite"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update_uu(Real eps) {
  register int i, mu;
  register site *s;

  FORALLDIR(mu) {
    FORALLSITES(i, s)
      scalar_mult_sum_matrix(&(s->mom[mu]), eps, &(s->link[mu]));
  }

  // Update plaquette determinants, DmuUmu and Fmunu with new links
  // (Needs to be done before calling gauge_force)
  compute_plaqdet();
  compute_Uinv();
  compute_DmuUmu();
  compute_Fmunu();
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update_step(double *fnorm, double *gnorm,
                Twist_Fermion **src, Twist_Fermion ***psim) {

  int step, iters = 0, n;
  Real final_rsq, eps = traj_length / (Real)nsteps[0], tr;
  node0_printf("eps %.4g\n", eps);

  // First u(t/2)
  update_uu(0.5 * eps);

  for (step = 0; step < nsteps[0]; step++) {
    // Inner steps p(t) u(t)
    tr = gauge_force(eps);
    *gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;

#ifndef PUREGAUGE
    for (n = 0; n < Nroot; n++) {
      // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
      iters += congrad_multi_field(src[n], psim[n], niter, rsqmin, &final_rsq);
      tr = fermion_force(eps, src[n], psim[n]);
      fnorm[n] += tr;
      if (tr > max_ff[n])
        max_ff[n] = tr;
    }
#endif

    if (step < nsteps[0] - 1)
      update_uu(eps);
    else                // Final u(t/2)
      update_uu(0.5 * eps);
  }
  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update() {
  int j, n, iters = 0;
  Real final_rsq;
  double startaction, endaction, change;
  Twist_Fermion **src = malloc(Nroot * sizeof(**src));
  Twist_Fermion ***psim = malloc(Nroot * sizeof(***psim));

  for (n = 0; n < Nroot; n++) {
    src[n] = malloc(sites_on_node * sizeof(Twist_Fermion));
    psim[n] = malloc(Norder * sizeof(Twist_Fermion*));
    for (j = 0; j < Norder; j++)
      psim[n][j] = malloc(sites_on_node * sizeof(Twist_Fermion));
  }

  // Refresh the momenta
  ranmom();

  // Set up the fermion variables, if needed
#ifndef PUREGAUGE
  // Compute g and src = (Mdag M)^(1 / 8) g
  for (n = 0; n < Nroot; n++)
    iters += grsource(src[n]);

  // Do a CG to get psim,
  // rational approximation to (Mdag M)^(-1 / 4) src = (Mdag M)^(-1 / 8) g
  for (j = 0; j < Norder; j++)
    shift[j] = shift4[j];
#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_leapfrog -- original action\n");
#endif
  // congrad_multi_field initializes psim
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi_field(src[n], psim[n], niter, rsqmin, &final_rsq);
#endif // ifndef PUREGAUGE

  // Find initial action
  startaction = action(src, psim);
  gnorm = 0.0;
  max_gf = 0.0;
  for (n = 0; n < Nroot; n++) {
    fnorm[n] = 0.0;
    max_ff[n] = 0.0;
  }

#ifdef HMC_ALGORITHM
  Real xrandom;   // For accept/reject test
  // Copy link field to old_link
  gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));
#endif
  // Do microcanonical updating
  iters += update_step(fnorm, &gnorm, src, psim);

  // Find ending action
  // Since update_step ended on a gauge update,
  // need to do conjugate gradient to get (Mdag M)^(-1 / 4) chi
#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_leapfrog -- new action\n");
#endif
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi_field(src[n], psim[n], niter, rsqmin, &final_rsq);
  endaction = action(src, psim);
  change = endaction - startaction;
#ifdef HMC_ALGORITHM
  // Reject configurations giving overflow
#ifndef HAVE_IEEEFP_H
  if (fabs((double)change) > 1e20) {
#else
  if (!finite((double)change)) {
#endif
    node0_printf("WARNING: Correcting Apparent Overflow: Delta S = %.4g\n",
                 change);
    change = 1.0e20;
  }

  // Decide whether to accept, if not, copy old link field back
  // Careful -- must generate only one random number for whole lattice
  if (this_node == 0)
    xrandom = myrand(&node_prn);
  broadcast_float(&xrandom);
  if (exp(-change) < (double)xrandom) {
    if (traj_length > 0.0) {
      gauge_field_copy(F_OFFSET(old_link[0]), F_OFFSET(link[0]));
      compute_plaqdet();
      compute_Uinv();
      compute_DmuUmu();
      compute_Fmunu();
    }
    node0_printf("REJECT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
  }
  else {
    node0_printf("ACCEPT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
  }
#else
  // Only print check if not doing HMC
  node0_printf("CHECK: delta S = %.4g\n", (double)(change));
#endif // ifdef HMC

  for (n = 0; n < Nroot; n++) {
    free(src[n]);
    for (j = 0; j < Norder; j++)
      free(psim[n][j]);
    free(psim[n]);
  }
  free(src);
  free(psim);

  if (traj_length > 0) {
    node0_printf("IT_PER_TRAJ %d\n", iters);
    node0_printf("MONITOR_FORCE_GAUGE    %.4g %.4g\n",
                 gnorm / (double)(2 * nsteps[0]), max_gf);
    for (n = 0; n < Nroot; n++) {
      node0_printf("MONITOR_FORCE_FERMION%d %.4g %.4g\n",
                   n, fnorm[n] / (double)(2 * nsteps[0]), max_ff[n]);
    }
    return iters;
  }
  else
    return -99;
}
// -----------------------------------------------------------------
