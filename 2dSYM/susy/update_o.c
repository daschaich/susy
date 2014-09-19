// -----------------------------------------------------------------
// Update lattice
// Omelyan integrator with Urbach, Jansen, Schindler, Wenger multiscale
// (CPC 174:87, 2006)

// Begin at "integral" time, with H and U evaluated at the same time
// For the final accept/reject, we already have a good solution to the CG
// The last update was of the momenta

// Uncomment to print out debugging messages
//#define UPDATE_DEBUG
#include "susy_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>         // For "finite"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void ranmom() {
  register int i, j, mu;
  register site *s;
  complex grn;

  FORALLSITES(i, s) {
    for (mu = 0; mu < NUMLINK; mu++) {
      clear_su3mat_f(&(s->mom[mu]));
      for (j = 0; j < NUMGEN; j++) {
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&(s->node_prn));
        grn.imag = gaussian_rand_no(&(s->node_prn));
#endif
        c_scalar_mult_add_su3mat_f(&(s->mom[mu]), &(Lambda[j]), &grn,
                                   &(s->mom[mu]));
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update_uu(Real eps) {
  register int i, mu;
  register site *s;

  for (mu = 0; mu < NUMLINK; mu++) {
    FORALLSITES(i, s)
      scalar_mult_add_su3_matrix_f(&(s->linkf[mu]), &(s->mom[mu]), eps,
                                   &(s->linkf[mu]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Omelyan version; ``dirty'' speeded-up version
double update_gauge_step(Real eps) {
  int nsw = nsteps[1], isw;
  double norm;

#ifdef UPDATE_DEBUG
  node0_printf("gauge %d steps %.4g dt\n", nsw, eps);
#endif
  norm = gauge_force(eps*LAMBDA);
  for (isw = 1; isw <= nsw; isw++) {
    update_uu(0.5 * eps);
    norm += gauge_force(eps * LAMBDA_MID);
    update_uu(0.5 * eps);
    if (isw < nsw)
      norm += gauge_force(eps * TWO_LAMBDA);

    else
      norm += gauge_force(eps * LAMBDA);
  }
#ifndef PUREGAUGE
  fermion_rep();
#endif
  return (norm / nsw);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update_step(Real *old_cg_time,Real *cg_time,Real *next_cg_time,
                double *fnorm, double *gnorm,
                Twist_Fermion *src, Twist_Fermion **psim) {

  int iters = 0, i_multi0;
  Real final_rsq, f_eps, g_eps;

  f_eps = traj_length / (Real)nsteps[0];
  g_eps = f_eps / (Real)(2.0 * nsteps[1]);

#ifndef PUREGAUGE
  fnorm[0] += fermion_force(f_eps * LAMBDA, src, psim);
#endif

  for (i_multi0 = 1; i_multi0 <= nsteps[0]; i_multi0++) {
    *gnorm += update_gauge_step(g_eps);

#ifndef PUREGAUGE
    // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
    iters += congrad_multi_field(src, psim, niter, rsqmin, &final_rsq);

    fnorm[0] += fermion_force(f_eps * LAMBDA_MID, src, psim);
#endif
    *gnorm += update_gauge_step(g_eps);

#ifndef PUREGAUGE
    // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
    iters += congrad_multi_field(src, psim, niter, rsqmin, &final_rsq);

    if (i_multi0 < nsteps[0])
      fnorm[0] += fermion_force(f_eps * TWO_LAMBDA, src, psim);
    else
      fnorm[0] += fermion_force(f_eps * LAMBDA, src, psim);
#endif
  }
  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update() {
  int i, iters = 0;
  Real final_rsq, cg_time[2], old_cg_time[2], next_cg_time[2];
  double gnorm = 0.0, fnorm[2] = {0.0, 0.0};
  double startaction = 0.0, endaction, change;
  Twist_Fermion *src, **psim;

  src = malloc(sites_on_node * sizeof(*src));
  psim = malloc(Norder * sizeof(*psim));
  for (i = 0; i < Norder; i++)
    psim[i] = malloc(sites_on_node * sizeof(Twist_Fermion));

  // Refresh the momenta
  // Higher rep code using fermion_rep:
  //   DIMFxDIMF link created from NCOLxNCOL linkf after each update
  ranmom();

  // Set up the fermion variables, if needed
#ifndef PUREGAUGE
  // Compute g and src = (Mdag M)^(1 / 8) g
  iters += grsource(src);

  // Do a CG to get psim,
  // rational approximation to (Mdag M)^(-1 / 4) src = (Mdag M)^(-1 / 8) g
  for (i = 0; i < Norder; i++)
    shift[i] = shift4[i];
#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_o -- original action\n");
#endif
  iters += congrad_multi_field(src, psim, niter, rsqmin, &final_rsq);
#endif // ifndef PUREGAUGE

  // Find initial action
  startaction = d_action(src, psim);

#ifdef HMC_ALGORITHM
  Real xrandom;   // For accept/reject test
  // Copy link field to old_link
  gauge_field_copy_f(F_OFFSET(linkf[0]), F_OFFSET(old_linkf[0]));
#endif
  // Do microcanonical updating
  iters += update_step(old_cg_time, cg_time, next_cg_time,
                       fnorm, &gnorm, src, psim);

  // Find ending action
  // Reuse data from update_step, don't need CG to get (Mdag M)^(-1) chi
  // If the final step were a gauge update, CG would be necessary
  endaction = d_action(src, psim);
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
      gauge_field_copy_f(F_OFFSET(old_linkf[0]), F_OFFSET(linkf[0]));
      fermion_rep();
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

  if (traj_length > 0) {
    node0_printf("IT_PER_TRAJ %d\n", iters);
    node0_printf("MONITOR_FORCE_GAUGE %.4g\n",
                 gnorm / (double)(2 * nsteps[0]));
    node0_printf("MONITOR_FORCE_FERMION0 %.4g\n",
                 fnorm[0] / (double)(2 * nsteps[0]));
    return iters;
  }
  else
    return -99;

  free(src);
  for (i = 0; i < Norder; i++)
    free(psim[i]);
  free(psim);
}
// -----------------------------------------------------------------
