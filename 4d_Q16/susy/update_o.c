// -----------------------------------------------------------------
// Update lattice
// Omelyan integrator multiscale following CPC 174:87 (2006)

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
void update_u(Real eps) {
  register int i, mu;
  register site *s;

  FORALLSITES(i, s) {
    FORALLDIR(mu)
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
// Omelyan version; ``dirty'' speeded-up version
double update_gauge_step(Real eps) {
  int n = nsteps[1], i;
  double norm;

#ifdef UPDATE_DEBUG
  node0_printf("gauge %d steps %.4g dt\n", n, eps);
#endif
  norm = gauge_force(eps * LAMBDA);
  for (i = 1; i <= n; i++) {
    update_u(0.5 * eps);
    norm += gauge_force(eps * LAMBDA_MID);
    update_u(0.5 * eps);
    if (i < n)
      norm += gauge_force(eps * TWO_LAMBDA);

    else
      norm += gauge_force(eps * LAMBDA);
  }
  return (norm / n);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update_step(Twist_Fermion **src, Twist_Fermion ***psim) {
  int iters = 0, i_multi0, n;
  Real final_rsq, f_eps, g_eps, tr;

  f_eps = traj_length / (Real)nsteps[0];
  g_eps = f_eps / (Real)(2.0 * nsteps[1]);

#ifndef PUREGAUGE
  for (n = 0; n < Nroot; n++) {
    // CG called before update_step
    tr = fermion_force(f_eps * LAMBDA, src[n], psim[n]);
    fnorm[n] += tr;
    if (tr > max_ff[n])
      max_ff[n] = tr;
  }
#endif

  for (i_multi0 = 1; i_multi0 <= nsteps[0]; i_multi0++) {
    tr = update_gauge_step(g_eps);
    gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;

#ifndef PUREGAUGE
    for (n = 0; n < Nroot; n++) {
      // Do conjugate gradient to get (Mdag M)^(-1 / 4Nroot) chi
      iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
      tr = fermion_force(f_eps * LAMBDA_MID, src[n], psim[n]);
      fnorm[n] += tr;
      if (tr > max_ff[n])
        max_ff[n] = tr;
    }
#endif
    tr = update_gauge_step(g_eps);
    gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;

#ifndef PUREGAUGE
    for (n = 0; n < Nroot; n++) {
      // Do conjugate gradient to get (Mdag M)^(-1 / 4Nroot) chi
      iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);

      if (i_multi0 < nsteps[0])
        tr = fermion_force(f_eps * TWO_LAMBDA, src[n], psim[n]);
      else
        tr = fermion_force(f_eps * LAMBDA, src[n], psim[n]);
      fnorm[n] += tr;
      if (tr > max_ff[n])
        max_ff[n] = tr;
    }
#endif
  }
  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update() {//src -> gsrc as src is already defined for EIG,etc.
              //Will check if possible to go back
  int j, n, iters = 0;
  Real final_rsq;
  double startaction, endaction, change;
#ifndef SMD_ALGORITHM
  Twist_Fermion **gsrc = malloc(sizeof(Twist_Fermion*) * Nroot);
#endif
  Twist_Fermion ***psim = malloc(sizeof(Twist_Fermion**) * Nroot);

  register int i,mu;
  register site *s; //For the accept/reject step

  for (n = 0; n < Nroot; n++) {
#ifndef SMD_ALGORITHM
    gsrc[n] = malloc(sizeof(Twist_Fermion) * sites_on_node);
#endif
    psim[n] = malloc(sizeof(Twist_Fermion*) * Norder);
    for (j = 0; j < Norder; j++)
      psim[n][j] = malloc(sizeof(Twist_Fermion) * sites_on_node);
  }

#ifdef SMD_ALGORITHM
  // Refresh the momenta
  smdmom(traj_length);
#else
  ranmom();
#endif

  // Set up the fermion variables, if needed
#ifndef PUREGAUGE
  // Compute g and gsrc = (Mdag M)^(1 / 8Nroot) g
  for (n = 0; n < Nroot; n++) {
#ifdef SMD_ALGORITHM
    iters += smdgrsource(gsrc[n], traj_length);
#else
    iters += grsource(gsrc[n]);
#endif
  }

  // Do a CG to get psim, rational approximation
  // to (Mdag M)^(-1 / 4Nroot) gsrc = (Mdag M)^(-1 / 8Nroot) g
  for (j = 0; j < Norder; j++)
    shift[j] = shift4[j];
#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_o -- original action\n");
#endif
  // congrad_multi initializes psim
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi(gsrc[n], psim[n], niter, rsqmin, &final_rsq);
#endif // ifndef PUREGAUGE

  // Find initial action
  startaction = action(gsrc, psim);
  gnorm = 0.0;
  max_gf = 0.0;
  for (n = 0; n < Nroot; n++) {
    fnorm[n] = 0.0;
    max_ff[n] = 0.0;
  }

  // Uncomment this block to test gauge invariance of action
  // by re-measuring after applying a random gauge transformation
  // at a single site in a lattice with at least L=4 in all directions
//  node0_printf("BEFORE GTRANS %.8g\n", startaction);
//  for (n = 0; n < Nroot; n++) {
//    random_gauge_trans(gsrc[n]);
//    congrad_multi(gsrc[n], psim[n], niter, rsqmin, &final_rsq);
//  }
//  startaction = action(gsrc, psim);
//  node0_printf("AFTER  GTRANS %.8g\n", startaction);
//  terminate(1);

#if defined HMC_ALGORITHM || defined SMD_ALGORITHM
  Real xrandom;   // For accept/reject test
  // Copy link field to old_link
  gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));
  FORALLSITES(i, s)
    FORALLDIR(mu)
      mat_copy(&(s->mom[mu]), &(s->old_mom[mu]));
#endif
  // Do microcanonical updating
  iters += update_step(gsrc, psim);

  // Find ending action
  // Reuse (Mdag M)^(-1 / 4Nroot) chi from update_step, don't need CG
  // If the final step were a gauge update, CG would be necessary
  endaction = action(gsrc, psim);
  change = endaction - startaction;
#if defined HMC_ALGORITHM || defined SMD_ALGORITHM
  // Reject configurations giving overflow
#ifndef HAVE_IEEEFP_H
  if (fabs((double)change) > 1e20) {
#else
  if (!finite((double)change)) {
#endif
    node0_printf("WARNING: Correcting `Apparent Overflow: Delta S = %.4g\n",
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
      FORALLSITES(i,s)
        FORALLDIR(mu)
          scalar_mult_matrix(&(s->old_mom[mu]), -1, &(s->mom[mu]));//Reverse momentum
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
  // Only print check if not doing HMC or SMD
  node0_printf("CHECK: delta S = %.4g\n", (double)(change));
#endif // ifdef HMC || SMD

  for (n = 0; n < Nroot; n++) {
#ifndef SMD_ALGORITHM
    free(gsrc[n]);
#endif
    for (j = 0; j < Norder; j++)
      free(psim[n][j]);
    free(psim[n]);
  }
#ifndef SMD_ALGORITHM
  free(gsrc);
#endif
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
