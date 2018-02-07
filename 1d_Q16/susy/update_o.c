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
void update_uu(Real eps) {
  register int i, j;
  register site *s;

  FORALLSITES(i, s) {
    scalar_mult_sum_matrix(&(s->mom), eps, &(s->link));
    for(j=0;j<NSCALAR;j++)
      scalar_mult_sum_matrix(&(s->mom_X[j]), eps, &(s->X[j]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Omelyan version; ``dirty'' speeded-up version
double update_bosonic_step(Real eps) {
  int nsw = nsteps[1], isw;
  double norm;

#ifdef UPDATE_DEBUG
  node0_printf("gauge %d steps %.4g dt\n", nsw, eps);
#endif
  norm = bosonic_force(eps * LAMBDA);
  for (isw = 1; isw <= nsw; isw++) {
    update_uu(0.5 * eps);
    norm += bosonic_force(eps * LAMBDA_MID);
    update_uu(0.5 * eps);
    if (isw < nsw)
      norm += bosonic_force(eps * TWO_LAMBDA);

    else
      norm += bosonic_force(eps * LAMBDA);
  }
  return (norm / nsw);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update_step(matrix **src[NFERMION], matrix ***psim[NFERMION]) {
  int iters = 0, i_multi0;
  Real f_eps, b_eps, tr;

  f_eps = traj_length / (Real)nsteps[0];
  b_eps = f_eps / (Real)(2.0 * nsteps[1]);

#ifndef PUREGAUGE
  int n;
  Real final_rsq;
  for (n = 0; n < Nroot; n++) {
    // CG called before update_step
    tr = fermion_force(f_eps * LAMBDA, src[n], psim[n]);
    fnorm[n] += tr;
    if (tr > max_ff[n])
      max_ff[n] = tr;
  }
#endif

  for (i_multi0 = 1; i_multi0 <= nsteps[0]; i_multi0++) {
    tr = update_bosonic_step(b_eps);
    gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;

#ifndef PUREGAUGE
    for (n = 0; n < Nroot; n++) {
      // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
      iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
      tr = fermion_force(f_eps * LAMBDA_MID, src[n], psim[n]);
      fnorm[n] += tr;
      if (tr > max_ff[n])
        max_ff[n] = tr;
    }
#endif
    tr = update_bosonic_step(b_eps);
    gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;

#ifndef PUREGAUGE
    for (n = 0; n < Nroot; n++) {
      // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
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
int update() {
  int j, k, n, iters = 0;
  double startaction, endaction, change;
  matrix **source[NFERMION], ***psim[NFERMION];
  for(k=0; k<NFERMION; k++) {
    source[k] = malloc(sizeof(matrix*) * Nroot);
    psim[k] = malloc(sizeof(matrix**) * Nroot);
    
    for (n = 0; n < Nroot; n++) {
      source[k][n] = malloc(sizeof(matrix) * sites_on_node);
      psim[k][n] = malloc(sizeof(matrix*) * Norder);
      for (j = 0; j < Norder; j++)
        psim[k][n][j] = malloc(sizeof(matrix) * sites_on_node);
    }
  }

  // Refresh the momenta
  ranmom();

  // Set up the fermion variables, if needed
#ifndef PUREGAUGE
  Real final_rsq;
  // Compute g and src = (Mdag M)^(1 / 8) g
  for (n = 0; n < Nroot; n++) {
    iters += grsource(src);
    for(j=0; j<NFERMION; j++) {
      FORALLSITES(i,s)
      mat_copy(&src[j][i], &source[j][n][i]);
    }
  }
  FIX THIS

  // Do a CG to get psim,
  // rational approximation to (Mdag M)^(-1 / 4) src = (Mdag M)^(-1 / 8) g
  for (j = 0; j < Norder; j++)
    shift[j] = shift4[j];
#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_o -- original action\n");
#endif
  // congrad_multi initializes psim
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
#endif // ifndef PUREGAUGE

  // Find initial action
  startaction = action(source, psim);
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
//    random_gauge_trans(src[n]);
//    congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
//  }
//  startaction = action(src, psim);
//  node0_printf("AFTER  GTRANS %.8g\n", startaction);
//  terminate(1);

#ifdef HMC_ALGORITHM
  Real xrandom;   // For accept/reject test
  // Copy link field to old_link
  copy_bosons(PLUS);
#endif
  // Do microcanonical updating
  iters += update_step(source, psim);

  // Find ending action
  // Reuse data from update_step, don't need CG to get (Mdag M)^(-1 / 4) chi
  // If the final step were a gauge update, CG would be necessary
  endaction = action(source, psim);
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
    if (traj_length > 0.0)
      copy_bosons(MINUS);
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
