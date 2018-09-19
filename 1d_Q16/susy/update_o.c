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
  register int i, j;
  register site *s;
  register Real t2, t3, t4, t5, t6, t7, t8;
  matrix tmat, tmat2, tmp_mom;

  // Calculate newU = exp(p).U
  // Go to eighth order in the exponential:
  //   exp(p) * U = (1 + p + p^2/2 + p^3/6 ...) * U
  //              = U + p*(U + (p/2)*(U + (p/3)*( ... )))
  // Take divisions out of site loop (can't be done by compiler)
  t2 = eps / 2.0;
  t3 = eps / 3.0;
  t4 = eps / 4.0;
  t5 = eps / 5.0;
  t6 = eps / 6.0;
  t7 = eps / 7.0;
  t8 = eps / 8.0;

  FORALLSITES(i, s) {
    uncompress_anti_hermitian(&(s->mom), &tmp_mom);
    mult_nn(&tmp_mom, &(s->link), &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t8, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t7, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t6, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t5, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t4, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t3, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t2, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_sum_matrix(&tmat, eps, &(s->link));

    for (j = 0; j < NSCALAR; j++)
      scalar_mult_sum_matrix(&(s->mom_X[j]), eps, &(s->X[j]));
  }
  // Update dot product of first eight scalar fields and gamma matrices
  build_Gamma_X();
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Omelyan version; ``dirty'' speeded-up version
double update_bosonic_step(Real eps) {
  int n = nsteps[1], i;
  double norm;
#ifdef UPDATE_DEBUG
  double td, td2;
#endif

#ifdef UPDATE_DEBUG
  node0_printf("gauge %d steps %.4g dt\n", n, eps);
#endif
  norm = bosonic_force(eps * LAMBDA);
  for (i = 1; i <= n; i++) {
    update_u(0.5 * eps);
    norm += bosonic_force(eps * LAMBDA_MID);
    update_u(0.5 * eps);
    if (i < n)
      norm += bosonic_force(eps * TWO_LAMBDA);

    else
      norm += bosonic_force(eps * LAMBDA);
  }

  // Reunitarize the gauge field and re-anti-hermitianize the scalars
#ifdef UPDATE_DEBUG
  td = check_unitarity();
  g_floatmax(&td);
#endif
  reunitarize();
  reantihermize();
#ifdef UPDATE_DEBUG
  td2 = check_unitarity();
  g_floatmax(&td2);
  node0_printf("Reunitarized after boson update step.  ");
  node0_printf("Max deviation %.2g changed to %.2g\n", td, td2);
#endif

  return (norm / n);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update_step(matrix ***src, matrix ****psim) {
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
#ifdef UPDATE_DEBUG
    node0_printf("Step %d - 1 Action %.4g Force %.4g\n", i_multi0,
                 action(src, psim), tr);
#endif
    bnorm += tr;
    if (tr > max_bf)
      max_bf = tr;

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
#ifdef UPDATE_DEBUG
    node0_printf("Step %d - 2 Action %.4g Force %.4g\n", i_multi0,
                 action(src, psim), tr);
#endif
    bnorm += tr;
    if (tr > max_bf)
      max_bf = tr;

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
  int n, iters = 0;
  double startaction, endaction, change;
  matrix ***source = malloc(sizeof(matrix**) * Nroot);
  matrix ****psim = malloc(sizeof(matrix***) * Nroot);

  // Refresh the momenta
  ranmom();

  // Set up the fermion variables, if needed
#ifndef PUREGAUGE
  register int i, j, k;
  register site *s;
  Real final_rsq;

  for (n = 0; n < Nroot; n++) {
    psim[n] = malloc(sizeof(matrix**) * Norder);
    for (j = 0; j < Norder; j++) {
      psim[n][j] = malloc(sizeof(matrix*) * NFERMION);
      for (k = 0; k < NFERMION; k++)
        psim[n][j][k] = malloc(sizeof(matrix) * sites_on_node);
    }

    source[n] = malloc(sizeof(matrix*) * NFERMION);
    for (k = 0; k < NFERMION; k++)
      source[n][k] = malloc(sizeof(matrix) * sites_on_node);
  }

  // Compute source = (Mdag M)^(1 / 8) g
  for (n = 0; n < Nroot; n++)
    iters += grsource(source[n]);

  // Do a CG to get psim,
  // rational approximation to (Mdag M)^(-1 / 4) src = (Mdag M)^(-1 / 8) g
  for (j = 0; j < Norder; j++)
    shift[j] = shift4[j];
#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_o -- original action\n");
#endif
  // congrad_multi initializes psim
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi(source[n], psim[n], niter, rsqmin, &final_rsq);
#endif // ifndef PUREGAUGE

  // Find initial action
  startaction = action(source, psim);
  bnorm = 0.0;
  max_bf = 0.0;
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
    if (traj_length > 0.0) {
      copy_bosons(MINUS);   // Restore link field from old_link
      build_Gamma_X();
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

#ifndef PUREGAUGE
  for (n = 0; n < Nroot; n++) {
    for (j = 0; j < Norder; j++) {
      for (k = 0; k < NFERMION; k++)
        free(psim[n][j][k]);
      free(psim[n][j]);
    }
    for (k = 0; k < NFERMION; k++)
      free(source[n][k]);
    free(source[n]);
    free(psim[n]);
  }
#endif
  free(source);
  free(psim);

  if (traj_length > 0) {
    node0_printf("IT_PER_TRAJ %d\n", iters);
    node0_printf("MONITOR_FORCE_GAUGE    %.4g %.4g\n",
        bnorm / (double)(2.0 * nsteps[0]), max_bf);
    for (n = 0; n < Nroot; n++) {
      node0_printf("MONITOR_FORCE_FERMION%d %.4g %.4g\n",
          n, fnorm[n] / (double)(2.0 * nsteps[0]), max_ff[n]);
    }
    return iters;
  }
  else
    return -99;
}
// -----------------------------------------------------------------
