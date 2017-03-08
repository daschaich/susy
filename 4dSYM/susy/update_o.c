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
#ifdef HYBRID

void update_uu(Real eps) {
    register int i, mu;
    register site *s;
    register Real t2, t3, t4, t5, t6, t7, t8;
    matrix tmat, tmat2;
    
    
    // -----------------------------------------------------------------
    // Calculate newU = exp(p).U
    // Here p is the traceless BUT NOT anti-hermitian lattice field.
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
    
    FORALLSITES(i, s){
        
        FORALLDIR(mu) {
            
            
            mult_nn(&(s->mom[mu]),&(s->link[mu]),&tmat);
            scalar_mult_add_matrix(&(s->link[mu]), &tmat, t8, &tmat2);
            
            mult_nn(&(s->mom[mu]),&tmat2,&tmat);
            scalar_mult_add_matrix(&(s->link[mu]), &tmat, t7, &tmat2);
            
            mult_nn(&(s->mom[mu]),&tmat2,&tmat);
            scalar_mult_add_matrix(&(s->link[mu]), &tmat, t6, &tmat2);
            
            mult_nn(&(s->mom[mu]),&tmat2,&tmat);
            scalar_mult_add_matrix(&(s->link[mu]), &tmat, t5, &tmat2);
            
            mult_nn(&(s->mom[mu]),&tmat2,&tmat);
            scalar_mult_add_matrix(&(s->link[mu]), &tmat, t4, &tmat2);
            
            mult_nn(&(s->mom[mu]),&tmat2,&tmat);
            scalar_mult_add_matrix(&(s->link[mu]), &tmat, t3, &tmat2);
            
            mult_nn(&(s->mom[mu]),&tmat2,&tmat);
            scalar_mult_add_matrix(&(s->link[mu]), &tmat, t2, &tmat2);
            
            
            mult_nn(&(s->mom[mu]),&tmat2,&tmat);
            scalar_mult_sum_matrix(&tmat, eps, &(s->link[mu]));
            
        }
        
    }
    
    // Update plaquette determinants, DmuUmu and Fmunu with new links
    // (Needs to be done before calling gauge_force)
    compute_plaqdet();
    compute_Uinv();
    compute_DmuUmu();
    compute_Fmunu();
}


#else

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

#endif 

// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Omelyan version; ``dirty'' speeded-up version
double update_gauge_step(Real eps) {
  int nsw = nsteps[1], isw;
  double norm;

#ifdef UPDATE_DEBUG
  node0_printf("gauge %d steps %.4g dt\n", nsw, eps);
#endif
  norm = gauge_force(eps * LAMBDA);
  for (isw = 1; isw <= nsw; isw++) {
    update_uu(0.5 * eps);
    norm += gauge_force(eps * LAMBDA_MID);
    update_uu(0.5 * eps);
    if (isw < nsw)
      norm += gauge_force(eps * TWO_LAMBDA);

    else
      norm += gauge_force(eps * LAMBDA);
  }
  return (norm / nsw);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update_step(double *fnorm, double *gnorm,
                Twist_Fermion **src, Twist_Fermion ***psim) {

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
    *gnorm += tr;
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
    tr = update_gauge_step(g_eps);
    *gnorm += tr;
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
  node0_printf("Calling CG in update_o -- original action\n");
#endif
  // congrad_multi initializes psim
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
#endif // ifndef PUREGAUGE

  // Find initial action
  startaction = action(src, psim);
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
  gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));
#endif
  // Do microcanonical updating
  iters += update_step(fnorm, &gnorm, src, psim);

  // Find ending action
  // Reuse data from update_step, don't need CG to get (Mdag M)^(-1) chi
  // If the final step were a gauge update, CG would be necessary
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
