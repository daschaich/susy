// -----------------------------------------------------------------
// Update the lattice using a leapfrog algorithm beginning at integral time,
// with H and U evaluated at the same time
// nsteps[num_masses+1] are defined in lattice.h and set in setup.c
#include "susy_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>    /* For "finite" */
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
        grn.real= gaussian_rand_no(&(s->site_prn));
        grn.imag= gaussian_rand_no(&(s->site_prn));
#else
        grn.real= gaussian_rand_no(&(s->node_prn));
        grn.imag= gaussian_rand_no(&(s->node_prn));
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
int update_step(Real *old_cg_time, Real *cg_time, Real *next_cg_time,
                double *fnorm, double *gnorm,
                Twist_Fermion *source, Twist_Fermion **psim) {

  int step, iters = 0;
  Real eps;
  Real final_rsq;
  void update_uu(Real eps);

  eps = traj_length / (Real)nsteps[0];

  node0_printf("eps %e\n",eps);

  for (step = 1; step <= nsteps[0]; step++) {
/*
// one step p(t/2) u(t) p(t/2)
    *gnorm += gauge_force(0.5 * eps);
    // Do conjugate gradient to get (Madj M)inverse * chi
    iters += congrad_multi_field( source,psim, niter,rsqmin, &final_rsq);
    fnorm[0] += fermion_force(0.5*eps,source,psim);

    update_uu(eps);
    fermion_rep();

    *gnorm += gauge_force(0.5 * eps);
    // Do conjugate gradient to get (Madj M)inverse * chi
    iters += congrad_multi_field(source, psim, niter, rsqmin, &final_rsq);
    fnorm[0] += fermion_force(0.5 * eps, source, psim);
*/


// one step u(t/2) p(t) u(t/2)
    update_uu(0.5 * eps);
    fermion_rep();

    *gnorm += gauge_force(eps);
    // Do conjugate gradient to get (Madj M)inverse * chi
    iters += congrad_multi_field(source, psim, niter, rsqmin, &final_rsq);
    fnorm[0] += fermion_force(eps, source, psim);
    update_uu(0.5 * eps);
    fermion_rep();
  }
  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update() {
  int i;
  register site *s;
  int iters=0;
  double source_norm;
  Real final_rsq, Real cg_time[2], old_cg_time[2], next_cg_time[2];

  Twist_Fermion *source;
  Twist_Fermion **psim;

  source = (Twist_Fermion *)malloc(sites_on_node * sizeof(Twist_Fermion));
  psim = (Twist_Fermion **)malloc(Norder * sizeof(Twist_Fermion *));
  for (i = 0; i < Norder; i++)
    psim[i] = (Twist_Fermion *)malloc(sites_on_node * sizeof(Twist_Fermion));

//#ifdef HMC_ALGORITHM
  double startaction = 0, endaction,change;
  Real xrandom;
//#endif

  double gnorm=0.0;
  double fnorm[2];
  fnorm[0]=fnorm[1]=0.0;

  /* refresh the momenta */
  ranmom();
  /* higher rep code:
     DIMFxDIMF link created from NCOLxNCOL linkf after each update,
     then DIMF gauge field is switched to antiperiodic b.c. if required
  */

  fermion_rep();
  grsource(source);
  // Do conjugate gradient to get psim=components of RHMC * source
  for (i = 0; i < Norder; i++)
    shift[i] = shift4[i];
  //node0_printf("Calling CG in update_o --original action\n");
  iters += congrad_multi_field(source, psim, niter, rsqmin, &final_rsq);

//#ifdef HMC_ALGORITHM
  startaction = d_action(source, psim);   // Find action
//  node0_printf("startaction= %.4g\n",startaction);
#ifdef HMC_ALGORITHM
  // Copy link field to old_link
  gauge_field_copy_f(F_OFFSET(linkf[0]), F_OFFSET(old_linkf[0]));
#endif
  // Do microcanonical updating
  iters += update_step(old_cg_time, cg_time, next_cg_time,
                       fnorm, &gnorm, source, psim);

  /* do conjugate gradient to get (Madj M)inverse * chi. */
  fermion_rep();
  iters += congrad_multi_field(source, psim, niter, rsqmin, &final_rsq);
  endaction = d_action(source, psim);     // Find action
  change = endaction - startaction;
  // Reject configurations giving overflow
#ifdef HMC_ALGORITHM
#ifndef HAVE_IEEEFP_H
  if (fabs((double)change)>1e20) {
#else
  if (!finite((double)change)) {
#endif
    node0_printf("WARNING: Correcting Apparent Overflow: Delta S = %e\n", change);
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
    node0_printf("REJECT: delta S = %.4g, start S = %.12g, end S = %.12g\n",
                 change, startaction, endaction);
  }
  else {
    node0_printf("ACCEPT: delta S = %.4g, start S = %.12g, end S = %.12g\n",
                 change, startaction, endaction);
  }
#else
  // Only print check if not doing HMC
  node0_printf("CHECK: delta S = %e\n", (double)(change));
#endif // HMC

  if (traj_length > 0) {
    node0_printf("IT_PER_TRAJ %d\n", iters );
    node0_printf("MONITOR_FORCE_GAUGE %e\n",gnorm/(double)(2*nsteps[0]) );
    node0_printf("MONITOR_FORCE_FERMION0 %e\n",fnorm[0]/(double)(2*nsteps[0]) );
    return iters;
  }

  else return(-99);

  free(source);
  for (i = 0; i < Norder; i++)
    free(psim[i]);
  free(psim);
}
// -----------------------------------------------------------------
