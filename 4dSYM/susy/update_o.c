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
    for (mu = XUP; mu < NUMLINK; mu++) {
      clear_su3mat_f(&(s->mom[mu]));
      for (j = 0; j < DIMF; j++) {
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

  for (mu = XUP; mu < NUMLINK; mu++) {
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
  Real final_rsq, f_eps, g_eps, tr;

  f_eps = traj_length / (Real)nsteps[0];
  g_eps = f_eps / (Real)(2.0 * nsteps[1]);

#ifndef PUREGAUGE
  tr = fermion_force(f_eps * LAMBDA, src, psim);
  fnorm[0] += tr;
  if (tr > max_ff[0])
    max_ff[0] = tr;
#endif

  for (i_multi0 = 1; i_multi0 <= nsteps[0]; i_multi0++) {
    tr = update_gauge_step(g_eps);
    *gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;

#ifndef PUREGAUGE
    // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
    iters += congrad_multi_field(src, psim, niter, rsqmin, &final_rsq);

    tr = fermion_force(f_eps * LAMBDA_MID, src, psim);
    fnorm[0] += tr;
    if (tr > max_ff[0])
      max_ff[0] = tr;
#endif
    tr = update_gauge_step(g_eps);
    *gnorm += tr;
    if (tr > max_gf)
      max_gf = tr;

#ifndef PUREGAUGE
    // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
    iters += congrad_multi_field(src, psim, niter, rsqmin, &final_rsq);

    if (i_multi0 < nsteps[0])
      tr = fermion_force(f_eps * TWO_LAMBDA, src, psim);
    else
      tr = fermion_force(f_eps * LAMBDA, src, psim);
    fnorm[0] += tr;
    if (tr > max_ff[0])
      max_ff[0] = tr;
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
  // congrad_multi_field initializes psim and calls fermion_rep()
  iters += congrad_multi_field(src, psim, niter, rsqmin, &final_rsq);
#endif // ifndef PUREGAUGE

  // Find initial action
  startaction = d_action(src, psim);

#if 0
  // Test gauge invariance of action
  // Gauge fixing doesn't include fermions...
//  gaugefix(TUP, 1.5, 5000, GAUGE_FIX_TOL, -1, -1);
  // Instead apply random gauge transformation at a single site
  su3_matrix_f Gmat, tmat;
  su3_matrix_f etamat, psimat[NUMLINK], chimat[NPLAQ];
  int a, b, j, x = 1, y = 1, z = 1, t = 1, s = node_index(x, y, z, t);
  complex tc;
  polar(&(lattice[10].mom[2]), &Gmat);

  // Confirm unitarity or check invariance when Gmat = I
//  mult_su3_na_f(&Gmat, &Gmat, &tmat);
//  dumpmat_f(&tmat);
//  su3mat_copy_f(&tmat, &Gmat);

  // Left side of local eta
  clear_su3mat_f(&etamat);
  // Construct G eta = sum_j eta^j G Lambda^j
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fsite.c[j];
    c_scalar_mult_add_su3mat_f(&etamat, &tmat, &tc, &etamat);
  }
  // Project out eta^j = -Tr[Lambda^j G eta]
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &etamat, &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fsite.c[j]);
  }

  // Right side of local eta
  clear_su3mat_f(&etamat);
  // Construct eta Gdag = sum_j eta^j Lambda^j Gdag
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = src[s].Fsite.c[j];
    c_scalar_mult_add_su3mat_f(&etamat, &tmat, &tc, &etamat);
  }
  // Project out eta^j = -Tr[eta Gdag Lambda^j]
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&etamat, &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fsite.c[j]);
  }

  // Left side of local links and psis; right side of local chis
  for (a = XUP; a < NUMLINK; a++) {
    mult_su3_nn_f(&Gmat, &(lattice[s].linkf[a]), &tmat);
    su3mat_copy_f(&tmat, &(lattice[s].linkf[a]));

    clear_su3mat_f(&(psimat[a]));
    for (j = 0; j < DIMF; j++) {
      mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
      tc = src[s].Flink[a].c[j];
      c_scalar_mult_add_su3mat_f(&(psimat[a]), &tmat, &tc, &(psimat[a]));
    }
    for (j = 0; j < DIMF; j++) {
      mult_su3_nn_f(&(Lambda[j]), &(psimat[a]), &tmat);
      tc = trace_su3_f(&tmat);
      CNEGATE(tc, src[s].Flink[a].c[j]);
    }

    for (b = a + 1; b < NUMLINK; b++) {
      i = plaq_index[a][b];
      clear_su3mat_f(&(chimat[i]));
      for (j = 0; j < DIMF; j++) {
        mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
        tc = src[s].Fplaq[i].c[j];
        c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
      }
      for (j = 0; j < DIMF; j++) {
        mult_su3_nn_f(&(chimat[i]), &(Lambda[j]), &tmat);
        tc = trace_su3_f(&tmat);
        CNEGATE(tc, src[s].Fplaq[i].c[j]);
      }
    }
  }
  // Right side of neighboring links and psis
  s = node_index(x - 1, y, z, t);
  mult_su3_na_f(&(lattice[s].linkf[0]), &Gmat, &tmat);
  su3mat_copy_f(&tmat, &(lattice[s].linkf[0]));
  clear_su3mat_f(&(psimat[0]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = src[s].Flink[0].c[j];
    c_scalar_mult_add_su3mat_f(&(psimat[0]), &tmat, &tc, &(psimat[0]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(psimat[0]), &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Flink[0].c[j]);
  }

  s = node_index(x, y - 1, z, t);
  mult_su3_na_f(&(lattice[s].linkf[1]), &Gmat, &tmat);
  su3mat_copy_f(&tmat, &(lattice[s].linkf[1]));
  clear_su3mat_f(&(psimat[1]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = src[s].Flink[1].c[j];
    c_scalar_mult_add_su3mat_f(&(psimat[1]), &tmat, &tc, &(psimat[1]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(psimat[1]), &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Flink[1].c[j]);
  }

  s = node_index(x, y, z - 1, t);
  mult_su3_na_f(&(lattice[s].linkf[2]), &Gmat, &tmat);
  su3mat_copy_f(&tmat, &(lattice[s].linkf[2]));
  clear_su3mat_f(&(psimat[2]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = src[s].Flink[2].c[j];
    c_scalar_mult_add_su3mat_f(&(psimat[2]), &tmat, &tc, &(psimat[2]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(psimat[2]), &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Flink[2].c[j]);
  }

  s = node_index(x, y, z, t - 1);
  mult_su3_na_f(&(lattice[s].linkf[3]), &Gmat, &tmat);
  su3mat_copy_f(&tmat, &(lattice[s].linkf[3]));
  clear_su3mat_f(&(psimat[3]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = src[s].Flink[3].c[j];
    c_scalar_mult_add_su3mat_f(&(psimat[3]), &tmat, &tc, &(psimat[3]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(psimat[3]), &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Flink[3].c[j]);
  }

  s = node_index(x + 1, y + 1, z + 1, t + 1);
  mult_su3_na_f(&(lattice[s].linkf[4]), &Gmat, &tmat);
  su3mat_copy_f(&tmat, &(lattice[s].linkf[4]));
  clear_su3mat_f(&(psimat[4]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = src[s].Flink[4].c[j];
    c_scalar_mult_add_su3mat_f(&(psimat[4]), &tmat, &tc, &(psimat[4]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(psimat[4]), &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Flink[4].c[j]);
  }

  // Left side of neighboring chis
  s = node_index(x - 1, y - 1, z, t);       // 01
  i = plaq_index[0][1];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fplaq[i].c[j]);
  }

  s = node_index(x - 1, y, z - 1, t);       // 02
  i = plaq_index[0][2];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fplaq[i].c[j]);
  }

  s = node_index(x - 1, y, z, t - 1);       // 03
  i = plaq_index[0][3];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fplaq[i].c[j]);
  }

  s = node_index(x, y + 1, z + 1, t + 1);   // 04
  i = plaq_index[0][4];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fplaq[i].c[j]);
  }

  s = node_index(x, y - 1, z - 1, t);       // 12
  i = plaq_index[1][2];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fplaq[i].c[j]);
  }

  s = node_index(x, y - 1, z, t - 1);       // 13
  i = plaq_index[1][3];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fplaq[i].c[j]);
  }

  s = node_index(x + 1, y, z + 1, t + 1);   // 14
  i = plaq_index[1][4];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fplaq[i].c[j]);
  }

  s = node_index(x, y, z - 1, t - 1);       // 23
  i = plaq_index[2][3];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fplaq[i].c[j]);
  }

  s = node_index(x + 1, y + 1, z, t + 1);   // 24
  i = plaq_index[2][4];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fplaq[i].c[j]);
  }

  s = node_index(x + 1, y + 1, z + 1, t);   // 34
  i = plaq_index[3][4];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = src[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, src[s].Fplaq[i].c[j]);
  }

  // Now ready to check!
  iters += congrad_multi_field(src, psim, niter, rsqmin, &final_rsq);
  startaction = d_action(src, psim);
  terminate(1);
#endif
  max_gf = 0.0;
  max_ff[0] = 0.0;

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
    node0_printf("MONITOR_FORCE_GAUGE    %.4g %.4g\n",
                 gnorm / (double)(2 * nsteps[0]), max_gf);
    node0_printf("MONITOR_FORCE_FERMION0 %.4g %.4g\n",
                 fnorm[0] / (double)(2 * nsteps[0]), max_ff[0]);
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
