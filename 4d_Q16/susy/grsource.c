// -----------------------------------------------------------------
// Routines for filling fields (momenta and pseudofermions)
// with gaussian random numbers
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Construct gaussian random momentum matrices
// as sum of U(N) generators with gaussian random coefficients
#ifdef SMD_ALGORITHM
void smdmom(Real eps) {
// mom -> c1*mom + c2*N(0,1), c1**2 + c2**2 == 1, c1 = exp(-g*eps)
  register int i, j, mu;
  register site *s;
  complex grn;
  Real c1, c2;

  if(friction < 0) //Easy way of getting ranmom for first few trajectories
    c1 = 0;
  else
    c1 = exp(-1 * eps * friction);
  c2 = sqrt(1-c1*c1);

  FORALLSITES(i, s) {
    FORALLDIR(mu) {
      scalar_mult_matrix(&(s->mom[mu]), c1, &(s->mom[mu]));
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        grn.real = c2 * gaussian_rand_no(&(s->site_prn));
        grn.imag = c2 * gaussian_rand_no(&(s->site_prn));
#else
        grn.real = c2 * gaussian_rand_no(&(s->node_prn));
        grn.imag = c2 * gaussian_rand_no(&(s->node_prn));
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(s->mom[mu]));
      }
    }
  }
}
#else
// -----------------------------------------------------------------
void ranmom(Real eps) {
  register int i, j, mu;
  register site *s;
  complex grn;

  FORALLSITES(i, s) {
    FORALLDIR(mu) {
      clear_mat(&(s->mom[mu]));
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&(s->node_prn));
        grn.imag = gaussian_rand_no(&(s->node_prn));
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(s->mom[mu]));
      }
    }
  }
}


#endif
// -----------------------------------------------------------------
// Construct a gaussian random vector R
// Return src = (Mdag M)^(1 / 8Nroot) R
// Need to invert despite the positive power, since it is fractional
// Return the number of iterations from the inversion
#if !defined PUREGAUGE && defined SMD_ALGORITHM
int smdgrsource(Twist_Fermion *src, Real eps) {
  register int i, j, mu;
  register site *s;
  int avs_iters, c1, c2;
  Real size_r;
  complex grn;
  Twist_Fermion **psim = malloc(sizeof(Twist_Fermion*) * Norder);
  Twist_Fermion *old_src = malloc(sizeof(Twist_Fermion) * sites_on_node);

  if(friction < 0) //Easy way of getting old behaviour for first few trajectories
    c1 = 0;
  else
    c1 = exp(-1 * eps * friction);
  c2 = sqrt(1-c1*c1);



  //FORALLSITES(i, s) {
  //  mat_copy(src[i].Fsite, old_src[i].Fsite);
  //  FORALLDIR(mu) mat_copy(src[i].Flink[mu], old_src[i].Flink[mu]);
  //  for (mu = 0; mu < NPLAQ; mu++) 
  //    mat_copy(src[i].Fplaq[mu], old_src[i].Fplaq[mu]);
  //}
  FORALLSITES(i, s) copy_TF(src, old_src);

  // Allocate psim (will be zeroed in congrad_multi)
  for (i = 0; i < Norder; i++)
    psim[i] = malloc(sizeof(Twist_Fermion) * sites_on_node);

  // Begin with pure gaussian random numbers
  FORALLSITES(i, s) {
    clear_TF(&(src[i]));
    for (j = 0; j < DIMF; j++) {                // Site fermions
#ifdef SITERAND
      grn.real = gaussian_rand_no(&(s->site_prn));
      grn.imag = gaussian_rand_no(&(s->site_prn));
#else
      grn.real = gaussian_rand_no(&node_prn);
      grn.imag = gaussian_rand_no(&node_prn);
#endif
      c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Fsite));
      FORALLDIR(mu) {                           // Link fermions
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Flink[mu]));
      }
      for (mu = 0; mu < NPLAQ; mu++) {         // Plaquette fermions
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Fplaq[mu]));
      }
    }
  }

#ifdef DEBUG_CHECK
  double source_norm = 0.0;
  FORALLSITES(i, s) {
//    if (i != 0)
//      clear_TF(&(src[i]));
//    if (i == 0)
//      dump_TF(&src[i]);
    source_norm += (double)magsq_TF(&(src[i]));
  }
  g_doublesum(&source_norm);
  node0_printf("source_norm in grsource %.4g\n", source_norm);

//  DSq(src, dest);
//  printf("\n\n TEST 1\n");
//  FORALLSITES(i, s) {
//    source_norm = (double)magsq_TF(&(psim[1][i]));
//    if (fabs(source_norm) > IMAG_TOL)
//      printf("%d %d %d %d %.4g\n",
//             s->x, s->y, s->z, s->t, source_norm);
//      FORALLDIR(mu)
//        dumpmat(&(s->link[mu]));
//    }
//  }

//  fermion_op(src, dest, PLUS);
//  printf("\n\n TEST 2\n");
//  FORALLSITES(i, s) {
//    source_norm = (double)magsq_TF(&(psim[0][i]));
//    if (fabs(source_norm) > IMAG_TOL)
//      printf("%d %d %d %d %.4g\n",
//             s->x, s->y, s->z, s->t, source_norm);
//      FORALLDIR(mu)
//        dumpmat(&(s->link[mu]));
//    }
//  }
#endif

  // Now compute (Mdag M)^(1 / 8Nroot).src
  for (i = 0; i < Norder; i++)
    shift[i] = shift8[i];

  avs_iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
#ifdef DEBUG_CHECK
  node0_printf("Iters for source %d\n", avs_iters);
#endif
  // Reconstruct (Mdag M)^(1 / 8Nroot).src from multi-mass CG solution psim
  FORALLSITES(i, s) { //c1 and c2 are the SMD coefficients
    scalar_mult_TF(&(src[i]), c2 * ampdeg8, &(src[i]));
    for (j = 0; j < Norder; j++)
      scalar_mult_sum_TF(&(psim[j][i]), c2 * amp8[j], &(src[i]));
    scalar_mult_sum_TF(&(old_src[i]), c1, &(src[i]));
  }

  for (i = 0; i < Norder; i++)
    free(psim[i]);
  free(psim);
  free(old_src);
  return avs_iters;
}
#endif
// -----------------------------------------------------------------
#if !defined PUREGAUGE && !defined SMD_ALGORITHM
int grsource(Twist_Fermion *src) {
  register int i, j, mu;
  register site *s;
  int avs_iters;
  Real size_r;
  complex grn;
  Twist_Fermion **psim = malloc(sizeof(Twist_Fermion*) * Norder);

  // Allocate psim (will be zeroed in congrad_multi)
  for (i = 0; i < Norder; i++)
    psim[i] = malloc(sizeof(Twist_Fermion) * sites_on_node);

  // Begin with pure gaussian random numbers
  FORALLSITES(i, s) {
    clear_TF(&(src[i]));
    for (j = 0; j < DIMF; j++) {                // Site fermions
#ifdef SITERAND
      grn.real = gaussian_rand_no(&(s->site_prn));
      grn.imag = gaussian_rand_no(&(s->site_prn));
#else
      grn.real = gaussian_rand_no(&node_prn);
      grn.imag = gaussian_rand_no(&node_prn);
#endif
      c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Fsite));
      FORALLDIR(mu) {                           // Link fermions
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Flink[mu]));
      }
      for (mu = 0; mu < NPLAQ; mu++) {         // Plaquette fermions
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Fplaq[mu]));
      }
    }
  }

#ifdef DEBUG_CHECK
  double source_norm = 0.0;
  FORALLSITES(i, s) {
//    if (i != 0)
//      clear_TF(&(src[i]));
//    if (i == 0)
//      dump_TF(&src[i]);
    source_norm += (double)magsq_TF(&(src[i]));
  }
  g_doublesum(&source_norm);
  node0_printf("source_norm in grsource %.4g\n", source_norm);

//  DSq(src, dest);
//  printf("\n\n TEST 1\n");
//  FORALLSITES(i, s) {
//    source_norm = (double)magsq_TF(&(psim[1][i]));
//    if (fabs(source_norm) > IMAG_TOL)
//      printf("%d %d %d %d %.4g\n",
//             s->x, s->y, s->z, s->t, source_norm);
//      FORALLDIR(mu)
//        dumpmat(&(s->link[mu]));
//    }
//  }

//  fermion_op(src, dest, PLUS);
//  printf("\n\n TEST 2\n");
//  FORALLSITES(i, s) {
//    source_norm = (double)magsq_TF(&(psim[0][i]));
//    if (fabs(source_norm) > IMAG_TOL)
//      printf("%d %d %d %d %.4g\n",
//             s->x, s->y, s->z, s->t, source_norm);
//      FORALLDIR(mu)
//        dumpmat(&(s->link[mu]));
//    }
//  }
#endif

  // Now compute (Mdag M)^(1 / 8Nroot).src
  for (i = 0; i < Norder; i++)
    shift[i] = shift8[i];

  avs_iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
#ifdef DEBUG_CHECK
  node0_printf("Iters for source %d\n", avs_iters);
#endif
  // Reconstruct (Mdag M)^(1 / 8Nroot).src from multi-mass CG solution psim
  FORALLSITES(i, s) {
    scalar_mult_TF(&(src[i]), ampdeg8, &(src[i]));
    for (j = 0; j < Norder; j++)
      scalar_mult_sum_TF(&(psim[j][i]), amp8[j], &(src[i]));
  }

  for (i = 0; i < Norder; i++)
    free(psim[i]);
  free(psim);
  return avs_iters;
}
#endif
