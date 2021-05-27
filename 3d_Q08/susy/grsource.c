// -----------------------------------------------------------------
// Routines for filling fields (momenta and pseudofermions)
// with gaussian random numbers
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Construct gaussian random momentum matrices
// as sum of U(N) generators with gaussian random coefficients
void ranmom() {
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
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Construct a gaussian random vector R
// Return src = (Mdag M)^(1 / 8Nroot) R
// Need to invert despite the positive power, since it is fractional
// Return the number of iterations from the inversion
#ifndef PUREGAUGE
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
#ifdef SITERAND                                 // Volume fermions (for 3d)
      grn.real = gaussian_rand_no(&(s->site_prn));
      grn.imag = gaussian_rand_no(&(s->site_prn));
#else
      grn.real = gaussian_rand_no(&node_prn);
      grn.imag = gaussian_rand_no(&node_prn);
#endif
      c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Fvolume));
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
//      printf("%d %d %d %.4g\n",
//             s->x, s->y, s->t, source_norm);
//      FORALLDIR(mu)
//        dumpmat(&(s->link[mu]));
//    }
//  }

//  fermion_op(src, dest, PLUS);
//  printf("\n\n TEST 2\n");
//  FORALLSITES(i, s) {
//    source_norm = (double)magsq_TF(&(psim[0][i]));
//    if (fabs(source_norm) > IMAG_TOL)
//      printf("%d %d %d %.4g\n",
//             s->x, s->y, s->t, source_norm);
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
// -----------------------------------------------------------------
