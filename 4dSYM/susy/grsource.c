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
      clear_mat_f(&(s->mom[mu]));
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&(s->node_prn));
        grn.imag = gaussian_rand_no(&(s->node_prn));
#endif
        c_scalar_mult_sum_mat_f(&(Lambda[j]), &grn, &(s->mom[mu]));
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Construct a gaussian random vector R, return src = (Mdag M)^{1 / 8} R
// Need to invert despite the positive power, since it is fractional
// Return the number of iterations from the inversion
int grsource(Twist_Fermion *src) {
  register int i, j, mu;
  register site *s;
  int avs_iters;
  Real size_r;
  Twist_Fermion **psim = malloc(Norder * sizeof(**psim));

  // Allocate psim (will be zeroed in congrad_multi_field)
  for (i = 0; i < Norder; i++)
    psim[i] = malloc(sites_on_node * sizeof(Twist_Fermion));

  // Begin with pure gaussian random numbers
  FORALLSITES(i, s) {
    clear_TF(&(src[i]));
    for (j = 0; j < DIMF; j++) {                // Site fermions
#ifdef SITERAND
      src[i].Fsite.c[j].real = gaussian_rand_no(&(s->site_prn));
      src[i].Fsite.c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
      src[i].Fsite.c[j].real = gaussian_rand_no(&node_prn);
      src[i].Fsite.c[j].imag = gaussian_rand_no(&node_prn);
#endif
      for (mu = 0; mu < NUMLINK; mu++) {        // Link fermions
#ifdef SITERAND
        src[i].Flink[mu].c[j].real = gaussian_rand_no(&(s->site_prn));
        src[i].Flink[mu].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
        src[i].Flink[mu].c[j].real = gaussian_rand_no(&node_prn);
        src[i].Flink[mu].c[j].imag = gaussian_rand_no(&node_prn);
#endif
      }
      for (mu = 0; mu < NPLAQ; mu++) {         // Plaquette fermions
#ifdef SITERAND
        src[i].Fplaq[mu].c[j].real = gaussian_rand_no(&(s->site_prn));
        src[i].Fplaq[mu].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
        src[i].Fplaq[mu].c[j].real = gaussian_rand_no(&node_prn);
        src[i].Fplaq[mu].c[j].imag = gaussian_rand_no(&node_prn);
#endif
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

//  fermion_op(src, psim[0], PLUS);
//  fermion_op(psim[0], psim[1], MINUS);
//  printf("\n\n TEST 1\n");
//  FORALLSITES(i, s) {
//    source_norm = (double)magsq_TF(&(psim[1][i]));
//    if (source_norm * source_norm > 0) {
//      printf("%d %d %d %d %.4g\n",
//             s->x, s->y, s->z, s->t, source_norm);
//      for (mu = 0; mu < NUMLINK; mu++)
//        dumpmat(&(s->link[mu]));
//    }
//  }

//  fermion_op(src, psim[0],PLUS);
//  printf("\n\n TEST 2\n");
//  FORALLSITES(i, s) {
//    source_norm = (double)magsq_TF(&(psim[0][i]));
//    if (source_norm * source_norm > 0) {
//      printf("%d %d %d %d %.4g\n",
//             s->x, s->y, s->z, s->t, source_norm);
//      for (mu = 0; mu < NUMLINK; mu++)
//        dumpmat(&(s->link[mu]));
//    }
//  }
#endif

  // We now compute (Mdag M)^{1 / 8}.src
  for (i = 0; i < Norder; i++)
    shift[i] = shift8[i];

  avs_iters = congrad_multi_field(src, psim, niter, rsqmin, &size_r);
#ifdef DEBUG_CHECK
  node0_printf("Iters for source %d\n", avs_iters);
#endif
  // Reconstruct (Mdag M)^{1 / 8}.src from multi-mass CG solution psim
  FORALLSITES(i, s) {
    scalar_mult_TF(&(src[i]), ampdeg8, &(src[i]));
    for (j = 0; j < Norder; j++)
      scalar_mult_add_TF(&(src[i]), &(psim[j][i]), amp8[j], &(src[i]));
  }

  for (i = 0; i < Norder; i++)
    free(psim[i]);
  free(psim);
  return avs_iters;
}
// -----------------------------------------------------------------
