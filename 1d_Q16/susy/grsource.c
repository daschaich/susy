// -----------------------------------------------------------------
// Routines for filling fields (momenta and pseudofermions)
// with gaussian random numbers
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Construct gaussian random momenta
// Either anti-hermitian matrices or conjugate to static diagonal phases
// Impose SITERAND rather than leaving this as a compile-time option
void ranmom() {
  register int i, j;
  register site *s;
  anti_hermitmat tah;

#ifndef SITERAND
  node0_printf("Error: Need to #define SITERAND\n");
  exit(1);
#endif

#ifndef UNGAUGED
#ifndef STATIC_GAUGE
  FORALLSITES(i, s)
    random_anti_hermitian(&(s->mom), &(s->site_prn));
#else
  // Careful -- must generate only one random array for whole lattice
  if (this_node == 0) {
    Real ave = 0.0;
    for (j = 0; j < NCOL; j++) {
      theta_mom[j] = gaussian_rand_no(&(lattice[0].site_prn));
      node0_printf("theta_mom[%d] = %.4g\n", j, theta_mom[j]);
      make_periodic(&(theta_mom[j]));
      ave += theta_mom[j];
    }

    // Determinant condition
    ave *= one_ov_N;
    for (j = 0; j < NCOL; j++)
      theta_mom[j] -= ave;
  }
  // Broadcast thetas from node0 to all other nodes
  broadcast_bytes((char *)&theta_mom, NCOL * sizeof(Real));
#endif
#endif

  // Scalar momenta
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      random_anti_hermitian(&tah, &(s->site_prn));
      uncompress_anti_hermitian(&tah, &(s->mom_X[j]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Construct a gaussian random vector R, return src = (Mdag M)^{1 / 8} R
// Need to invert despite the positive power, since it is fractional
// Return the number of iterations from the inversion
#ifndef PUREGAUGE
int grsource(matrix *src[NFERMION]) {
  register int i, j, k;
  register site *s;
  int avs_iters;
  Real size_r;
  complex grn;
  matrix ***psim = malloc(sizeof(matrix**) * Norder);

  // Allocate psim (will be zeroed in congrad_multi)
  for (j = 0; j < Norder; j++) {
    psim[j] = malloc(sizeof(matrix*) * NFERMION);
    for (i = 0; i < NFERMION; i++)
      psim[j][i] = malloc(sizeof(matrix) * sites_on_node);
  }

  // Begin with pure gaussian random numbers
  FORALLSITES(i, s) {
    for (k = 0; k < NFERMION; k++) {
      clear_mat(&(src[k][i]));
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[k][i]));
      }
    }
  }

#ifdef DEBUG_CHECK
  double source_norm = 0.0;
  FORALLSITES(i, s) {
    for (k = 0; k < NFERMION; k++) {
//      if (i != 0)
//        clear_mat(&(src[k][i]));
//      if (i == 0)
//        dumpmat(&src[k][i]);
      source_norm += (double)realtrace_nn(&(src[k][i]), &(src[k][i]));
    }
  }
  g_doublesum(&source_norm);
  node0_printf("source_norm in grsource %.4g\n", source_norm);
#endif

  // We now compute (Mdag M)^{1 / 8}.src
  for (i = 0; i < Norder; i++)
    shift[i] = shift8[i];

  avs_iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
#ifdef DEBUG_CHECK
  node0_printf("Iters for source %d\n", avs_iters);
#endif
  // Reconstruct (Mdag M)^{1 / 8}.src from multi-mass CG solution psim
  FORALLSITES(i, s) {
    for (k = 0; k < NFERMION; k++) {
      scalar_mult_matrix(&(src[k][i]), ampdeg8, &(src[k][i]));
      for (j = 0; j < Norder; j++)
        scalar_mult_sum_matrix(&(psim[j][k][i]), amp8[j], &(src[k][i]));
    }
  }

  for (i = 0; i < Norder; i++) {
    for (k = 0; k < NFERMION; k++)
      free(psim[i][k]);
    free(psim[i]);
  }
  free(psim);
  return avs_iters;
}
#endif
// -----------------------------------------------------------------
