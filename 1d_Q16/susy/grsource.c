// -----------------------------------------------------------------
// Routines for filling fields (momenta and pseudofermions)
// with gaussian random numbers
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Construct gaussian random momentum matrices
// as sum of SU(N) generators with gaussian random coefficients
void ranmom() {
  register int i, j;
  register site *s;
  complex grn;

  FORALLSITES(i, s) {
      clear_mat(&(s->mom));
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&(s->node_prn));
        grn.imag = gaussian_rand_no(&(s->node_prn));
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(s->mom));
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
  matrix **psim[NFERMION];

  // Allocate psim (will be zeroed in congrad_multi)
  for (j = 0; j < NFERMION; j++) {
    psim[j] = malloc(sizeof(matrix*) * Norder);
  
  for (i = 0; i < Norder; i++)
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
      //    if (i != 0)
      //      clear_mat(&(src[k][i]));
      //    if (i == 0)
      //      dumpmat(&src[k][i]);
      source_norm += (double)realtrace_nn(&(src[k][i]), &(src[k][i]));
    }
  }
  g_doublesum(&source_norm);
  node0_printf("source_norm in grsource %.4g\n", source_norm);

//  DSq(src, dest);
//  printf("\n\n TEST 1\n");
//  FORALLSITES(i, s) {
//  for (k = 0; k < NFERMION; k++) {
//    source_norm = (double)realtrace_nn(&(dest[k][i]), &(dest[k][i]));
//    if (fabs(source_norm) > IMAG_TOL)
//      printf("%d %d %.4g\n", s->t, k, source_norm);
//    }
//  }

//  fermion_op(src, dest, PLUS);
//  printf("\n\n TEST 2\n");
//  FORALLSITES(i, s) {
//  for (k = 0; k < NFERMION; k++) {
//    source_norm = (double)realtrace_nn(&(dest[k][i]), &(dest[k][i]));
//    if (fabs(source_norm) > IMAG_TOL)
//      printf("%d %d %.4g\n", s->t, k, source_norm);
//    }
//  }
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
        scalar_mult_sum_matrix(&(psim[k][j][i]), amp8[j], &(src[k][i]));
    }
  }
  
  for (k = 0; k < NFERMION; k++) {
    for (i = 0; i < Norder; i++)
      free(psim[k][i]);
    free(psim[k]);
  }
  return avs_iters;
}
#endif
// -----------------------------------------------------------------
