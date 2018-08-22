// -----------------------------------------------------------------
// Measure Ward identity involving eta.psi_a fermion bilinear
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// g_rand is random gaussian only for N components of site fermion
// and all components of link fermions
// src = Mdag g_rand
// Adapted from grsource to mimic pbp calculation
void bilinear_src(matrix *g_rand[NFERMION], matrix *src[NFERMION]) {
  register int i, j, k, mu;
  register site *s;
  complex grn;
  
  FORALLSITES(i, s) {
    for (k = 0; k < NFERMION; k++) {
      clear_TF(&(g_rand[k][i]));
      // Source all the fermions
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(g_rand[k][i]));
      }
    }
  }
  
  // Set up src = Mdag g_rand
  fermion_op(g_rand, src, MINUS);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure Ward identity involving eta.psi_a fermion bilinear,
//   Q sum_a tr(eta U_a Ubar_a) = tr(sum_b Dbar_b U_b sum_a U_a Ubar_a)
//                                 - sum_a tr(eta psi_a Ubar_a)
// Return total number of iterations
// Assume compute_DmuUmu() has already been run
int bilinearWard() {
  if (nsrc < 1)   // Only run if there are inversions to do
    return 0;
  register int i;
  register site *s;
  int mu, isrc, iters, tot_iters = 0, sav = Norder;
  Real size_r, norm;
  double sum = 0.0;
  double_complex tc, StoL, LtoS, ave = cmplx(0.0, 0.0);
  matrix tmat;
  Twist_Fermion *g_rand, *src, **psim;

  g_rand = malloc(sizeof *g_rand * sites_on_node);
  src = malloc(sizeof *src * sites_on_node);

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(Twist_Fermion*));
  psim[0] = malloc(sizeof(Twist_Fermion) * sites_on_node);
  shift[0] = 0;

  // Normalization: sum over NUMLINK but divide by volume
  // and divide by 2kappa as discussed on 15--17 December 2013
  norm = 2.0 * kappa * (Real)volume;

  for (isrc = 0; isrc < nsrc; isrc++) {
    // Make random source g_rand, now including trace piece
    // Hit it with Mdag to get src, invert to get M^{-1} g_rand
    // congrad_multi initializes psim
    bilinear_src(g_rand, src, DIMF);
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    tot_iters += iters;
#ifdef DEBUG_CHECK
    dump_TF(&(psim[0][10]));    // Check what components are non-zero
#endif

    // Now construct bilinear sum_a psi_a Udag_a eta
    // All fields are on the same site, no gathers
    // Negative sign from generator normalization
    StoL = cmplx(0.0, 0.0);
    LtoS = cmplx(0.0, 0.0);
    FORALLSITES(i, s) {
      FORALLDIR(mu) {
        mult_na(&(psim[0][i].Flink[mu]), &(s->link[mu]), &tmat);
        tc = complextrace_an(&(g_rand[i].Fsite), &tmat);
        CSUM(StoL, tc);
#ifdef DEBUG_CHECK
        printf("StoL[%d](%d) %d (%.4g, %.4g)\n",
               i, mu, isrc, tc.real, tc.imag);
#endif

        mult_an(&(s->link[mu]), &(psim[0][i].Fsite), &tmat);
        tc = complextrace_an(&(g_rand[i].Flink[mu]), &tmat);
        CSUM(LtoS, tc);
#ifdef DEBUG_CHECK
        printf("LtoS[%d](%d) %d (%.4g, %.4g)\n",
               i, mu, isrc, tc.real, tc.imag);
#endif
      }
    }
    g_dcomplexsum(&StoL);
    g_dcomplexsum(&LtoS);
    CDIVREAL(StoL, norm, StoL);
    CDIVREAL(LtoS, norm, LtoS);
    node0_printf("susy %.6g %.6g %.6g %.6g ( %d of %d ) %d\n",
                 StoL.real, StoL.imag, LtoS.real, LtoS.imag,
                 isrc + 1, nsrc, iters);

    CSUM(ave, StoL);
    CDIF(ave, LtoS);
  }
  // Normalize by number of sources (already averaged over volume)
  CDIVREAL(ave, 2.0 * (Real)nsrc, ave);

  // Now add gauge piece, including plaquette determinant term
  // Accumulate sum_a U_a Udag_a in tmat
  // Multiply by DmuUmu and trace
  FORALLSITES(i, s) {
    mult_na(&(s->link[0]), &(s->link[0]), &tmat);   // Initialize
    for (mu = 1; mu < NUMLINK; mu++)
      mult_na_sum(&(s->link[mu]), &(s->link[mu]), &tmat);
    tc = complextrace_nn(&(DmuUmu[i]), &tmat);
    // Make sure trace really is real
    if (fabs(tc.imag) > IMAG_TOL) {
      printf("node%d WARNING: Im(sum[%d]) = %.4g > %.4g\n",
             this_node, i, tc.imag, IMAG_TOL);
    }
    sum += tc.real;
  }
  g_doublesum(&sum);

  // Average gauge term over volume, print along with difference
  sum /= (Real)volume;
  node0_printf("SUSY %.6g %.6g %.6g %.6g ( ave over %d )\n",
               ave.real, ave.imag, sum, sum - ave.real, nsrc);

  // Reset multi-mass CG and clean up
  Norder = sav;
  free(g_rand);
  free(src);
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------
