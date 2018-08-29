// -----------------------------------------------------------------
// Measure Ward identity involving eta.psi_a fermion bilinear
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// g_rand is random gaussian fermion source, then src = Mdag g_rand
// Adapted from grsource to mimic QCD pbp calculation
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
// Measure a few fermion bilinears (with k = j + NCHIRAL_FERMION):
//   psi[k](t) - BC_+ * U(t) psi[k](t+1) Udag(t)        [Dplus]
//   psi[j](t) - BC_- * Udag(t-1) psi[j](t-1) U(t-1)    [Dminus]
//   TODO: THE REST...
// Return total number of iterations
int bilinearWard() {
  if (nsrc < 1)   // Only run if there are inversions to do
    return 0;
  register int i;
  register site *s;
  int mu, isrc, iters, tot_iters = 0, sav = Norder;
  Real size_r, norm;
  double sum = 0.0;
  double_complex tc, StoL, LtoS, ave = cmplx(0.0, 0.0);
  matrix tmat, *g_rand[NFERMION], *src[NFERMION];
  matrix ***psim = malloc(sizeof(matrix**));
  msg_tag *tag[NCHIRAL_FERMION];

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  shift[0] = 0;
  psim[0] = malloc(sizeof(matrix*) * NFERMION);
  for (i = 0; i < NFERMION; i++) {
    psim[0][i] = malloc(sizeof(matrix) * sites_on_node);
    g_rand[i] = malloc(sizeof(matrix) * sites_on_node);
    src[i] = malloc(sizeof(matrix) * sites_on_node);
  }

  for (isrc = 0; isrc < nsrc; isrc++) {
    // Make random source g_rand, now including trace piece
    // Hit it with Mdag to get src, invert to get M^{-1} g_rand
    // congrad_multi initializes psim
    bilinear_src(g_rand, src);
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    tot_iters += iters;

    // TODO: EDITS NEEDED FROM HERE...
    // Now construct bilinears
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
    node0_printf("bilin %.6g %.6g %.6g %.6g ( %d of %d ) %d\n",
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
  node0_printf("BILIN %.6g %.6g %.6g %.6g ( ave over %d )\n",
               ave.real, ave.imag, sum, sum - ave.real, nsrc);
  // ...TO HERE EDITS NEEDED TODO

  // Reset multi-mass CG and clean up
  Norder = sav;
  for (i = 0; i < NFERMION; i++) {
    free(g_rand[i]);
    free(src[i]);
    free(psim[0][i]);
  }
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------
