// -----------------------------------------------------------------
// Measure Ward identity involving eta.psi_a fermion bilinear
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Separate random gaussian fermion sources
// for the two sets of NCHIRAL_FERMION indices j and k
// Then invert on each, src_j = Mdag g_rand_j
// Adapted from grsource to mimic QCD pbp calculation
void bilinear_src(matrix *g_rand[NFERMION], matrix *src[NFERMION]) {

  register int i, j, k, n;
  register site *s;
  complex grn;

  FORALLSITES(i, s) {
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      k = j + NCHIRAL_FERMION;
      // Clear all g_rand_* components
      clear_mat(&(g_rand[j][i]));
      clear_mat(&(g_rand[k][i]));

      // Source first (second) half of g_rand_j (g_rand_k)
      for (n = 0; n < DIMF; n++) {
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[n]), &grn, &(g_rand[j][i]));

#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[n]), &grn, &(g_rand[k][i]));
      }
    }
  }

  // Set up src = Mdag g_rand
  fermion_op(g_rand, src, MINUS);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure a few fermion bilinears (with k = j + NCHIRAL_FERMION):
//   psidag[j] {psi[k](t) - BC_+ * U(t) psi[k](t+1) Udag(t)}        [Dplus]
//   psidag[k] {psi[j](t) - BC_- * Udag(t-1) psi[j](t-1) U(t-1)}   [Dminus]
//   TODO: THE REST...
//   TODO: CHECK ORDER OF MATRIX MULTS...
// Use tempmat for temporary storage
// Return total number of iterations
int bilinear() {
  if (nsrc < 1)   // Only run if there are inversions to do
    return 0;
  register int i;
  register site *s;
  int j, k, m, n, isrc, iters, tot_iters = 0, sav = Norder;
  Real size_r, norm = (Real)nt;
  double_complex tc, Yukawa;
  double_complex ave_Yukawa = cmplx(0.0, 0.0);
  matrix *g_rand[NFERMION], *src[NFERMION], tmat;
  matrix ***psim = malloc(sizeof(matrix**));
  msg_tag *tag[NCHIRAL_FERMION];

  // Hack a basic CG out of the multi-mass CG
  // congrad_multi will initializes psim
  Norder = 1;
  shift[0] = 0;
  psim[0] = malloc(sizeof(matrix*) * NFERMION);
  for (i = 0; i < NFERMION; i++) {
    psim[0][i] = malloc(sizeof(matrix) * sites_on_node);
    g_rand[i] = malloc(sizeof(matrix) * sites_on_node);
    src[i] = malloc(sizeof(matrix) * sites_on_node);
  }

  build_Gamma_X();

  for (isrc = 0; isrc < nsrc; isrc++) {
    // Make random sources g_rand_*, hit with Mdag to get src_*
    bilinear_src(g_rand, src);

    // Invert to get psim[k] = M_{kj}^{-1} g_rand_j
    // Contract psim[k] with g_rand_j to isolate Dplus bilinear
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    tot_iters += iters;

    Yukawa = cmplx(0.0, 0.0);
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      m = j + NCHIRAL_FERMION;
      for (k = 0; k < NCHIRAL_FERMION; k++) {
        n = k + NCHIRAL_FERMION;
        FORALLSITES(i, s) {
          mult_na(&(Gamma_X[k][j][i]), &(psim[0][n][i]), &tmat);
          tc = complextrace_an(&(g_rand[j][i]), &tmat);
          CDIF(Yukawa, tc);

          mult_nn(&(psim[0][n][i]), &(Gamma_X[j][k][i]), &tmat);
          tc = complextrace_an(&(g_rand[j][i]), &tmat);
          CDIF(Yukawa, tc);

          mult_na(&(Gamma_X[k][j][i]), &(psim[0][k][i]), &tmat);
          tc = complextrace_an(&(g_rand[m][i]), &tmat);
          CDIF(Yukawa, tc);

          mult_nn(&(psim[0][k][i]), &(Gamma_X[j][k][i]), &tmat);
          tc = complextrace_an(&(g_rand[m][i]), &tmat);
          CDIF(Yukawa, tc);
        }
      }

      // Last 2 gammas are diagonal
      k = NSCALAR - 2;
      n = NSCALAR - 1;
      FORALLSITES(i, s) {
        mult_na(&(s->X[k]), &(psim[0][m][i]), &tmat);
        mult_nn_sum(&(psim[0][m][i]), &(s->X[k]), &tmat);
        tc = complextrace_an(&(g_rand[j][i]), &tmat);
        CDIF(Yukawa, tc);

        mult_na(&(s->X[k]), &(psim[0][j][i]), &tmat);
        mult_nn_sum(&(psim[0][j][i]), &(s->X[k]), &tmat);
        tc = complextrace_an(&(g_rand[m][i]), &tmat);
        CDIF(Yukawa, tc);

        mult_na(&(s->X[n]), &(psim[0][m][i]), &tmat);
        mult_nn_sum(&(psim[0][m][i]), &(s->X[n]), &tmat);
        tc = complextrace_an(&(g_rand[j][i]), &tmat);
        CDIF(Yukawa, tc);

        mult_na(&(s->X[n]), &(psim[0][j][i]), &tmat);
        mult_nn_sum(&(psim[0][j][i]), &(s->X[n]), &tmat);
        tc = complextrace_an(&(g_rand[m][i]), &tmat);
        CDIF(Yukawa, tc);

      }
    }
    CMULREAL(Yukawa, 0.5, Yukawa);

    g_dcomplexsum(&Yukawa);

    // Average over Nt and print estimate from each stochastic source
    CDIVREAL(Yukawa, norm, Yukawa);
    node0_printf("bilin %.6g %.6g ( %d of %d ) %d\n",
                 Yukawa.real, Yukawa.imag,
                 isrc + 1, nsrc, iters);

    CSUM(ave_Yukawa, Yukawa);
  }

  // Normalize by number of sources (already averaged over Nt)
  CDIVREAL(ave_Yukawa, (Real)nsrc, ave_Yukawa);

  node0_printf("BILIN %.6g %.6g ( ave over %d )\n",
               ave_Yukawa.real, ave_Yukawa.imag, nsrc);

  // Reset multi-mass CG and clean up
  Norder = sav;
  for (i = 0; i < NFERMION; i++) {
    free(psim[0][i]);
    free(g_rand[i]);
    free(src[i]);
  }
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------
