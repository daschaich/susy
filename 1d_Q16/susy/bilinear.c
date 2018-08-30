// -----------------------------------------------------------------
// Measure Ward identity involving eta.psi_a fermion bilinear
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// g_rand is random gaussian fermion source, then src = Mdag g_rand
// Adapted from grsource to mimic QCD pbp calculation
void bilinear_src(matrix *g_rand[NFERMION], matrix *src[NFERMION]) {
  register int i, j, k;
  register site *s;
  complex grn;

  FORALLSITES(i, s) {
    for (k = 0; k < NFERMION; k++) {
      clear_mat(&(g_rand[k][i]));
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
int bilinear() {
  if (nsrc < 1)   // Only run if there are inversions to do
    return 0;
  register int i;
  register site *s;
  int j, k, isrc, iters, tot_iters = 0, sav = Norder;
  Real size_r, norm = (Real)nt;
  double_complex tc, bilinear_Dplus, bilinear_Dminus;
  double_complex ave_Dplus = cmplx(0.0, 0.0), ave_Dminus = cmplx(0.0, 0.0);
  matrix tmat, tmat2, *g_rand[NFERMION], *src[NFERMION];
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

    // Now construct bilinears
    // Negative sign from generator normalization
    bilinear_Dplus = cmplx(0.0, 0.0);
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      k = j + NCHIRAL_FERMION;
      tag[j] = start_gather_field(psim[0][k], sizeof(matrix),
                                  TUP, EVENANDODD, gen_pt[j]);
    }
    
    for (j = 0; j < NCHIRAL_FERMION; j++) {     // Overwrite dest
      k = j + NCHIRAL_FERMION;
      wait_gather(tag[j]);
      FORALLSITES(i, s) {
        mult_nn(&(s->link), (matrix *)(gen_pt[j][i]), &tmat);
        scalar_mult_na(&tmat, &(s->link), s->bc[0], &(tmat2));
        dif_matrix(&(psim[0][k][i]), &(tmat2));
        tc = complextrace_an(&(g_rand[j][i]), &tmat2);
        CSUM(bilinear_Dplus, tc);
      }
      cleanup_gather(tag[j]);
    }
    g_dcomplexsum(&bilinear_Dplus);
    
    
    bilinear_Dminus = cmplx(0.0, 0.0);
    // Include negative sign in product that is then gathered
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      k = j + NCHIRAL_FERMION;
      FORALLSITES(i, s) {
        mult_an(&(s->link), &(psim[0][j][i]), &tmat);
        scalar_mult_nn(&tmat, &(s->link), -1.0, &(tempmat[i]));
      }
      tag[0] = start_gather_field(tempmat, sizeof(matrix),
                               TDOWN, EVENANDODD, gen_pt[0]);
      
      wait_gather(tag[0]);
      FORALLSITES(i, s) {           // Overwrite dest
        scalar_mult_add_matrix(&(psim[0][j][i]), (matrix *)(gen_pt[0][i]),
                               s->bc[1], &tmat);
        tc = complextrace_an(&(g_rand[k][i]), &tmat);
        CSUM(bilinear_Dminus, tc);
      }
      cleanup_gather(tag[0]);
    }
    g_dcomplexsum(&bilinear_Dminus);
    
    CDIVREAL(bilinear_Dplus, norm, bilinear_Dplus);
    CDIVREAL(bilinear_Dminus, norm, bilinear_Dminus);
    node0_printf("bilin %.6g %.6g %.6g %.6g ( %d of %d ) %d\n",
                 bilinear_Dplus.real, bilinear_Dplus.imag,
                 bilinear_Dminus.real, bilinear_Dminus.imag,
                 isrc + 1, nsrc, iters);

    CSUM(ave_Dplus, bilinear_Dplus);
    CDIF(ave_Dminus, bilinear_Dminus);
  }
  // Normalize by number of sources (already averaged over volume)
  CDIVREAL(ave_Dplus, (Real)nsrc, ave_Dplus);
  CDIVREAL(ave_Dminus, (Real)nsrc, ave_Dminus);

  node0_printf("BILIN %.6g %.6g %.6g %.6g ( ave over %d )\n",
               ave_Dplus.real, ave_Dplus.imag, ave_Dminus.real, ave_Dminus.imag, nsrc);

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
