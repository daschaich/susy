// -----------------------------------------------------------------
// Measure Ward identity involving eta.psi_a fermion bilinear
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Separate random gaussian fermion sources
// for the two sets of NCHIRAL_FERMION indices j and k
// Then invert on each, src_j = Mdag g_rand_j
// Adapted from grsource to mimic QCD pbp calculation
void bilinear_src(matrix *g_rand_j[NFERMION], matrix *src_j[NFERMION],
                  matrix *g_rand_k[NFERMION], matrix *src_k[NFERMION]) {

  register int i, j, k, n;
  register site *s;
  complex grn;

  FORALLSITES(i, s) {
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      k = j + NCHIRAL_FERMION;
      // Clear all g_rand_* components
      clear_mat(&(g_rand_j[j][i]));
      clear_mat(&(g_rand_j[k][i]));
      clear_mat(&(g_rand_k[j][i]));
      clear_mat(&(g_rand_k[k][i]));

      // Source first (second) half of g_rand_j (g_rand_k)
      for (n = 0; n < DIMF; n++) {
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[n]), &grn, &(g_rand_j[j][i]));

#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[n]), &grn, &(g_rand_k[k][i]));
      }
    }
  }

  // Set up src = Mdag g_rand
  fermion_op(g_rand_j, src_j, MINUS);
  fermion_op(g_rand_k, src_k, MINUS);
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
  int j, k, isrc, iters, tot_iters = 0, sav = Norder;
  Real size_r, norm = (Real)nt;
  double_complex tc, bilinear_Dplus, bilinear_Dminus;
  double_complex ave_Dplus = cmplx(0.0, 0.0), ave_Dminus = cmplx(0.0, 0.0);
  matrix *g_rand_j[NFERMION], *src_j[NFERMION];
  matrix *g_rand_k[NFERMION], *src_k[NFERMION], tmat, tmat2;
  matrix ***psim = malloc(sizeof(matrix**));
  msg_tag *tag[NCHIRAL_FERMION];

  // Hack a basic CG out of the multi-mass CG
  // congrad_multi will initializes psim
  Norder = 1;
  shift[0] = 0;
  psim[0] = malloc(sizeof(matrix*) * NFERMION);
  for (i = 0; i < NFERMION; i++) {
    psim[0][i] = malloc(sizeof(matrix) * sites_on_node);
    g_rand_j[i] = malloc(sizeof(matrix) * sites_on_node);
    g_rand_k[i] = malloc(sizeof(matrix) * sites_on_node);
    src_j[i] = malloc(sizeof(matrix) * sites_on_node);
    src_k[i] = malloc(sizeof(matrix) * sites_on_node);
  }

  for (isrc = 0; isrc < nsrc; isrc++) {
    // Make random sources g_rand_*, hit with Mdag to get src_*
    bilinear_src(g_rand_j, g_rand_k, src_j, src_k);

    // Invert to get psim[k] = M_{kj}^{-1} g_rand_j
    // Contract psim[k] with g_rand_j to isolate Dplus bilinear
    iters = congrad_multi(src_j, psim, niter, rsqmin, &size_r);
    tot_iters += iters;

    bilinear_Dplus = cmplx(0.0, 0.0);
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      k = j + NCHIRAL_FERMION;
      tag[j] = start_gather_field(psim[0][k], sizeof(matrix),
                                  TUP, EVENANDODD, gen_pt[j]);
    }

    for (j = 0; j < NCHIRAL_FERMION; j++) {
      k = j + NCHIRAL_FERMION;
      wait_gather(tag[j]);
      FORALLSITES(i, s) {   // psim[k] in gen_pt[j]
        mult_nn(&(s->link), (matrix *)(gen_pt[j][i]), &tmat);
        scalar_mult_na(&tmat, &(s->link), s->bc[0], &tmat2);
        dif_matrix(&(psim[0][k][i]), &tmat2);
        tc = complextrace_an(&(g_rand_j[j][i]), &tmat2);
        CSUM(bilinear_Dplus, tc);
      }
      cleanup_gather(tag[j]);
    }
    g_dcomplexsum(&bilinear_Dplus);

    // Invert to get psim[j] = M_{jk}^{-1} g_rand_k
    // Contract psim[j] with g_rand_k to isolate Dminus bilinear
    iters = congrad_multi(src_k, psim, niter, rsqmin, &size_r);
    tot_iters += iters;

    bilinear_Dminus = cmplx(0.0, 0.0);
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      k = j + NCHIRAL_FERMION;
      FORALLSITES(i, s) {
        mult_an(&(s->link), &(psim[0][j][i]), &tmat);
        scalar_mult_nn(&tmat, &(s->link), -1.0, &(tempmat[i]));
      } // Include negative sign in product that is then gathered
      tag[0] = start_gather_field(tempmat, sizeof(matrix),
                                  TDOWN, EVENANDODD, gen_pt[0]);

      wait_gather(tag[0]);
      FORALLSITES(i, s) {   // Subtraction via negative sign above
        scalar_mult_add_matrix(&(psim[0][j][i]), (matrix *)(gen_pt[0][i]),
                               s->bc[1], &tmat);
        tc = complextrace_an(&(g_rand_k[k][i]), &tmat);
        CSUM(bilinear_Dminus, tc);
      }
      cleanup_gather(tag[0]);
    }
    g_dcomplexsum(&bilinear_Dminus);

    // Average over Nt and print estimate from each stochastic source
    CDIVREAL(bilinear_Dplus, norm, bilinear_Dplus);
    CDIVREAL(bilinear_Dminus, norm, bilinear_Dminus);
    node0_printf("bilin %.6g %.6g %.6g %.6g ( %d of %d ) %d\n",
                 bilinear_Dplus.real, bilinear_Dplus.imag,
                 bilinear_Dminus.real, bilinear_Dminus.imag,
                 isrc + 1, nsrc, iters);

    CSUM(ave_Dplus, bilinear_Dplus);
    CSUM(ave_Dminus, bilinear_Dminus);
  }

  // Normalize by number of sources (already averaged over Nt)
  CDIVREAL(ave_Dplus, (Real)nsrc, ave_Dplus);
  CDIVREAL(ave_Dminus, (Real)nsrc, ave_Dminus);

  node0_printf("BILIN %.6g %.6g %.6g %.6g ( ave over %d )\n",
               ave_Dplus.real, ave_Dplus.imag,
               ave_Dminus.real, ave_Dminus.imag, nsrc);

  // Reset multi-mass CG and clean up
  Norder = sav;
  for (i = 0; i < NFERMION; i++) {
    free(psim[0][i]);
    free(g_rand_j[i]);
    free(g_rand_k[i]);
    free(src_j[i]);
    free(src_k[i]);
  }
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------
