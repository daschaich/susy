// -----------------------------------------------------------------
// Measure Ward identity involving eta.psi_a fermion bilinear
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// g_rand is random gaussian only for N components of site fermion
// and all components of link fermions
// src = Mdag g_rand
// Adapted from grsource to mimic pbp calculation
void bilinear_src(Twist_Fermion *g_rand, Twist_Fermion *src, int N) {
  register int i, j, mu;
  register site *s;
  complex grn;

#ifdef DEBUG_CHECK  // Test that fermion_op connects proper components
  FORALLSITES(i, s) {
    clear_TF(&(g_rand[i]));       // Zero plaquette fermions
    clear_TF(&(src[i]));          // To be safe/explicit
  }
  int index = node_index(1, 1, 1, 1);
  g_rand[index].Flink[0].c[0].real = 1;
  fermion_op(g_rand, src, PLUS);
  FORALLSITES(i, s) {
    if (magsq_TF(&(src[i])) > IMAG_TOL) {
      printf("(%d %d %d %d)\n", s->x, s->y, s->z, s->t);
      dump_TF(&(src[i]));
    }
  }
  terminate(1);
#endif

  FORALLSITES(i, s) {
    clear_TF(&(g_rand[i]));
    // Source either all or traceless site fermions, depending on N
    // The last Lambda[DIMF - 1] (N = DIMF) is proportional to the identity
    // The others (N = DIMF - 1) are traceless
    for (j = 0; j < N; j++) {                   // Site fermions
#ifdef SITERAND
      grn.real = gaussian_rand_no(&(s->site_prn));
      grn.imag = gaussian_rand_no(&(s->site_prn));
#else
      grn.real = gaussian_rand_no(&node_prn);
      grn.imag = gaussian_rand_no(&node_prn);
#endif
      c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(g_rand[i].Fsite));
    }

    // Source all link and plaquette fermions
    FORALLDIR(mu) {                           // Link fermions
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(g_rand[i].Flink[mu]));
      }
    }
    for (mu = 0; mu < NPLAQ; mu++) {         // Plaquette fermions
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(g_rand[i].Fplaq[mu]));
      }
    }
  }

  // Set up src = Mdag g_rand
  fermion_op(g_rand, src, MINUS);
  FORALLSITES(i, s)
    scalar_mult_sum_TF(&(g_rand[i]), fmass, &(src[i]));

#ifdef DEBUG_CHECK
//  dump_TF(&(src[10]));

  // Have checked that each component of g_rand is unit-normalized
  double source_norm;
  FORALLSITES(i, s)
    source_norm += (double)magsq_TF(&(g_rand[i]));
  g_doublesum(&source_norm);
  source_norm /= (volume * N);
  node0_printf("average bilinear_src component magSq %.4g\n", source_norm);
#endif
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

  g_rand = malloc(sites_on_node * sizeof(*g_rand));
  src = malloc(sites_on_node * sizeof(*src));

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(Twist_Fermion));
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
