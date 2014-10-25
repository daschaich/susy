// -----------------------------------------------------------------
// Measure fermion bilinear and Konishi superpartner correlator
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// g_rand is random gaussian only for N components of site fermion
// and all components of link fermions
// src = Mdag g_rand
// Adapted from grsource to mimic pbp calculation
void bilinsrc(Twist_Fermion *g_rand, Twist_Fermion *src, int N) {
  register int i, j, mu;
  register site *s;

  FORALLSITES(i, s) {
    clear_TF(&(g_rand[i]));       // Zero plaquette fermions
    clear_TF(&(src[i]));          // To be safe/explicit

    // Source either all or traceless site fermions, depending on N
    // The last Lambda[DIMF - 1] (N = DIMF) is proportional to the identity
    // The others (N = DIMF - 1) are traceless
    for (j = 0; j < N; j++) {                 // Site fermions
#ifdef SITERAND
      g_rand[i].Fsite.c[j].real = gaussian_rand_no(&(s->site_prn));
      g_rand[i].Fsite.c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
      g_rand[i].Fsite.c[j].real = gaussian_rand_no(&node_prn);
      g_rand[i].Fsite.c[j].imag = gaussian_rand_no(&node_prn);
#endif
    }

    // Source all link fermions
    for (j = 0; j < DIMF; j++) {
      for (mu = 0; mu < NUMLINK; mu++) {        // Link fermions
#ifdef SITERAND
        g_rand[i].Flink[mu].c[j].real = gaussian_rand_no(&(s->site_prn));
        g_rand[i].Flink[mu].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
        g_rand[i].Flink[mu].c[j].real = gaussian_rand_no(&node_prn);
        g_rand[i].Flink[mu].c[j].imag = gaussian_rand_no(&node_prn);
#endif
      }
    }
  }

  // Set up src = Mdag g_rand
  fermion_op(g_rand, src, MINUS);
  FORALLSITES(i, s)
    scalar_mult_add_TF(&(src[i]), &(g_rand[i]), fmass, &(src[i]));

#ifdef DEBUG_CHECK
//  dump_TF(&(src[10]));

  // Have checked that each component of g_rand is unit-normalized
  double source_norm;
  FORALLSITES(i, s)
    source_norm += (double)magsq_TF(&(g_rand[i]));
  g_doublesum(&source_norm);
  source_norm /= (volume * N);
  node0_printf("average bilinsrc component magSq %.4g\n", source_norm);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure the fermion bilinear
//   sum_a tr(eta Ubar_a psi_a) - tr(eta)/N tr(Ubar_a psi_a)
// Return total number of iterations
int d_bilinear() {
  if (nsrc < 1)   // Only run if there are inversions to do
    return 0;
  register int i;
  register site *s;
  int mu, isrc, iters, tot_iters = 0;
  Real size_r, norm;
  double_complex tc, StoL, LtoS, ave = cmplx(0.0, 0.0);
  su3_vector tvec;
  Twist_Fermion *g_rand, *src, **psim;

  fermion_rep();    // Make sure site->link are up to date
  g_rand = malloc(sites_on_node * sizeof(*g_rand));
  src = malloc(sites_on_node * sizeof(*src));

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(Twist_Fermion));
  shift[0] = 0;             // Reset in update or grsource

  // Normalization: sum over NUMLINK but divide by volume
  // and divide by 2kappa as discussed on 15--17 December 2013
  norm = 2.0 * kappa * (Real)volume;

  for (isrc = 0; isrc < nsrc; isrc++) {
    // Make random source g_rand without trace piece
    // Hit it with Mdag to get src, invert to get M^{-1} g_rand
    // (The psim are zeroed in congrad_multi_field)
    bilinsrc(g_rand, src, DIMF - 1);
    iters = congrad_multi_field(src, psim, niter, rsqmin, &size_r);
    tot_iters += iters;
//    dump_TF(&(psim[0][10]));    // Check what components are non-zero

    // Now construct bilinear -- all on same site, no gathers
    StoL = cmplx(0.0, 0.0);
    LtoS = cmplx(0.0, 0.0);
    FORALLSITES(i, s) {
      // First site-to-link, g^dag . sum_a psi_a Vdag_a
      for (mu = 0; mu < NUMLINK; mu++) {
        mult_su3_vec_adj_mat(&(psim[0][i].Flink[mu]), &(s->link[mu]), &tvec);
        tc = su3_dot(&(g_rand[i].Fsite), &tvec);
        CSUM(StoL, tc);
      }

      // Now link-to-site, sum_a g_a^dag . eta Vdag_a
      // Omit last component of site propagator to avoid trace piece
      psim[0][i].Fsite.c[DIMF - 1] = cmplx(0.0, 0.0);
      for (mu = 0; mu < NUMLINK; mu++) {
        mult_adj_su3_mat_vec(&(s->link[mu]), &(psim[0][i].Fsite), &tvec);
        tc = su3_dot(&(g_rand[i].Flink[mu]), &tvec);
        CSUM(LtoS, tc);
      }
    }
    g_dcomplexsum(&StoL);
    g_dcomplexsum(&LtoS);

    // Normalize and print
    CDIVREAL(StoL, norm, StoL);   // Sum over NUMLINK
    CDIVREAL(LtoS, norm, LtoS);
    node0_printf("bilin %.6g %.6g %.6g %.6g ( %d of %d ) %d\n",
                 StoL.real, StoL.imag, LtoS.real, LtoS.imag,
                 isrc + 1, nsrc, iters);

    CSUM(ave, StoL);
    CDIF(ave, LtoS);
  }
  // Normalize by number of sources (already averaged over volume)
  CDIVREAL(ave, 2.0 * (Real)nsrc, ave);
  node0_printf("BILIN %.6g %.6g ( ave over %d )\n", ave.real, ave.imag, nsrc);

  // Reset multi-mass CG and clean up
  Norder = DEGREE;
  free(g_rand);
  free(src);
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure the supersymmetry transformation
//   Q sum_a tr(eta Ubar_a U_a) = sum_{a, b} tr(DbUb Ubar_a U_a)
//                              - sum_a tr(eta Ubar_a psi_a)
// Return total number of iterations
int d_susyTrans() {
  if (nsrc < 1)   // Only run if there are inversions to do
    return 0;
  register int i;
  register site *s;
  int mu, isrc, iters, tot_iters = 0;
  Real size_r, norm;
  double sum = 0.0;
  double_complex tc, StoL, LtoS, ave = cmplx(0.0, 0.0);
  su3_vector tvec;
  su3_matrix_f tmat1, tmat2;
  Twist_Fermion *g_rand, *src, **psim;

  fermion_rep();    // Make sure site->link are up to date
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
    // (The psim are zeroed in congrad_multi_field)
    bilinsrc(g_rand, src, DIMF);
    iters = congrad_multi_field(src, psim, niter, rsqmin, &size_r);
    tot_iters += iters;
//    dump_TF(&(psim[0][10]));    // Check what components are non-zero

    // Now construct bilinear -- all on same site, no gathers
    StoL = cmplx(0.0, 0.0);
    LtoS = cmplx(0.0, 0.0);
    FORALLSITES(i, s) {
      // First site-to-link, g^dag . sum_a psi_a Vdag_a
      for (mu = 0; mu < NUMLINK; mu++) {
        mult_su3_vec_adj_mat(&(psim[0][i].Flink[mu]), &(s->link[mu]), &tvec);
        tc = su3_dot(&(g_rand[i].Fsite), &tvec);
        CSUM(StoL, tc);
      }

      // Now link-to-site, sum_a g_a^dag . eta Vdag_a
      for (mu = 0; mu < NUMLINK; mu++) {
        mult_adj_su3_mat_vec(&(s->link[mu]), &(psim[0][i].Fsite), &tvec);
        tc = su3_dot(&(g_rand[i].Flink[mu]), &tvec);
        CSUM(LtoS, tc);
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

  // Now add gauge piece
  // Accumulate sum_a U_a Udag_a in tmat1
  // Multiply by DmuUmu into tmat2 and trace
  compute_DmuUmu();     // Compute sum_b [U_b Udag_b - Udag_b U_b]
  FORALLSITES(i, s) {
    mult_su3_na_f(&(s->linkf[0]), &(s->linkf[0]), &tmat1);
    for (mu = 1; mu < NUMLINK; mu++) {
      mult_su3_na_f(&(s->linkf[mu]), &(s->linkf[mu]), &tmat2);
      add_su3_matrix_f(&tmat1, &tmat2, &tmat1);
    }
    mult_su3_nn_f(&(s->DmuUmu), &tmat1, &tmat2);
    tc = trace_su3_f(&tmat2);
    // Make sure trace really is real
    if (abs(tc.imag) > IMAG_TOL) {
      printf("node%d WARNING: sum[%d]) = %.4g > %.4g)\n",
             this_node, i, tc.imag, IMAG_TOL);
    }
    sum += tc.real;
  }
  g_doublesum(&sum);

  // Average second term over volume, print along with difference
  sum /= (Real)volume;
  node0_printf("SUSY %.6g %.6g %.6g %.6g ( ave over %d )\n",
               ave.real, ave.imag, sum, sum - ave.real, nsrc);

  // Reset multi-mass CG and clean up
  Norder = DEGREE;
  free(g_rand);
  free(src);
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------
