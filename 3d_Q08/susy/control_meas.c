// -----------------------------------------------------------------
// Main procedure for N=4 SYM measurements on saved configurations,
// including scalar correlators, Wilson loops, Ward identity violations
// and discrete R symmetry observables
#define CONTROL
#include "susy_includes.h"

int main(int argc, char *argv[]) {
  int prompt, dir, j;
  double ss_plaq, st_plaq, dtime, plpMod = 0.0;
  double linktr[NUMLINK], linktr_ave, linktr_width;
  double link_det[NUMLINK], det_ave, det_width;
  double ave_eigs[NCOL], eig_widths[NCOL], min_eigs[NCOL], max_eigs[NCOL];
  complex plp = cmplx(99.0, 99.0);

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();
  setup_lambda();
  epsilon();
  setup_PtoV();
  setup_FQ();

#ifdef SMEAR
  double max_plaq = 0.0;
#endif

#ifdef WLOOP
  register int i;
  register site *s;
#endif

#ifdef PL_CORR
  // Set up Fourier transform for Polyakov loop correlator
  int key[3] = {1, 1, 0};
  int restrict[4];
  Real space_vol = (Real)(nx * ny);
  setup_restrict_fourier(key, restrict);
#endif

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Uncomment this block to print plaquette and determinant distributions
  // Be sure to uncomment PLAQ_DIST and DET_DIST, and run in serial
//  local_plaquette(&ss_plaq, &st_plaq);
//  node0_printf("\n");
//  if (G < IMAG_TOL)
//    G = 999;
//  compute_DmuUmu();
//  dtime += dclock();
//  node0_printf("\nTime = %.4g seconds\n", dtime);
//  fflush(stdout);
//  return 0;

  // Check: compute initial plaquette and bosonic action
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("START %.8g %.8g %.8g ", ss_plaq, st_plaq, ss_plaq + st_plaq);
  ss_plaq = gauge_action(NODET);
  node0_printf("%.8g\n", ss_plaq / (double)volume);

  // Do "local" measurements to check configuration
  // Tr[Udag.U] / N
  linktr_ave = link_trace(linktr, &linktr_width,
                          link_det, &det_ave, &det_width);
  node0_printf("FLINK");
  FORALLDIR(dir)
    node0_printf(" %.6g", linktr[dir]);
  node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);
  node0_printf("FLINK_DET");
  FORALLDIR(dir)
    node0_printf(" %.6g", link_det[dir]);
  node0_printf(" %.6g %.6g\n", det_ave, det_width);

  // Polyakov loop and plaquette measurements
  // Format: GMES Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq
  plp = ploop(TUP, NODET, &plpMod);
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("GMES %.8g %.8g 0 %.8g %.8g ",
               plp.real, plp.imag, ss_plaq, st_plaq);

  // Bosonic action (printed twice by request)
  // Might as well spit out volume average of Polyakov loop modulus
  ss_plaq = gauge_action(NODET) / (double)volume;
//  ss_plaq = gauge_action(YESDET) / (double)volume;
  node0_printf("%.8g ", ss_plaq);
  node0_printf("%.8g\n", plpMod);
  node0_printf("BACTION %.8g\n", ss_plaq);

  // Full and polar-projected Wilson lines in all basis dirs
  node0_printf("LINES      ");
  FORALLDIR(dir) {
    plp = ploop(dir, NODET, &plpMod);
    node0_printf(" %.6g %.6g", plp.real, plp.imag);
  }
  node0_printf("\nLINES_POLAR");
  FORALLDIR(dir) {
    plp = ploop(dir, YESDET, &plpMod);
    node0_printf(" %.6g %.6g", plp.real, plp.imag);
  }
  node0_printf("\n");

  // Test polar-projected Polyakov loop eigenvalues
//  plp = ploop_eig(TUP, YESDET, &plpMod);

#ifdef SMEAR
  // Optionally smear before main measurements
  // NO_SMEAR sets Nsmear = 0
  if (smearflag == STOUT_SMEAR) {
    node0_printf("Doing %d stout smearing steps ", Nsmear);
    node0_printf("with rho=%.4g\n", alpha);
  }
  else if (smearflag == APE_SMEAR) {
    node0_printf("Doing %d polar-projected APE smearing steps ", Nsmear);
    node0_printf("with alpha=%.4g\n", alpha);
  }

  // Check minimum plaquette in addition to averages
  node0_printf("BEFORE ");
  max_plaq = local_plaquette(&ss_plaq, &st_plaq);     // Prints out MIN_PLAQ
  node0_printf(" %.8g %.8g %.8g\n", ss_plaq, st_plaq, max_plaq);

  // Overwrite s->link
  if (smearflag == STOUT_SMEAR)
    stout_smear(Nsmear, alpha);
  else if (smearflag == APE_SMEAR)
    APE_smear(Nsmear, alpha, YESDET);
  node0_printf("AFTER  ");
  max_plaq = local_plaquette(&ss_plaq, &st_plaq);      // Prints out MIN_PLAQ
  node0_printf(" %.8g %.8g %.8g\n", ss_plaq, st_plaq, max_plaq);

  // Update plaquette determinants, DmuUmu and Fmunu with smeared links
  compute_plaqdet();
  compute_Uinv();
  compute_DmuUmu();
  compute_Fmunu();
#endif

  // Main measurements
  // Plaquette determinant
  measure_det();

  // Monitor widths of plaquette and plaquette determinant distributions
  widths();

  // Monitor scalar eigenvalues
  // Format: SCALAR_EIG # ave width min max
  scalar_eig(NODET, ave_eigs, eig_widths, min_eigs, max_eigs);
  for (j = 0; j < NCOL; j++) {
    node0_printf("UUBAR_EIG %d %.6g %.6g %.6g %.6g\n",
                 j, ave_eigs[j], eig_widths[j], min_eigs[j], max_eigs[j]);
  }
  scalar_eig(YESDET, ave_eigs, eig_widths, min_eigs, max_eigs);
  for (j = 0; j < NCOL; j++) {
    node0_printf("POLAR_EIG %d %.6g %.6g %.6g %.6g\n",
                 j, ave_eigs[j], eig_widths[j], min_eigs[j], max_eigs[j]);
  }

#ifdef PL_CORR
  // Polyakov loop correlator
  // Compute ploop_corr and take the absolute square of its FFT
  ploop_c();      // Computes Polyakov loop at each spatial site
  restrict_fourier(F_OFFSET(ploop_corr),
                   F_OFFSET(fft1), F_OFFSET(fft2),
                   sizeof(complex), FORWARDS);

  FORALLSITES(i, s)
    CMULJ_((s->ploop_corr), (s->ploop_corr), (s->print_var));

  // Invert FFT in place
  restrict_fourier(F_OFFSET(print_var),
                   F_OFFSET(fft1), F_OFFSET(fft2),
                   sizeof(complex), BACKWARDS);

  // Divide by volume for correct inverse FFT
  // Divide by another volume to average convolution
  FORALLSITES(i, s)
    CDIVREAL((s->print_var), space_vol * space_vol, (s->print_var));

  print_var2("PLCORR");
#endif

#ifdef CORR
  // Konishi and SUGRA correlators
  konishi();
  // correlator_r(); Not doing it for 3D, for now at least
#endif

#ifdef BILIN
  // Ward identity involving eta.psi_a fermion bilinear
  int avm_iters = bilinearWard();
#endif

#ifdef CORR
  // R symmetry transformations -- uses LAPACK invert
  rsymm();

  // Measure density of monopole world lines in non-diagonal cubes
  monopole();
#endif

#ifdef WLOOP
  // Gauge fix to Coulomb gauge
  // This lets us easily access arbitrary displacements
  // Gauge fixing arguments explained in generic/gaugefix.c
  // With first argument outside XUP, ..., TUP,
  // first four links are included in gauge-fixing condition
  if (fixflag == COULOMB_GAUGE_FIX) {
    plaquette(&ss_plaq, &st_plaq);    // To be printed below
    node0_printf("Fixing to Coulomb gauge...\n");
    double gtime = -dclock();
    gaugefix(TUP, 1.5, 5000, GAUGE_FIX_TOL, -1, -1);
    gtime += dclock();
    node0_printf("GFIX time = %.4g seconds\n", gtime);
    node0_printf("BEFORE %.8g %.8g\n", ss_plaq, st_plaq);
    plaquette(&ss_plaq, &st_plaq);
    node0_printf("AFTER  %.8g %.8g\n", ss_plaq, st_plaq);
  }
  hvy_pot(NODET);
  hvy_pot(YESDET);

  // Save and restore links overwritten by polar projection
  FORALLSITES(i, s)
    mat_copy(&(s->link[TUP]), &(s->mom[TUP]));
  hvy_pot_polar();
  FORALLSITES(i, s)
    mat_copy(&(s->mom[TUP]), &(s->link[TUP]));
#endif

  node0_printf("RUNNING COMPLETED\n");
#ifdef BILIN
  node0_printf("CG iters for measurements: %d\n", avm_iters);
  // total_iters is accumulated in the multiCG itself
  // Should equal total for bilinear measurements
  node0_printf("total_iters = %d\n", total_iters);
#endif
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  normal_exit(0);         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
