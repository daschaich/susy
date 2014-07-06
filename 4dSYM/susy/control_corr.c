// -----------------------------------------------------------------
// Main procedure for N=4 SYM correlator and bilinear measurements
#define CONTROL
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt, mu;
  double dssplaq, dstplaq, dtime;
  complex plp = cmplx(99, 99);

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
  setup_PtoP();
  setup_FQ();

#ifdef WLOOP
  register int i;
  register site *s;
#endif

#ifdef PL_CORR
  // Set up Fourier transform for Polyakov loop correlator
  int key[4] = {1, 1, 1, 0};
  int restrict[4];
  Real space_vol = (Real)(nx * ny * nz);
  setup_restrict_fourier(key, restrict);
#endif

  // Load input and run (loop removed)
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Check: compute initial plaquette and bosonic action
  d_plaquette(&dssplaq, &dstplaq);
  node0_printf("START %.8g %.8g %.8g ", dssplaq, dstplaq, dssplaq + dstplaq);
  dssplaq = d_gauge_action();
  node0_printf("%.8g\n", dssplaq / (double)volume);

  // Do "local" measurements to check evolution
  // Polyakov loop measurement
  plp = ploop();

  // Tr[Udag.U] / N and plaquette measurements
  d_link();
  d_plaquette(&dssplaq, &dstplaq);

  // Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq
  node0_printf("GMES %.8g %.8g 0 %.8g %.8g ",
               plp.real, plp.imag, dssplaq, dstplaq);

  // Bosonic action (printed twice by request)
  dssplaq = d_gauge_action();
  node0_printf("%.8g\n", dssplaq / (double)volume);
  node0_printf("BACTION %.8g\n", dssplaq / (double)volume);

#ifdef STOUT
#define MIN_PLAQ
  // Optionally smear before main measurements
  node0_printf("Doing %d stout smearing steps with rho=%.4g...\n",
               Nstout, rho);

  // Check minimum plaquette in addition to averages
  node0_printf("BEFORE ");
  d_plaquette_lcl(&dssplaq, &dstplaq);    // Prints out MIN_PLAQ
  node0_printf(" %.8g %.8g\n", dssplaq, dstplaq);

  // Overwrites s->linkf, saving original values in thin_link field
  block_stout(Nstout, rho);
  node0_printf("AFTER  ");
  d_plaquette_lcl(&dssplaq, &dstplaq);    // Prints out MIN_PLAQ
  node0_printf(" %.8g %.8g\n", dssplaq, dstplaq);
#endif

  // Main measurements
  // Plaquette determinant
  measure_det();

#ifdef PL_CORR
  // Polyakov loop correlator
  // Compute ploop_corr and take the absolute square of its FFT
  ploop_c();
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

  print_var3("PLCORR");
#endif

#ifdef CORR
  // Konishi and SUGRA correlators
  d_correlator();
  d_correlator_r();
#endif

#ifdef BILIN
  // Fermion bilinear and supersymmetry transformation vevs
  int avm_iters = d_bilinear();
  avm_iters += d_susyTrans();
#endif

  // R symmetry transformations -- use find_det and adjugate
  rsymm();

  // Measure density of monopole world lines in non-diagonal cubes
//  monopole();

#ifdef WLOOP
  // First calculate a few Wilson loops more directly, using explicit paths
  // Save and restore all links overwritten by polar projection
  hvy_pot_loop();
  FORALLSITES(i, s) {
    for (mu = XUP; mu < NUMLINK; mu++)
      su3mat_copy_f(&(s->linkf[mu]), &(s->mom[mu]));
  }
  hvy_pot_polar_loop();
  FORALLSITES(i, s) {
    for (mu = XUP; mu < NUMLINK; mu++)
      su3mat_copy_f(&(s->mom[mu]), &(s->linkf[mu]));
  }

  // Now gauge fix to easily access arbitrary displacements
  if (fixflag == COULOMB_GAUGE_FIX) {
    d_plaquette(&dssplaq, &dstplaq);    // To be printed below
    node0_printf("Fixing to Coulomb gauge...\n");
    double gtime = -dclock();

    // Gauge fixing arguments explained in generic/gaugefix.c
    // With first argument outside XUP, ..., TUP,
    // first four links are included in gauge-fixing condition
    gaugefix(NODIR, 1.5, 5000, GAUGE_FIX_TOL, -1, -1, 0, NULL, NULL);
    gtime += dclock();
    node0_printf("GFIX time = %.4g seconds\n", gtime);
    node0_printf("BEFORE %.8g %.8g\n", dssplaq, dstplaq);
    d_plaquette(&dssplaq, &dstplaq);
    node0_printf("AFTER  %.8g %.8g\n", dssplaq, dstplaq);
  }
  else if (fixflag == NO_GAUGE_FIX) { // Braces suppress compiler warning
    node0_printf("Gauge fixing skipped\n");
  }
  else {
    node0_printf("ERROR: only COULOMB_GAUGE_FIX ");
    node0_printf("and NO_GAUGE_FIX supported\n");
    terminate(1);
  }
  hvy_pot();

  // Save and restore links overwritten by polar projection
  FORALLSITES(i, s)
    su3mat_copy_f(&(s->linkf[DIR_5]), &(s->mom[DIR_5]));
  hvy_pot_polar();
  FORALLSITES(i, s)
    su3mat_copy_f(&(s->mom[DIR_5]), &(s->linkf[DIR_5]));
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
  return 0;
}
// -----------------------------------------------------------------
