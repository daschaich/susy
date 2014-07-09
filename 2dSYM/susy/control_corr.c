// -----------------------------------------------------------------
// Main procedure for N=4 SYM correlator and bilinear measurements
#define CONTROL
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt;
  double dplaq, dtime;
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

#ifdef PL_CORR
  // Set up Fourier transform for Polyakov loop correlator
  register int i;
  register site *s;
  int key[4] = {1, 1, 1, 0};
  int restrict[4];
  Real space_vol = (Real)nx;
  setup_restrict_fourier(key, restrict);
#endif

  // Load input and run (loop removed)
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Check: compute initial plaquette and bosonic action
  d_plaquette(&dplaq);
  node0_printf("START %.8g ", dplaq);
  dplaq = d_gauge_action();
  node0_printf("%.8g\n", dplaq / (double)volume);

  // Do "local" measurements to check evolution
  // Polyakov loop measurement
  plp = ploop();

  // Tr[Udag.U] / N and plaquette measurements
  d_link();
  d_plaquette(&dplaq);

  // Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq
  node0_printf("GMES %.8g %.8g 0 %.8g ",
               plp.real, plp.imag, dplaq);

  // Bosonic action (printed twice by request)
  dplaq = d_gauge_action();
  node0_printf("%.8g\n", dplaq / (double)volume);
  node0_printf("BACTION %.8g\n", dplaq / (double)volume);

  // Main measurements
#ifdef DET
  // Plaquette determinant
  measure_det();

  // Flavor transformation observables -- uses find_det
  d_flavor();
#endif

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
#endif

#ifdef BILIN
  // Fermion bilinear and supersymmetry transformation vevs
  int avm_iters = d_bilinear();
  avm_iters += d_susyTrans();
#endif

#ifdef WLOOP
  // Gauge-fixed Wilson loops
  if (fixflag == COULOMB_GAUGE_FIX) {
    d_plaquette(&dplaq);    // To be printed below
    node0_printf("Fixing to Coulomb gauge...\n");
    double gtime = -dclock();

    // Gauge fixing arguments explained in generic/gaugefix.c
    gaugefix(TUP, 1.5, 500, GAUGE_FIX_TOL, -1, -1);
    gtime += dclock();
    node0_printf("GFIX time = %.4g seconds\n", gtime);
    node0_printf("BEFORE %.8g\n", dplaq);
    d_plaquette(&dplaq);
    node0_printf("AFTER  %.8g\n", dplaq);
  }
  else if (fixflag == NO_GAUGE_FIX) { // Braces suppress compiler warning
    node0_printf("Gauge fixing skipped\n");
  }
  else {
    node0_printf("ERROR: only COULOMB_GAUGE_FIX ");
    node0_printf("and NO_GAUGE_FIX supported\n");
    terminate(1);
  }

  // We only consider the fundamental links
  // The irrep links would require separate routines
  // with the correct data structures (su3_matrix instead of su3_matrix_f)
  hvy_pot();

  // Polar projection overwrites temporal links
  // Need to save and restore if checking effect of gauge transformation
  FORALLSITES(i, s)
    su3mat_copy_f(&(s->linkf[TUP]), &(s->mom[TUP]));

  hvy_pot_polar();

  // Restore non-unitary temporal links
  FORALLSITES(i, s)
    su3mat_copy_f(&(s->mom[TUP]), &(s->linkf[TUP]));
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
