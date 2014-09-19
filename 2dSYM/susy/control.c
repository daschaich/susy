// -----------------------------------------------------------------
// Main procedure for N=4 SYM evolution and measurements
#define CONTROL
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int traj_done, prompt, s_iters, avs_iters = 0, avm_iters = 0, Nmeas = 0;
  Real f_eps, g_eps;
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
  setup_rhmc();

#ifdef WLOOP
  register int i, mu;
  register site *s;
#endif

#ifdef PL_CORR
  // Set up Fourier transform for Polyakov loop correlator
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
  d_link();

  // Perform warmup trajectories
  f_eps = traj_length / (Real)nsteps[0];
  g_eps = f_eps / (Real)(2 * nsteps[1]);
  node0_printf("f_eps %.4g g_eps %.4g\n", f_eps, g_eps);
  for (traj_done = 0; traj_done < warms; traj_done++)
    update();
  node0_printf("WARMUPS COMPLETED\n");

  // Perform trajectories with measurements
  // But WITHOUT reunitarizing!
  for (traj_done = 0; traj_done < trajecs; traj_done++) {
    s_iters = update();
    avs_iters += s_iters;

    // Do "local" measurements every trajectory!
    // Polyakov loop measurement
    plp = ploop();

    // Tr[Udag.U] / N and plaquette measurements
    d_link();
    d_plaquette(&dplaq);
//    d_plaquette_frep(&dplaq_frep);

    // Re(Polyakov) Im(Poyakov) cg_iters plaq
    node0_printf("GMES %.8g %.8g %d %.8g ",
                 plp.real, plp.imag, s_iters, dplaq);

    // Bosonic action (printed twice by request)
    dplaq = d_gauge_action();
    node0_printf("%.8g\n", dplaq / (double)volume);
    node0_printf("BACTION %.8g\n", dplaq / (double)volume);

    // Less frequent measurements every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
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
#endif

#ifdef BILIN
      // Ward identity violations
      Nmeas++;
      avm_iters += d_susyTrans();

      // Don't run d_bilinear() for now
      // In the future it may be useful
      // to compare U(1) vs. SU(N) contributions
//      avm_iters += d_bilinear();
#endif

#ifdef CORR
  // R symmetry transformations -- use find_det and adjugate
      rsymm();

  // Measure density of monopole world lines in non-diagonal cubes
//      monopole();
#endif

#ifdef WLOOP
      // Gauge-fixed Wilson loops
      // Save un-fixed links to be saved if requested
      if (fixflag == COULOMB_GAUGE_FIX) {
        d_plaquette(&dplaq);    // To be printed below
        FORALLSITES(i, s) {
          for (mu = XUP; mu < NUMLINK; mu++)
            su3mat_copy_f(&(s->linkf[mu]), &(s->mom[mu]));
        }

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
      hvy_pot();

      // Save and restore links overwritten by polar projection
      // Don't use mom[TUP], which is already storing the un-fixed links
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->linkf[TUP]), &(s->f_U[TUP]));
      hvy_pot_polar();
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->f_U[TUP]), &(s->linkf[TUP]));

      // Restore the un-fixed links to be saved if requested
      if (fixflag == COULOMB_GAUGE_FIX) {
        FORALLSITES(i, s) {
          for (mu = XUP; mu < NUMLINK; mu++)
            su3mat_copy_f(&(s->mom[mu]), &(s->linkf[mu]));
        }
      }
#endif
    }
    fflush(stdout);
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final plaquette and bosonic action
  d_plaquette(&dplaq);
  node0_printf("STOP %.8g ", dplaq);
  dplaq = d_gauge_action();
  node0_printf("%.8g\n", dplaq / (double)volume);

  node0_printf("Average CG iters for steps: %.4g\n",
               (double)avs_iters / trajecs);
  if (Nmeas > 0) {
    node0_printf("Average CG iters for measurements: %.4g\n",
                 (double)avm_iters / Nmeas);
  }
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  // total_iters is accumulated in the multiCG itself
  // Should equal total for steps plus measurements
  node0_printf("total_iters = %d\n\n", total_iters);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);
  return 0;
}
// -----------------------------------------------------------------
