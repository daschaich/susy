// -----------------------------------------------------------------
// Main procedure for N=4 SYM evolution and measurements
#define CONTROL
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt, dir;
  int traj_done, s_iters, avs_iters = 0, avm_iters = 0, Nmeas = 0;
  Real f_eps, g_eps;
  double dssplaq, dstplaq, dtime, plpMod = 0.0;
  double linktr[NUMLINK], linktr_ave, linktr_width;
  double link_det[NUMLINK], det_ave, det_width;
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
  setup_PtoP();
  setup_FQ();
  setup_rhmc();

#if defined(WLOOP) || defined(SMEAR)
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

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Check: compute initial plaquette and bosonic action
  d_plaquette(&dssplaq, &dstplaq);
  node0_printf("START %.8g %.8g %.8g ", dssplaq, dstplaq, dssplaq + dstplaq);
  dssplaq = d_gauge_action(NODET);
  node0_printf("%.8g\n", dssplaq / (double)volume);
  // N>2 determinant of traceless part of U.Udag crashes with unit config
  if (startflag != FRESH) {
    linktr_ave = d_link(linktr, &linktr_width, link_det, &det_ave, &det_width);
    node0_printf("FLINK");
    for (dir = XUP; dir < NUMLINK; dir++)
      node0_printf(" %.6g", linktr[dir]);
    node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);
    node0_printf("FLINK_DET");
    for (dir = XUP; dir < NUMLINK; dir++)
      node0_printf(" %.6g", link_det[dir]);
    node0_printf(" %.6g %.6g\n", det_ave, det_width);
  }

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
    // Tr[Udag.U] / N
    linktr_ave = d_link(linktr, &linktr_width, link_det, &det_ave, &det_width);
    node0_printf("FLINK");
    for (dir = XUP; dir < NUMLINK; dir++)
      node0_printf(" %.6g", linktr[dir]);
    node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);

    // Polyakov loop and plaquette measurements
    // Format: GMES Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq
    plp = ploop(&plpMod);
    d_plaquette(&dssplaq, &dstplaq);
    node0_printf("GMES %.8g %.8g %d %.8g %.8g ",
                 plp.real, plp.imag, s_iters, dssplaq, dstplaq);

    // Bosonic action (printed twice by request)
    // Might as well spit out volume average of Polyakov loop modulus
    dssplaq = d_gauge_action(NODET);
    node0_printf("%.8g ", dssplaq / (double)volume);
    node0_printf("%.8g\n", plpMod);
    node0_printf("BACTION %.8g\n", dssplaq / (double)volume);

    // Less frequent measurements every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
#ifdef SMEAR
#define MIN_PLAQ
      // Optionally smear before less frequent measurements
      node0_printf("Doing %d det-divided APE smearing steps ", Nsmear);
      node0_printf("with alpha=%.4g\n", alpha);

      // Check minimum plaquette in addition to averages
      node0_printf("BEFORE ");
      d_plaquette_lcl(&dssplaq, &dstplaq);    // Prints out MIN_PLAQ
      node0_printf(" %.8g %.8g\n", dssplaq, dstplaq);

      // Overwrite s->linkf
      // Save unsmeared links in Uinv (mom and f_U both already used)
      for (dir = XUP; dir < NUMLINK; dir++) {
        FORALLSITES(i, s)
          su3mat_copy_f(&(s->linkf[dir]), &(Uinv[dir][i]));
      }
      APE_smear(Nsmear, alpha, YESDET);
      node0_printf("AFTER  ");
      d_plaquette_lcl(&dssplaq, &dstplaq);    // Prints out MIN_PLAQ
      node0_printf(" %.8g %.8g\n", dssplaq, dstplaq);
#endif

      // Plaquette determinant
      measure_det();

      // Monitor widths of plaquette and plaquette determinant distributions
      widths();

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

      print_var3("PLCORR");
#endif

#ifdef CORR
      // Konishi and SUGRA correlators
      setup_P();
      d_konishi();
      d_correlator_r();
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
      rsymm(NODET);

      // Measure density of monopole world lines in non-diagonal cubes
      monopole();
#endif

#ifdef WLOOP
      // Gauge fix to easily access arbitrary displacements
      // Save un-fixed links to be written to disk if requested
      if (fixflag == COULOMB_GAUGE_FIX) {
        d_plaquette(&dssplaq, &dstplaq);    // To be printed below
        FORALLSITES(i, s) {
          for (dir = XUP; dir < NUMLINK; dir++)
            su3mat_copy_f(&(s->linkf[dir]), &(s->mom[dir]));
        }

        node0_printf("Fixing to Coulomb gauge...\n");
        double gtime = -dclock();

        // Gauge fixing arguments explained in generic/gaugefix.c
        // With first argument outside XUP, ..., TUP,
        // first four links are included in gauge-fixing condition
        gaugefix(TUP, 1.5, 5000, GAUGE_FIX_TOL, -1, -1);
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
      hvy_pot(NODET);
      hvy_pot(YESDET);

      // Save and restore links overwritten by polar projection
      // Don't use mom[TUP], which is already storing the un-fixed links
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->linkf[TUP]), &(s->f_U[TUP]));
      hvy_pot_polar();
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->f_U[TUP]), &(s->linkf[TUP]));

      // Restore the un-fixed links to be written to disk if requested
      if (fixflag == COULOMB_GAUGE_FIX) {
        FORALLSITES(i, s) {
          for (dir = XUP; dir < NUMLINK; dir++)
            su3mat_copy_f(&(s->mom[dir]), &(s->linkf[dir]));
        }
      }
#endif

#ifdef SMEAR
      // Restore unsmeared links from Uinv
      for (dir = XUP; dir < NUMLINK; dir++) {
        FORALLSITES(i, s)
          su3mat_copy_f(&(Uinv[dir][i]), &(s->linkf[dir]));
      }
#endif
    }
    fflush(stdout);
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final plaquette and bosonic action
  d_plaquette(&dssplaq, &dstplaq);
  node0_printf("STOP %.8g %.8g %.8g ",
               dssplaq, dstplaq, dssplaq + dstplaq);
  dssplaq = d_gauge_action(NODET);
  node0_printf("%.8g\n", dssplaq / (double)volume);

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
    save_lattice(saveflag, savefile);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
