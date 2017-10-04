// -----------------------------------------------------------------
// Main procedure for N=4 SYM evolution and measurements
#define CONTROL
#include "susy_includes.h"

int main(int argc, char *argv[]) {
  int prompt, dir, j;
  int traj_done, s_iters, avs_iters = 0, avm_iters = 0, Nmeas = 0;
  Real f_eps, g_eps;
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
  setup_PtoP();
  setup_FQ();
  setup_rhmc();

#ifdef SMEAR
  register int i;
  register site *s;
  double max_plaq = 0.0;
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
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("START %.8g %.8g %.8g ", ss_plaq, st_plaq, ss_plaq + st_plaq);
  ss_plaq = gauge_action(NODET);
  node0_printf("%.8g\n", ss_plaq / (double)volume);
  // N>2 determinant of traceless part of U.Udag crashes with unit config
  if (startflag != FRESH) {
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
    linktr_ave = link_trace(linktr, &linktr_width,
                            link_det, &det_ave, &det_width);
    node0_printf("FLINK");
    FORALLDIR(dir)
      node0_printf(" %.6g", linktr[dir]);
    node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);

    // Polyakov loop and plaquette measurements
    // Format: GMES Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq
    plp = ploop(TUP, NODET, &plpMod);
    plaquette(&ss_plaq, &st_plaq);
    node0_printf("GMES %.8g %.8g %d %.8g %.8g ",
                 plp.real, plp.imag, s_iters, ss_plaq, st_plaq);

    // Bosonic action (printed twice by request)
    // Might as well spit out volume average of Polyakov loop modulus
    ss_plaq = gauge_action(NODET);
    node0_printf("%.8g ", ss_plaq / (double)volume);
    node0_printf("%.8g\n", plpMod);
    node0_printf("BACTION %.8g\n", ss_plaq / (double)volume);

    // Full and polar-projected Wilson lines in all five basis dirs
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

    // Less frequent measurements every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
#ifdef SMEAR
      // Optionally smear before less frequent measurements
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
      max_plaq = local_plaquette(&ss_plaq, &st_plaq);   // Prints out MIN_PLAQ
      node0_printf(" %.8g %.8g %.8g\n", ss_plaq, st_plaq, max_plaq);

      // Overwrite s->link
      // Save unsmeared links in UpsiU (mom and f_U both already used)
      FORALLSITES(i, s) {
        FORALLDIR(dir)
          mat_copy(&(s->link[dir]), &(UpsiU[dir][i]));
      }
      if (smearflag == STOUT_SMEAR)
        stout_smear(Nsmear, alpha);
      else if (smearflag == APE_SMEAR)
        APE_smear(Nsmear, alpha, YESDET);
      node0_printf("AFTER  ");
      max_plaq = local_plaquette(&ss_plaq, &st_plaq);   // Prints out MIN_PLAQ
      node0_printf(" %.8g %.8g %.8g\n", ss_plaq, st_plaq, max_plaq);

      // Update plaqdet, Uinv, DmuUmu and Fmunu with smeared links
      compute_plaqdet();
      compute_Uinv();
      compute_DmuUmu();
      compute_Fmunu();
#endif

#ifdef PL_CORR
      // Polyakov loop correlator
      // Compute ploop_corr and take the absolute square of its FFT
      ploop_c();      // Compute Polyakov loop at each spatial site
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
      konishi();
      correlator_r();
#endif

#ifdef BILIN
      // Ward identity involving eta.psi_a fermion bilinear
      Nmeas++;
      avm_iters += bilinearWard();
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
      // Save and restore original links in mom before fixing gauge
      if (fixflag == COULOMB_GAUGE_FIX) {
#ifndef SMEAR
        register int i;
        register site *s;
#endif
        plaquette(&ss_plaq, &st_plaq);    // To be printed below
        FORALLSITES(i, s) {
          FORALLDIR(dir)
            mat_copy(&(s->link[dir]), &(s->mom[dir]));
        }

        node0_printf("Fixing to Coulomb gauge...\n");
        double gtime = -dclock();

        // Gauge fixing arguments explained in generic/gaugefix.c
        // With first argument outside XUP, ..., TUP,
        // first four links are included in gauge-fixing condition
        gaugefix(TUP, 1.5, 5000, GAUGE_FIX_TOL, -1, -1);
        gtime += dclock();
        node0_printf("GFIX time = %.4g seconds\n", gtime);
        node0_printf("BEFORE %.8g %.8g\n", ss_plaq, st_plaq);
        plaquette(&ss_plaq, &st_plaq);
        node0_printf("AFTER  %.8g %.8g\n", ss_plaq, st_plaq);
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
        mat_copy(&(s->link[TUP]), &(s->f_U[TUP]));
      hvy_pot_polar();
      FORALLSITES(i, s)
        mat_copy(&(s->f_U[TUP]), &(s->link[TUP]));

      // Restore the un-fixed links to be written to disk if requested
      if (fixflag == COULOMB_GAUGE_FIX) {
        FORALLSITES(i, s) {
          FORALLDIR(dir)
            mat_copy(&(s->mom[dir]), &(s->link[dir]));
        }
      }
#endif

#ifdef SMEAR
      // Restore unsmeared links from UpsiU
      FORALLSITES(i, s) {
        FORALLDIR(dir)
          mat_copy(&(UpsiU[dir][i]), &(s->link[dir]));
      }

      // Recompute unsmeared plaqdet, Uinv, DmuUmu and Fmunu
      compute_plaqdet();
      compute_Uinv();
      compute_DmuUmu();
      compute_Fmunu();
#endif
    }
    fflush(stdout);
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final plaquette and bosonic action
  // Reset DmuUmu and Fmunu in case they were smeared above
  compute_DmuUmu();
  compute_Fmunu();
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("STOP %.8g %.8g %.8g ",
               ss_plaq, st_plaq, ss_plaq + st_plaq);
  ss_plaq = gauge_action(NODET);
  node0_printf("%.8g\n", ss_plaq / (double)volume);

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
