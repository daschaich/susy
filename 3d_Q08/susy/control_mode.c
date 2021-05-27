// -----------------------------------------------------------------
// Main procedure for N=4 SYM stochastic mode number calculations
// Adapted from adjoint SU(N) code by Georg Bergner

// This Giusti--Luescher method (arXiv:0812.3638) needs coefficients for the
// minmax polynomial approximation to the step function
// These are computed offline and copied into mode_coeffs.c
// The serial (GSL-based) code is in the mode_polynomial directory
//   https://www.gnu.org/software/gsl/
// It finds the smallest step_order that satisfies the input epsilon and delta
// It is based on code kindly shared by Agostino Patella
//   http://inspirehep.net/author/profile/A.Patella.1
#define CONTROL
#include "susy_includes.h"

#ifndef MODE
#error "Don't use control_mode unless compiling with -DMODE!"
#endif
#ifdef CHEB
#error "-DCHEB and -DMODE will clobber each other!"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt, dir, k;
  double ss_plaq, st_plaq, dtime, plpMod = 0.0;
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
  node0_printf("%.8g ", ss_plaq);
  node0_printf("%.8g\n", plpMod);
  node0_printf("BACTION %.8g\n", ss_plaq);

  // Load coefficients for minmax step function of given order
  // "star" is actually (Omega / Omega_*)^2, so take the square root
  coefficients();
  starSq = star;
  star = sqrt(starSq);
  node0_printf("Step function order %d epsilon %.4g delta %.4g star %.4g\n",
               step_order, step_eps, delta, star);

  // Run Giusti--Luescher stochastic mode number computation
  compute_mode();

  // Print results
  for (k = 0; k < numOmega; k++)
    node0_printf("MODE %.4g %.8g %.4g\n", Omega[k], mode[k], err[k]);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);

  normal_exit(0);         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
