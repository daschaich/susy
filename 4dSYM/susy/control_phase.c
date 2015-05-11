// -----------------------------------------------------------------
// Main procedure for N=4 SYM pfaffian phase
#define CONTROL
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt;
  double dssplaq, dstplaq, dtime;
  complex plp = cmplx(99, 99);
#ifndef PHASE
  node0_printf("Don't use control_phase unless compiling with -DPHASE!\n");
  terminate(1);
#endif

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

  // Load input and run (loop removed)
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

  // Do "local" measurements to check evolution
  // Polyakov loop measurement
  plp = ploop();

  // Tr[Udag.U] / N and plaquette measurements
  d_link(0);
  d_plaquette(&dssplaq, &dstplaq);

  // Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq
  node0_printf("GMES %.8g %.8g 0 %.8g %.8g ",
               plp.real, plp.imag, dssplaq, dstplaq);

  // Bosonic action (printed twice by request)
  dssplaq = d_gauge_action(NODET);
  node0_printf("%.8g\n", dssplaq / (double)volume);
  node0_printf("BACTION %.8g\n", dssplaq / (double)volume);

#if 0
  // Optionally fix to Coulomb gauge to check gauge invariance
  d_plaquette(&dssplaq, &dstplaq);    // To be printed below
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
#endif

  // Main measurement: pfaffian phase
  d_phase();

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------
