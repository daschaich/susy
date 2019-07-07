// -----------------------------------------------------------------
// Main procedure for BFSS/BMN pfaffian phase
#define CONTROL
#include "susy_includes.h"

int main(int argc, char *argv[]) {
  int prompt, j;
  double b_act, dtime, Xtr[NSCALAR], Xtr_ave, Xtr_width;
  complex plp = cmplx(99.0, 99.0);
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
  setup_gamma();

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Compute bosonic action, scalar squares and Polyakov loop to check config
  b_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]), &(Xtr[3]));
  node0_printf("START %.8g\n", b_act / (double)nt);

  Xtr_ave = scalar_trace(Xtr, &Xtr_width);
  node0_printf("SCALAR SQUARES");
  for (j = 0; j < NSCALAR; j++)
    node0_printf(" %.6g", Xtr[j]);
  node0_printf(" %.6g %.6g\n", Xtr_ave, Xtr_width);

  // Print Polyakov loop eigenvalues followed by trace itself
  // Format: GMES Re(Polyakov) Im(Poyakov) cg_iters
  plp = ploop_eig();
  node0_printf("GMES %.8g %.8g 0 ", plp.real, plp.imag);

  // More details of bosonic action (to match expected GMES output format)
  node0_printf("%.8g %.8g %.8g %.8g\n",
               b_act / (double)nt, Xtr[0] / (double)nt,
               Xtr[1] / (double)nt, Xtr[2] / (double)nt);

  // Main measurement: pfaffian
  phase();

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  normal_exit(0);         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
