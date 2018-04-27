// -----------------------------------------------------------------
// Main procedure for BFSS/BMN measurements on saved configurations
#define CONTROL
#include "susy_includes.h"

int main(int argc, char *argv[]) {
  int prompt, dir, j;
  double b_act, dtime, Xtr[NSCALAR], Xtr_ave, Xtr_width;
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

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Compute bosonic action, scalar squares and Polyakov loop to check config
  b_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]));
  node0_printf("START %.8g\n", b_act / (double)nt);

  Xtr_ave = scalar_trace(Xtr, &Xtr_width);
  node0_printf("SCALAR SQUARES");
  for (j = 0; j < NSCALAR; j++)
    node0_printf(" %.6g", Xtr[j]);
  node0_printf(" %.6g %.6g\n", Xtr_ave, Xtr_width);

  // Format: GMES Re(Polyakov) Im(Poyakov) cg_iters
  plp = ploop();
  node0_printf("GMES %.8g %.8g 0 ", plp.real, plp.imag);

  // More details of bosonic action (to match expected GMES output format)
  node0_printf("%.8g %.8g %.8g %.8g\n",
               b_act / (double)nt, Xtr[0] / (double)nt,
               Xtr[1] / (double)nt, Xtr[2] / (double)nt);

  // Main measurements
  // Monitor scalar eigenvalues
  // Format: SCALAR_EIG # ave width min max
  scalar_eig(ave_eigs, eig_widths, min_eigs, max_eigs);
  for (j = 0; j < NCOL; j++) {
    node0_printf("SCALAR_EIG %d %.6g %.6g %.6g %.6g\n",
                 j, ave_eigs[j], eig_widths[j], min_eigs[j], max_eigs[j]);
  }

#ifdef CORR
  // Konishi and SUGRA
  konishi();
#endif

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
