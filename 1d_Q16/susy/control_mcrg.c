// -----------------------------------------------------------------
// Main procedure for BFSS/BMN RG blocking
// and measurements of blocked observables
#define CONTROL
#include "susy_includes.h"

int main(int argc, char *argv[]) {
  register int i;
  register site *s;
  int prompt, dir, j = nt, bl, blmax;
  double dtime, Xtr[NSCALAR], Xtr_ave, Xtr_width;
  double ave_eigs[NCOL], eig_widths[NCOL], min_eigs[NCOL], max_eigs[NCOL];
  complex plp = cmplx(99.0, 99.0);

#ifndef MCRG
  node0_printf("Don't use control_mcrg unless compiling with -DMCRG!\n");
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

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Works even if we can only block down to odd j > 4
  blmax = 0;
  while (j % 2 == 0 && j > 2) {    // While j is even
    j /= 2;
    blmax++;
  }

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

  // ---------------------------------------------------------------
  // First print unblocked observables
#ifdef CORR
  // Unblocked Konishi and SUGRA operators
  blocked_ops(0, 0);
#endif

  // Unblocked scalar eigenvalues
  // Format: SCALAR_EIG bl # ave width min max
  scalar_eig(ave_eigs, eig_widths, min_eigs, max_eigs);
  for (j = 0; j < NCOL; j++) {
    node0_printf("SCALAR_EIG 0 %d %.6g %.6g %.6g %.6g\n",
                 j, ave_eigs[j], eig_widths[j], min_eigs[j], max_eigs[j]);
  }
  // ---------------------------------------------------------------



  // ---------------------------------------------------------------
  // Loop over blocking levels (automatically determined)
  // Can probably merge in the unblocked checks above with a little work
  j = 1;
  for (bl = 1; bl <= blmax; bl++) {
    j *= 2;
    node0_printf("Blocking %d gives L = %d\n", bl, nx / j);
    block_mcrg(bl);

#ifdef CORR
    // Blocked Konishi and SUGRA correlators
    blocked_ops(0, bl);
#endif

    // Blocked Polyakov and Wilson loops
    blocked_ploop(0, bl);
    blocked_rsymm(0, bl);

    // Blocked scalar eigenvalues
    scalar_eig(ave_eigs, eig_widths, min_eigs, max_eigs);
    for (j = 0; j < NCOL; j++) {
      node0_printf("SCALAR_EIG %d ", bl);
      node0_printf("%d %.6g %.6g %.6g %.6g\n",
                   j, ave_eigs[j], eig_widths[j], min_eigs[j], max_eigs[j]);
    }
  }
  // ---------------------------------------------------------------

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
