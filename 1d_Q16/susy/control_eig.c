// -----------------------------------------------------------------
// Main procedure for BFSS/BMN eigenvalues
#define CONTROL
#include "susy_includes.h"

int main(int argc, char *argv[]) {
  int j, prompt;
  double b_act, dtime, Xtr[NSCALAR], Xtr_ave, Xtr_width;
  complex plp = cmplx(99.0, 99.0);
  int ivec, total_iters = 0;
#ifndef EIG
  node0_printf("Don't use control_eig unless compiling with -DEIG!\n");
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

  // Main measurement: PRIMME eigenvalues
  // Calculate and print smallest eigenvalues,
  // checking |D^dag D phi - lambda phi|^2
  total_iters = make_evs(Nvec, eigVec, eigVal, 1);

  // Calculate and print largest eigenvalues, for tuning RHMC
  // Don't need to compute many here...
  // (Initial eigVal/eigVec alloc moved to setup.c)
  if (Nvec > 12) {
    free(eigVal);
    for (ivec = 0; ivec < Nvec; ivec++)
      free(eigVec[ivec]);
    free(eigVec);

    Nvec = 12;
    eigVal = malloc(sizeof *eigVal * Nvec);
    eigVec = malloc(sizeof *eigVec * Nvec);
    for (ivec = 0; ivec < Nvec; ivec++) {
      eigVec[ivec] = malloc(sizeof(matrix*) * NFERMION);
      for (j = 0; j < NFERMION; j++)
        FIELD_ALLOC(eigVec[ivec][j], matrix);
    }
  }
  total_iters += make_evs(Nvec, eigVec, eigVal, -1);
  //eig();
  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
