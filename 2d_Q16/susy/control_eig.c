// -----------------------------------------------------------------
// Main procedure for N=(8,8) SYM eigenvalues
#define CONTROL
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt, dir;
  double dplaq, dtime, plpMod = 0.0;
  double linktr[NUMLINK], linktr_ave, linktr_width;
  double link_det[NUMLINK], det_ave, det_width;
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

  // Check: compute initial plaquette and bosonic action
  plaquette(&dplaq);
  node0_printf("START %.8g ", dplaq);
  dplaq = gauge_action(NODET);
  node0_printf("%.8g\n", dplaq / (double)volume);

  // Do "local" measurements to check configuration
  // Tr[Udag.U] / N
  linktr_ave = link_trace(linktr, &linktr_width,
                          link_det, &det_ave, &det_width);
  node0_printf("FLINK");
  for (dir = XUP; dir < NUMLINK; dir++)
    node0_printf(" %.6g", linktr[dir]);
  node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);
  node0_printf("FLINK_DET");
  for (dir = XUP; dir < NUMLINK; dir++)
    node0_printf(" %.6g", link_det[dir]);
  node0_printf(" %.6g %.6g\n", det_ave, det_width);

  // Polyakov loop and plaquette measurements
  // Format: GMES Re(Polyakov) Im(Poyakov) cg_iters plaq
  plp = ploop(TUP, NODET, &plpMod);
  plaquette(&dplaq);
  node0_printf("GMES %.8g %.8g 0 %.8g ",
               plp.real, plp.imag, dplaq);

  // Bosonic action (printed twice by request)
  // Might as well spit out volume average of Polyakov loop modulus
  dplaq = gauge_action(NODET) / (double)volume;
  node0_printf("%.8g ", dplaq);
  node0_printf("%.8g\n", plpMod);
  node0_printf("BACTION %.8g\n", dplaq);

  // Main measurement: PRIMME eigenvalues
  // Allocate eigenvectors
  eigVal = malloc(Nvec * sizeof(*eigVal));
  eigVec = malloc(Nvec * sizeof(*eigVec));
  for (ivec = 0; ivec < Nvec; ivec++)
    eigVec[ivec] = malloc(sites_on_node * sizeof(Twist_Fermion));

  // Calculate and print smallest eigenvalues,
  // checking |D^dag D phi - lambda phi|^2
  total_iters = make_evs(Nvec, eigVec, eigVal, 1);

  // Check matrix elements of D with DDdag eigenmodes
  // The eigenvalues should be paired, with each pair producing
  // positive/negative matrix elements
  // In principle, one could tighten eig_tol until all pairs are found
  // For now we just print them all out to check offline
  check_Dmat(Nvec, eigVec);

  // Calculate and print largest eigenvalues, for tuning RHMC
  // Don't need to compute many here...
  if (Nvec > 12) {
    Nvec = 12;
    free(eigVal);
    eigVal = malloc(Nvec * sizeof(*eigVal));

    for (ivec = 0; ivec < Nvec; ivec++)
      free(eigVec[ivec]);
    free(eigVec);

    eigVec = malloc(Nvec * sizeof(*eigVec));
    for (ivec = 0; ivec < Nvec; ivec++)
      eigVec[ivec] = malloc(sites_on_node * sizeof(Twist_Fermion));
  }
  total_iters += make_evs(Nvec, eigVec, eigVal, -1);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);

  free(eigVal);
  for (ivec = 0; ivec < Nvec; ivec++)
    free(eigVec[ivec]);
  free(eigVec);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
