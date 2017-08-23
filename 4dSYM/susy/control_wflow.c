// -----------------------------------------------------------------
// Main procedure for N=4 SYM Wilson flow 
#define CONTROL
#include "susy_includes.h"

int main(int argc, char *argv[]) {

#ifndef WFLOW
  node0_printf("Don't use control_wflow unless compiling with -DWFLOW!\n");
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
  for (dir = XUP; dir < NUMLINK; dir++)
    node0_printf(" %.6g", linktr[dir]);
  node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);
  node0_printf("FLINK_DET");
  for (dir = XUP; dir < NUMLINK; dir++)
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

  // Main Wilson flow routine : 
  
  wflow();  
  
  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  g_sync();         // Needed by at least some clusters
  return 0; 
}
