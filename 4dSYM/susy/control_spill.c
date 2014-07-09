// -----------------------------------------------------------------
// Main procedure for N=4 SYM lattice printing
#define CONTROL
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt, s, x, y, z, t, mu, i, j;
  Real re, im;
  double dssplaq, dstplaq;

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();
  setup_lambda();
  setup_PtoP();
  setup_FQ();

  // Load input and run (loop removed)
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }

  // Serial code!
  if (this_node != 0) {
    node0_printf("ERROR: run this thing in serial!\n");
    terminate(1);
  }

  // Optionally gauge fix to Coulomb gauge
  if (fixflag == COULOMB_GAUGE_FIX) {
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
  }
  else if (fixflag == NO_GAUGE_FIX) { // Braces suppress compiler warning
    node0_printf("Gauge fixing skipped\n");
  }
  else {
    node0_printf("ERROR: only COULOMB_GAUGE_FIX ");
    node0_printf("and NO_GAUGE_FIX supported\n");
    terminate(1);
  }

  // Spill lattice in format expected by serial code
  // (cf. read_in.cpp, loop_over_lattice and << overloading for Umatrix
  for (t = 0; t < nt; t++) {
    for (z = 0; z < nz; z++) {
      for (y = 0; y < ny; y++) {
        for (x = 0; x < nx; x++) {
          s = node_index(x, y, z, t);
          for (mu = 0; mu < NUMLINK; mu++) {
            for (i = 0; i < NCOL; i++) {
              for (j = 0; j < NCOL; j++) {
                re = lattice[s].linkf[mu].e[i][j].real;
                im = lattice[s].linkf[mu].e[i][j].imag;
                printf("%g\t%g\t", re, im);
              }
            }
            printf("\n");
          }
        }
      }
    }
  }
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------
