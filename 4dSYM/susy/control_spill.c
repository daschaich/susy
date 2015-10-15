// -----------------------------------------------------------------
// Main procedure for N=4 SYM lattice printing
#define CONTROL
#include "susy_includes.h"

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
