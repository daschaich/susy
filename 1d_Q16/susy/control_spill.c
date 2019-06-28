// -----------------------------------------------------------------
// Main procedure for BFSS/BMN lattice printing
#define CONTROL
#include "susy_includes.h"

int main(int argc, char *argv[]) {
  int prompt, s, t, i, j;
  Real re, im;

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();
  setup_lambda();

  // Load input and run (loop removed)
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }

  // Serial code!
  if (this_node != 0) {
    printf("ERROR: run this thing in serial!\n");
    terminate(1);
  }

  // Spill lattice in format expected by serial code
  // (cf. read_in.cpp, loop_over_lattice and << overloading for Umatrix)
  for (t = 0; t < nt; t++) {
    s = node_index(t);
    for (i = 0; i < NCOL; i++) {
      for (j = 0; j < NCOL; j++) {
        re = lattice[s].link.e[i][j].real;
        im = lattice[s].link.e[i][j].imag;
        printf("%g\t%g\t", re, im);
      }
    }
    printf("\n");
  }
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------
