// -----------------------------------------------------------------
// Main procedure for N=4 SYM lattice printing
#define CONTROL
#include "susy_includes.h"

int main(int argc, char *argv[]) {
  int prompt, s, x, y, t, mu, i, j;
  Real re, im;

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();
  setup_lambda();
  setup_PtoV();
  setup_FQ();

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
  // Also need to add first line: nx\t nt\t kappa\t f_eps\t N
  for (t = 0; t < nt; t++) {
    for (y = 0; y < ny; y++) {
      for (x = 0; x < nx; x++) {
        s = node_index(x, y, z, t);
        FORALLDIR(mu) {
          for (i = 0; i < NCOL; i++) {
            for (j = 0; j < NCOL; j++) {
              re = lattice[s].link[mu].e[i][j].real;
              im = lattice[s].link[mu].e[i][j].imag;
              printf("%g\t%g\t", re, im);
            }
          }
          printf("\n");
        }
      }
    }
  }
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------
