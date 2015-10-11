// -----------------------------------------------------------------
// Konishi and SUGRA operators integrated over lattice volume
#include "susy_includes.h"

void blocked_ops(int Nsmear, int block) {
  register int i;
  register site *s;
  int a, b, j, bl = 2;
  double norm, OK[N_K], OS[N_K];

  // Initialize Konishi and SUGRA operators
  for (j = 0; j < N_K; j++) {
    OK[j] = 0.0;
    OS[j] = 0.0;
  }

  // Compute traces of bilinears of scalar field interpolating ops
  // Needs to know how much to shift hermitian part of paolor decomposition
  for (j = 1; j < block; j++)
    bl *= 2;
  if (block <= 0)
    bl = 1;
  compute_Ba(bl);

  // Now sum operators over lattice volume
  FORALLSITES(i, s) {
    for (a = 0; a < NUMLINK; a++) {
      for (j = 0; j < N_K; j++) {
        OK[j] += traceBB[j][a][a][i];   // Konishi

        // Now SUGRA, averaged over 20 off-diagonal components
        for (b = 0; b < NUMLINK; b++) {
          if (a == b)
            continue;
          OS[j] += 0.05 * traceBB[j][a][b][i];
        }
      }
    }
  }
  for (j = 0; j < N_K; j++) {
    g_doublesum(&(OK[j]));
    g_doublesum(&(OS[j]));
  }

  // Print operators summed over volume
  // Have to divide by number of blocked lattices: 2^(4block) = bl^4
  norm = (Real)(bl * bl * bl * bl);
  node0_printf("OK %d %d", Nsmear, block);    // Konishi
  for (j = 0; j < N_K; j++)
    node0_printf(" %.8g", OK[j] / norm);
  node0_printf("\n");
  node0_printf("OS %d %d", Nsmear, block);    // SUGRA
  for (j = 0; j < N_K; j++)
    node0_printf(" %.8g", OS[j] / norm);
  node0_printf("\n");
}
// -----------------------------------------------------------------
