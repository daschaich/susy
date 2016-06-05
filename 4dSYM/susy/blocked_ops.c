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
  compute_Ba();

  // Now sum operators over lattice volume
  FORALLSITES(i, s) {
    FORALLDIR(a) {
      for (j = 0; j < N_K; j++) {
        OK[j] += traceBB[j][a][a][i];   // Konishi

        FORALLDIR(b) {
          // Add Tr phi_6^2 = (1 / 5) Tr \sum_{a, b} A_a A_b
          // to cancel SUGRA mixing in log-polar Konishi
          if (j == 0)
            OK[j] += 0.2 * traceAA[a][b][i];

          // Now SUGRA, averaged over 20 off-diagonal components
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

  // Print each operator summed over volume
  // Have to divide by number of blocked lattices: 2^(4block) = bl^4
  // Format: TAG  smearing  block  a  op[a]  subtracted[a]
  for (j = 1; j < block; j++)
    bl *= 2;
  if (block <= 0)
    bl = 1;
  norm = (Real)(bl * bl * bl * bl);
  for (j = 0; j < N_K; j++)    // Konishi
    node0_printf("OK %d %d %d %.8g\n", Nsmear, block, j, OK[j] / norm);
  for (j = 0; j < N_K; j++)    // SUGRA
    node0_printf("OS %d %d %d %.8g\n", Nsmear, block, j, OS[j] / norm);
}
// -----------------------------------------------------------------
