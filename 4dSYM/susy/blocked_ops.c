// -----------------------------------------------------------------
// Konishi and SUGRA operators integrated over lattice volume
#include "susy_includes.h"

void blocked_ops(int Nsmear, int block) {
  register int i;
  register site *s;
  int a, b, j, mu, nu, bl = 2;
  Real norm, tr, OK[numK], OS[numK];    // Konishi and SUGRA operators

  // Initialize Konishi and SUGRA operators
  // SUGRA will be symmetric by construction, so ignore nu <= mu
  for (j = 0; j < numK; j++) {
    OK[j] = 0.0;
    OS[j] = 0.0;
  }

  // Compute traces of bilinears of scalar field interpolating ops
  compute_Ba();

  // Now sum operators over lattice volume
  for (a = 0; a < NUMLINK; a++) {
    FORALLSITES(i, s) {
      for (j = 0; j < numK; j++)
        OK[j] += traceBB[j][a][a][i];

      // Now SUGRA summed over ten components with a < b
      for (b = a + 1; b < NUMLINK; b++) {
        for (j = 0; j < numK; j++)
          OS[j] += traceBB[j][a][b][i];
      }
    }
  }
  for (j = 0; j < numK; j++) {
    g_doublesum(&(OK[j]));
    g_doublesum(&(OS[j]));
  }

  // Print out integrated operators
  // Have to divide by number 2^(4block) = bl^4 of blocked lattices
  for (a = 1; a < block; a++)
    bl *= 2;
  if (block <= 0)
    bl = 1;
  norm = (Real)(bl * bl * bl * bl);
  node0_printf("OK %d %d", Nsmear, block);
  for (j = 0; j < numK; j++)
    node0_printf(" %.8g", OK[j] / norm);
  node0_printf("\n");

  // SUGRA, averaging over ten components with a < b
  norm = (Real)(10.0 * bl * bl * bl * bl);
  node0_printf("OS %d %d", Nsmear, block);
  for (j = 0; j < numK; j++)
    node0_printf(" %.8g", OS[j] / norm);
  node0_printf("\n");
}
// -----------------------------------------------------------------
