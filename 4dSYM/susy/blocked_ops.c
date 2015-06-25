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

      for (b = 0; b < NUMLINK; b++) {
        for (j = 0; j < numK; j++) {
          // Now SUGRA with mu--nu trace subtraction
          // Symmetric and traceless by construction so ignore nu <= mu
          for (mu = 0; mu < NDIMS ; mu++) {
            for (nu = mu + 1; nu < NDIMS ; nu++) {
              tr = P[mu][a] * P[nu][b] + P[nu][a] * P[mu][b];
              OS[j] += 0.5 * tr * traceBB[j][a][b][i];
            }
          }
        }
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

  // SUGRA, averaging over six components with mu < nu
  norm = (Real)(6.0 * bl * bl * bl * bl);
  node0_printf("OS %d %d", Nsmear, block);
  for (j = 0; j < numK; j++)
    node0_printf(" %.8g", OS[j] / norm);
  node0_printf("\n");
}
// -----------------------------------------------------------------
