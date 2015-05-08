// -----------------------------------------------------------------
// Konishi and SUGRA operators integrated over lattice volume
#include "susy_includes.h"

void blocked_ops(int Nstout, int block) {
  register int i;
  register site *s;
  int a, b, mu, nu, bl = 2, numK = 2;
  Real norm, tr, OK[numK], OS[NDIMS][NDIMS]; // Konishi and SUGRA operators

  // Initialize Konishi and SUGRA operators
  // SUGRA will be symmetric by construction, so ignore nu <= mu
  for (mu = 0; mu < numK; mu++)
    OK[mu] = 0.0;
  for (mu = 0; mu < NDIMS; mu++) {
    for (nu = mu + 1; nu < NDIMS; nu++)
      OS[mu][nu] = 0.0;
  }

  // Compute at each site B_a = U_a Udag_a - volume average
  // as well as traceBB[a][b] = tr[B_a(x) B_b(x)]
  // Now stored in the site structure
  compute_Ba();

  // Now sum operators over lattice volume
  for (a = 0; a < NUMLINK; a++) {
    FORALLSITES(i, s) {
      OK[0] += traceBB[a][a][i];
      for (b = 0; b < NUMLINK; b++) {
        OK[1] += traceBB[a][a][i] * traceBB[b][b][i];

        // Now SUGRA with mu--nu trace subtraction
        // Symmetric and traceless by construction so ignore nu <= mu
        for (mu = 0; mu < NDIMS ; mu++) {
          for (nu = mu + 1; nu < NDIMS ; nu++) {
            tr = P[mu][a] * P[nu][b] + P[nu][a] * P[mu][b];
            OS[mu][nu] += 0.5 * tr * traceBB[a][b][i];
          }
        }
      }
    }
  }
  for (mu = 0; mu < numK; mu++)
    g_doublesum(&(OK[mu]));
  for (mu = 0; mu < NDIMS; mu++) {
    for (nu = mu + 1; nu < NDIMS; nu++)
      g_doublesum(&(OS[mu][nu]));
  }

  // Print out integrated operators
  // Have to divide by number 2^(4block) = bl^4 of blocked lattices
  for (a = 1; a < block; a++)
    bl *= 2;
  if (block <= 0)
    bl = 1;
  norm = (Real)(bl * bl * bl * bl);
  node0_printf("OK %d %d", Nstout, block);
  for (mu = 0; mu < numK; mu++)
    node0_printf(" %.8g", OK[mu] / norm);
  node0_printf("\n");

  // SUGRA, averaging over six components with mu < nu
  norm = (Real)(6.0 * bl * bl * bl * bl);
  tr = 0.0;
  for (mu = 0; mu < NDIMS; mu++) {
    for (nu = mu + 1; nu < NDIMS; nu++)
      tr += OS[mu][nu];
  }
  node0_printf("OS %d %d %.8g\n", Nstout, block, tr / norm);
}
// -----------------------------------------------------------------
