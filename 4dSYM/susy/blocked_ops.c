// -----------------------------------------------------------------
// Konishi and SUGRA operators averaged over lattice volume
// Blocking actually doesn't affect the computation!
#include "susy_includes.h"

void blocked_ops(int block) {
  register int i;
  register site *s;
  int a, b, c, d, mu, nu, bl = 2;
  Real norm, tr, OK[NDIMS], OS[NDIMS][NDIMS]; // Konishi and SUGRA operators

  // Initialize Konishi and SUGRA operators
  // SUGRA will be symmetric by construction, so ignore nu <= mu
  for (mu = 0; mu < NDIMS; mu++) {
    OK[mu] = 0.0;
    for (nu = mu + 1; nu < NDIMS; nu++)
      OS[mu][nu] = 0.0;
  }

  // Compute at each site B_a = U_a Udag_a - volume average
  // as well as traceBB[mu][nu] = tr[B_mu(x) B_nu(x)]
  // Now stored in the site structure
  compute_Bmu();

  // Now sum operators over lattice volume
  FORALLSITES(i, s) {
    for (a = 0; a < NUMLINK; a++) {
      for (b = 0; b < NUMLINK; b++) {
        OK[0] += s->traceBB[a][b];        // Konishi is easy
        for (c = 0; c < NUMLINK; c++) {
          OK[1] += s->traceBBB[a][b][c];  // Another Konishi
          for (d = 0; d < NUMLINK; d++) {
            OK[2] += s->traceBBBB[a][b][c][d];    // A third Konishi
            OK[3] += s->traceBB[a][b] * s->traceBB[c][d];     // Double trace
          }
        }

        // Now SUGRA with mu--nu trace subtraction
        // Symmetric and traceless by construction so ignore nu <= mu
        for (mu = 0; mu < NDIMS ; mu++) {
          for (nu = mu + 1; nu < NDIMS ; nu++) {
            tr = P[mu][a] * P[nu][b] + P[nu][a] * P[mu][b];
            OS[mu][nu] += 0.5 * tr * s->traceBB[a][b];
          }
        }
      }
    }

    // Check that OK[1] vanishes site by site while OK[3] = 2OK[2]
    if (OK[1] > IMAG_TOL) {
      printf("node%d WARNING: sum(tr(BBB)) = %.4g after i = %d\n",
             this_node, OK[1], i);
    }
    tr = 2.0 * OK[2] - OK[3];
    if (tr > 1.0e-8) {
      printf("node%d: 2sum(tr(BBBB)) - sum(tr(BB)*tr(BB)) = %.4g - %.4g ",
             this_node, 2.0 * OK[2], OK[3]);
      printf("= %.4g after i = %d\n", tr, i);
    }
  }
  OK[0] /= 5.0;   // Remove from site loop
  OK[1] /= (5.0 * sqrt(5.0));
  OK[2] /= 25.0;
  OK[3] /= 25.0;
  for (mu = 0; mu < NDIMS; mu++) {
    g_doublesum(&(OK[mu]));
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
  node0_printf("OK %d", block);
  for (mu = 0; mu < NDIMS; mu++)
    node0_printf(" %.8g", OK[mu] / norm);
  node0_printf("\n");

  // SUGRA, averaging over six components with mu < nu
  norm = (Real)(6.0 * bl * bl * bl * bl);
  tr = 0.0;
  for (mu = 0; mu < NDIMS; mu++) {
    for (nu = mu + 1; nu < NDIMS; nu++)
      tr += OS[mu][nu];
  }
  node0_printf("OS %d %.8g\n", block, tr / norm);
}
// -----------------------------------------------------------------
