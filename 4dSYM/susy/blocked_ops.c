// -----------------------------------------------------------------
// Konishi and SUGRA operators averaged over lattice volume
// Blocking actually doesn't affect the computation!
#include "susy_includes.h"

void blocked_ops(int block) {
  register int i;
  register site *s;
  int a, b, mu, nu;
  Real norm, tr, sub;
  Real OK, OS[NDIMS][NDIMS];    // Konishi and SUGRA operators

  // Initialize Konishi and SUGRA operators
  // SUGRA will be symmetric by construction, so ignore nu < mu
  // For just the operator could also ignore diagonal,
  // but for now let's check that it vanishes as it should
  OK = 0.0;
  for (mu = 0; mu < NDIMS; mu++) {
    for (nu = mu; nu < NDIMS; nu++)
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
        OK += s->traceBB[a][b];    // Konishi is easy

        // Compute mu--nu trace to be subtracted
        sub = P[0][a] * P[0][b] * s->traceBB[a][b];
        for (mu = 1; mu < NDIMS; mu++)
          sub += P[mu][a] * P[mu][b] * s->traceBB[a][b];

        // Now SUGRA with mu--nu trace subtraction
        // Symmetric by construction so ignore nu < mu
        for (mu = 0; mu < NDIMS ; mu++) {
          OS[mu][mu] -= 0.25 * sub;
          for (nu = mu; nu < NDIMS ; nu++) {
            tr = P[mu][a] * P[nu][b] + P[mu][a] * P[nu][b];
            OS[mu][nu] += 0.5 * tr * s->traceBB[a][b];
          }
        }
      }
    }
  }
  OK /= 5.0;   // Remove from site loop
  g_doublesum(&OK);
  for (mu = 0; mu < NDIMS; mu++) {
    for (nu = mu; nu < NDIMS; nu++)
      g_doublesum(&OS[mu][nu]);
  }

  // Print out averaged operators, averaging over all 6 SUGRA components
  norm = (Real)(volume);
  node0_printf("OK %d %.8g\n", block, OK / norm);

  // SUGRA, averaging over six components with mu < nu
  norm = (Real)(6.0 * nx * ny * nz * volume);
  // Debugging check that trace was successfully subtracted
  tr = OS[0][0] + OS[1][1] + OS[2][2] + OS[3][3];
  if (tr > IMAG_TOL) {
    node0_printf("WARNING: Non-zero trace %.4g + %.4g + %.4g + %.4g = %.4g\n",
                 OS[0][0], OS[1][1], OS[2][2], OS[3][3], tr);
  }
  tr = 0.0;
  for (mu = 0; mu < NDIMS; mu++) {
    for (nu = mu + 1; nu < NDIMS; nu++)
      tr += OS[mu][nu];
  }
  node0_printf("OS %d %.8g\n", block, tr / norm);
}
// -----------------------------------------------------------------
