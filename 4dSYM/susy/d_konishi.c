// -----------------------------------------------------------------
#include "susy_includes.h"

// Measure the Konishi and SUGRA operators on each timeslice
void d_konishi() {
  register int i;
  register site *s;
  int a, b, mu, nu, t, j;
  Real tr, norm;
  double *OK[numK], *OS[numK];      // Konishi and SUGRA operators

  // Allocate and initialize Konishi and SUGRA operators
  // SUGRA will be symmetric by construction, so ignore nu < mu
  for (j = 0; j < numK; j++) {
    OK[j] = malloc(nt * sizeof(*OK[j]));
    OS[j] = malloc(nt * sizeof(*OS[j]));
  }
  for (t = 0; t < nt; t++) {
    for (j = 0; j < numK; j++) {
      OK[j][t] = 0.0;
      OS[j][t] = 0.0;
    }
  }

  // Compute traces of bilinears of scalar field interpolating ops
  compute_Ba();

  // Now form the zero momentum projected operators (summing across nodes)
  // Average SUGRA over ten components with mu <= nu
  for (a = 0; a < NUMLINK; a++) {
    FORALLSITES(i, s) {
      t = s->t;
      for (j = 0; j < numK; j++)
        OK[j][t] += traceBB[j][a][a][i];

      for (b = 0; b < NUMLINK; b++) {
        for (j = 0; j < numK; j++) {
          // Now SUGRA with mu--nu trace subtraction
          // Symmetric and traceless by construction so ignore nu < mu
          for (mu = 0; mu < NDIMS; mu++) {
            for (nu = mu + 1; nu < NDIMS; nu++) {
              tr = P[mu][a] * P[nu][b] + P[nu][a] * P[mu][b];
              OS[j][t] += 0.5 * tr * traceBB[j][a][b][i];
            }
          }
        }
      }
    }
  }
  // Normalization removed from site loop, followed by sum over nodes
  norm = (Real)(nx * ny * nz);
  for (t = 0; t < nt; t++) {
    for (j = 0; j < numK; j++) {
      OK[j][t] /= norm;
      g_doublesum(&(OK[j][t]));

      OS[j][t] /= (6.0 * norm);
      g_doublesum(&(OS[j][t]));
    }
  }

  // Just print operators for offline vev subtraction and correlator analysis
  for (t = 0; t < nt; t++) {
    node0_printf("KONISHI %d", t);
    for (j = 0; j < numK; j++)
      node0_printf(" %.16g", OK[j][t] - vevK[j]);
    node0_printf("\n");
  }

  for (t = 0; t < nt; t++) {
    node0_printf("SUGRA %d", t);
    for (j = 0; j < numK; j++)
      node0_printf(" %.16g", OS[j][t]);
    node0_printf("\n");
  }

  for (j = 0; j < numK; j++) {
    free(OK[j]);
    free(OS[j]);
  }
}
// -----------------------------------------------------------------
