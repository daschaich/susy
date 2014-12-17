// -----------------------------------------------------------------
// Measure the Konishi and SUGRA scalar field correlation functions
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA correlators, projected to zero spatial momentum
void d_correlator() {
  register int i;
  register site *s;
  int a, b, c, d, mu, nu, t;
  Real tr, sub, norm;
  Real *OK[NDIMS], *OS[NDIMS][NDIMS];       // Konishi and SUGRA operators

  // Allocate and initialize Konishi and SUGRA operators
  // SUGRA will be symmetric by construction, so ignore nu < mu
  for (mu = 0; mu < NDIMS; mu++) {
    OK[mu] = malloc(nt * sizeof(*OK[mu]));
    for (nu = mu; nu < NDIMS; nu++)
      OS[mu][nu] = malloc(nt * sizeof(*OS[mu][nu]));
  }
  for (t = 0; t < nt; t++) {
    for (mu = 0; mu < NDIMS; mu++) {
      OK[mu][t] = 0.0;
      for (nu = mu; nu < NDIMS; nu++)
        OS[mu][nu][t] = 0.0;
    }
  }

  // Compute at each site B_a = U_a Udag_a - trace
  // as well as traceBB[mu][nu] = tr[B_mu(x) B_nu(x)]
  // Now stored in the site structure
  compute_Bmu();

  // Now form the zero momentum projected operators (summing across nodes)
  FORALLSITES(i, s) {
    t = s->t;
    for (a = 0; a < NUMLINK; a++) {
      for (b = 0; b < NUMLINK; b++) {
        // Four possible Konishi operators
        // All fairly easy and normalized below
        OK[0][t] += s->traceBB[a][b];
        OK[1][t] += s->traceCC[a][b];
        for (c = 0; c < NUMLINK; c++) {
          for (d = 0; d < NUMLINK; d++) {
            OK[2][t] += s->traceBB[a][b] * s->traceBB[c][d];
            OK[3][t] += s->traceCC[a][b] * s->traceCC[c][d];
          }
        }

        // Compute mu--nu trace to be subtracted from SUGRA
        sub = P[0][a] * P[0][b] * s->traceBB[a][b];
        for (mu = 1; mu < NDIMS; mu++)
          sub += P[mu][a] * P[mu][b] * s->traceBB[a][b];

        // Now SUGRA with mu--nu trace subtraction
        // Symmetric by construction so ignore nu < mu
        for (mu = 0; mu < NDIMS ; mu++) {
          OS[mu][mu][t] -= 0.25 * sub;
          for (nu = mu; nu < NDIMS ; nu++) {
            tr = P[mu][a] * P[nu][b] + P[nu][a] * P[mu][b];
            OS[mu][nu][t] += 0.5 * tr * s->traceBB[a][b];
          }
        }
      }
    }
  }
  // Normalization removed from site loop, followed by sum over nodes
  norm = (Real)(nx * ny * nz);
  for (t = 0; t < nt; t++) {
    OK[0][t] /= (5.0 * norm);
    OK[1][t] /= (5.0 * norm);
    OK[2][t] /= (25.0 * norm);
    OK[3][t] /= (25.0 * norm);
    for (mu = 0; mu < NDIMS; mu++) {
      g_doublesum(&(OK[mu][t]));
      for (nu = mu; nu < NDIMS; nu++) {
        OS[mu][nu][t] /= norm;
        g_doublesum(&(OS[mu][nu][t]));
      }
    }
  }

  // Just print operators for offline vev subtraction and correlator analysis
  // 4 Konishi operators
  for (t = 0; t < nt; t++) {
    node0_printf("KONISHI %d", t);
    for (mu = 0; mu < NDIMS; mu++)
      node0_printf(" %.8g", OK[mu][t]);
    node0_printf("\n");
  }

  // SUGRA, averaging over ten components with mu <= nu
  for (t = 0; t < nt; t++) {
    // Debugging check that trace was successfully subtracted
//    node0_printf("Trace check %d: %.4g + %.4g + %.4g + %.4g = %.4g\n",
//                 t, OS[0][0][t], OS[1][1][t], OS[2][2][t], OS[3][3][t],
//                 OS[0][0][t] + OS[1][1][t] + OS[2][2][t] + OS[3][3][t]);
    tr = 0.0;
    for (mu = 0; mu < NDIMS; mu++) {
      for (nu = mu; nu < NDIMS; nu++)
        tr += OS[mu][nu][t];
    }
    node0_printf("SUGRA %d %.8g\n", t, tr * 0.1);
  }

  for (mu = 0; mu < NDIMS; mu++) {
  free(OK[mu]);
    for (nu = mu; nu < NDIMS; nu++)
      free(OS[mu][nu]);
  }
}
// -----------------------------------------------------------------
