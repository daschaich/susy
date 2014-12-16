// -----------------------------------------------------------------
// Measure the Konishi and SUGRA scalar field correlation functions
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA correlators, projected to zero spatial momentum
void d_correlator() {
  register int i;
  register site *s;
  int a, b, c, d, mu, nu, t, tt;
  Real norm, tr, sub;
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
  // Don't yet normalize by (Nt / vol)
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
  for (t = 0; t < nt; t++) {
    OK[0][t] /= 5.0;   // Remove from site loop
    OK[1][t] /= 5.0;
    OK[2][t] /= 25.0;
    OK[3][t] /= 25.0;
    for (mu = 0; mu < NDIMS; mu++) {
      g_doublesum(&(OK[mu][t]));
      for (nu = mu; nu < NDIMS; nu++)
        g_doublesum(&(OS[mu][nu][t]));
    }
  }

  // Try subtracting volume average from Konishi
  for (mu = 0; mu < NDIMS; mu++) {
    tr = OK[mu][0];
    for (t = 1; t < nt; t++)
      tr += OK[mu][t];
    tr /= nt;
    for (t = 0; t < nt; t++)
      OK[mu][t] -= tr;
  }

  // Form and print out correlators, normalized by Nt / vol^2
  // (Averaging over tt removes one factor of Nt from normalization)
  // Konishi
  // 16 Konishi correlators
  norm = (Real)(nx * ny * nz * volume);
  for (t = 0; t <= (int)(nt / 2); t++) {
    node0_printf("KONISHI %d", t);
    for (mu = 0; mu < NDIMS; mu++) {
      for (nu = 0; nu < NDIMS; nu++) {
        tr = OK[mu][0] * OK[nu][t];
        for (tt = 1; tt < nt; tt++)
          tr += OK[mu][tt] * OK[nu][(t + tt) % nt];
        node0_printf(" %.8g", tr / norm);
      }
    }
    node0_printf("\n");
  }

  // SUGRA, averaging over ten components with mu <= nu
  norm = (Real)(10.0 * nx * ny * nz * volume);
  for (t = 0; t <= (int)(nt / 2); t++) {
    // Debugging check that trace was successfully subtracted
//    node0_printf("Trace check %d: %.4g + %.4g + %.4g + %.4g = %.4g\n",
//                 t, OS[0][0][t], OS[1][1][t], OS[2][2][t], OS[3][3][t],
//                 OS[0][0][t] + OS[1][1][t] + OS[2][2][t] + OS[3][3][t]);
    tr = 0.0;
    for (mu = 0; mu < NDIMS; mu++) {
      for (nu = mu; nu < NDIMS; nu++) {
        for (tt = 0; tt < nt; tt++)
          tr += OS[mu][nu][tt] * OS[mu][nu][(t + tt) % nt];
      }
    }
    node0_printf("SUGRA %d %.8g\n", t, tr / norm);
  }

  for (mu = 0; mu < NDIMS; mu++) {
  free(OK[mu]);
    for (nu = mu; nu < NDIMS; nu++)
      free(OS[mu][nu]);
  }
}
// -----------------------------------------------------------------
