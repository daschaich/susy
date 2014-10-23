// -----------------------------------------------------------------
// Measure the Konishi and SUGRA scalar field correlation functions
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA correlators, projected to zero spatial momentum
void d_correlator() {
  register int i;
  register site *s;
  int mu, nu, t, tt;
  Real norm, corr, sub;
  Real *OK, *OS[NDIMS][NDIMS];    // Konishi and SUGRA operators

  // Allocate and initialize Konishi and SUGRA operators
  // SUGRA will be symmetric by construction, so ignore nu < mu
  OK = malloc(nt * sizeof(*OK));
  for (mu = 0; mu < NDIMS; mu++) {
    for (nu = mu; nu < NDIMS; nu++)
      OS[mu][nu] = malloc(nt * sizeof(*OS[mu][nu]));
  }
  for (t = 0; t < nt; t++) {
    OK[t] = 0.0;
    for (mu = 0; mu < NDIMS; mu++) {
      for (nu = mu; nu < NDIMS; nu++)
        OS[mu][nu][t] = 0.0;
    }
  }

  // Compute at each site B_a = U_a Udag_a - volume average
  // as well as traceBB[mu][nu] = tr[B_mu(x) B_nu(x)]
  // Now stored in the site structure
  compute_Bmu();

  // Now form the zero momentum projected operators (summing across nodes)
  // Don't yet normalize by (Nt / vol)
  FORALLSITES(i, s) {
    t = s->t;
    sub = 0.5 * (s->traceBB[0][0] + s->traceBB[1][1]);
    for (mu = 0; mu < NUMLINK; mu++) {
      OK[t] += s->traceBB[mu][mu];
      OS[mu][mu][t] += s->traceBB[mu][mu] - sub;
    }
    OS[0][1][t] += s->traceBB[0][1];
  }
  for (t = 0; t < nt; t++) {
    g_doublesum(&OK[t]);
    for (mu = 0; mu < NDIMS; mu++) {
      for (nu = mu; nu < NDIMS; nu++)
        g_doublesum(&OS[mu][nu][t]);
    }
  }

  // Form and print out correlators, normalized by Nt / vol^2
  // (Averaging over tt removes one factor of Nt from normalization)
  // Konishi
  norm = (Real)(nx * volume);
  for (t = 0; t <= (int)(nt / 2); t++) {
    corr = OK[0] * OK[t];
    for (tt = 1; tt < nt; tt++)
      corr += OK[tt] * OK[(t + tt) % nt];
    node0_printf("KONISHI %d %.8g\n", t, corr / norm);
  }

  // SUGRA, averaging over three components with mu <= nu
  norm = (Real)(3.0 * nx * volume);
  for (t = 0; t <= (int)(nt / 2); t++) {
    // Debugging check that trace was successfully subtracted
//    node0_printf("Trace check %d: %.4g + %.4g = %.4g\n",
//                 t, OS[0][0][t], OS[1][1][t], OS[0][0][t] + OS[1][1][t]);

    corr = 0.0;
    for (mu = 0; mu < NDIMS; mu++) {
      for (nu = mu; nu < NDIMS; nu++) {
        for (tt = 0; tt < nt; tt++)
          corr += OS[mu][nu][tt] * OS[mu][nu][(t + tt) % nt];
      }
    }
    node0_printf("SUGRA %d %.8g\n", t, corr / norm);
  }

  free(OK);
  for (mu = 0; mu < NDIMS; mu++) {
    for (nu = mu; nu < NDIMS; nu++)
      free(OS[mu][nu]);
  }
}
// -----------------------------------------------------------------
