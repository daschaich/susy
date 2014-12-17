// -----------------------------------------------------------------
// Measure the Konishi and SUGRA scalar field correlation functions
// Use general gathers, but combine Konishi and SUGRA into single vector
// Then only need one general gather per displacement
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA correlators as functions of (x, t)
void d_correlator_r() {
  register int i;
  register site *s;
  int mu, nu, index, x_dist, t_dist, t_start, len = 4, d[NDIMS] = {0, 0};
  Real corr, sub, OK = 0.0, normK = 1.0 / (Real)volume, normS = normK / 3.0;
  Real *ops = malloc(sites_on_node * len * sizeof(Real*));
  msg_tag *mtag;

  node0_printf("d_correlator_r: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  // Initialize Konishi and SUGRA operators
  // On each site, components [0, 2] of ops are the SUGRA, last is Konishi
  FORALLSITES(i, s) {
    index = i * len;
    for (mu = 0; mu < len; mu++) {
      ops[index] = 0.0;
      index++;
    }
  }

  // Compute at each site B_a = U_a Udag_a - volume average
  // as well as traceBB[mu][nu] = tr[B_mu(x) B_nu(x)]
  // Now stored in the site structure
  compute_Bmu();

  // Construct the operators
  FORALLSITES(i, s) {
    for (mu = 0; mu < NUMLINK; mu++) {
      // Konishi is last entry
      index = i * len + len - 1;
      for (nu = 0; nu < NUMLINK; nu++) {
        ops[index] += s->traceBB[mu][nu];
        OK += s->traceBB[mu][nu];
      }
    }

    // Subtract trace from first and third SUGRA components
    sub = 0.5 * (s->traceBB[0][0] + s->traceBB[1][1]);
    index = i * len;
    ops[index] += s->traceBB[0][0] - sub;
    index = i * len + 2;
    ops[index] += s->traceBB[1][1] - sub;
    index = i * len + 1;
    ops[index] += s->traceBB[0][1];
  }

  // Try subtracting volume average from Konishi
  OK /= (Real)(volume);
  g_doublesum(&OK);
  FORALLSITES(i, s) {
    index = i * len + len - 1;
    ops[index] -= OK;
  }

  // Construct and print correlators
  // Use general gathers, at least for now
  for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
    d[XUP] = x_dist;

    // Don't need negative t_dist when x_dist = 0
    if (x_dist > 0)
      t_start = -MAX_X;
    else
      t_start = 0;

    for (t_dist = t_start; t_dist <= MAX_T; t_dist++) {
      d[TUP] = t_dist;
      mtag = start_general_gather_field(ops, len * sizeof(Real), d,
                                        EVENANDODD, gen_pt[0]);
      wait_general_gather(mtag);

      // Konishi
      corr = 0.0;
      FORALLSITES(i, s) {
        index = i * len + len - 1;
        corr += ops[index] * ((Real *)gen_pt[0][i])[len - 1];
      }
      g_doublesum(&corr);
      node0_printf("CORR_K %d %d %.6g\n", x_dist, t_dist, corr * normK);

      // SUGRA, averaging over three components with mu <= nu
      corr = 0.0;
      for (mu = 0; mu < 3; mu++) {
        FORALLSITES(i, s) {
          index = i * len + mu;
          corr += ops[index] * ((Real *)gen_pt[0][i])[mu];
        }
      }
      g_doublesum(&corr);
      node0_printf("CORR_S %d %d %.6g\n", x_dist, t_dist, corr * normS);
      cleanup_general_gather(mtag);
    } // t_dist
  } // x dist
  free(ops);
}
// -----------------------------------------------------------------
