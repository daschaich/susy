// -----------------------------------------------------------------
// Measure the Konishi and SUGRA scalar field correlation functions
// Use general gathers, but combine Konishi and SUGRA into single vector
// Then only need one general gather per displacement
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA correlators as functions of (x, y, z, t)
void d_correlator_r() {
  register int i;
  register site *s;
  int a, b, mu, nu, index, x_dist, y_dist, z_dist, t_dist;
  int y_start, z_start, t_start, len = 11;
  int d[NDIMS] = {0, 0, 0, 0};
  Real tr, corr, sub;
  Real normK = 1.0 / (Real)volume, normS = 0.1 / (Real)volume;
  Real *ops = malloc(sites_on_node * len * sizeof(Real*));
  msg_tag *mtag;

  node0_printf("d_correlator_r: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  // Initialize Konishi and SUGRA operators
  // On each site, components [0, 9] of ops are the SUGRA, last is Konishi
  FORALLSITES(i, s) {
    index = i * len;
    for (a = 0; a < len; a++) {
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
    for (a = 0; a < NUMLINK; a++) {
      for (b = 0; b < NUMLINK; b++) {
        // Konishi is easy -- normalized below
        index = i * len + len - 1;
        ops[index] += s->traceBB[a][b];

        // Compute mu--nu trace to be subtracted
        sub = P[0][a] * P[0][b] * s->traceBB[a][b];
        for (mu = 1; mu < NDIMS ; mu++)
          sub += P[mu][a] * P[mu][b] * s->traceBB[a][b];

        // Now SUGRA with mu--nu trace subtraction
        index = i * len;
        for (mu = 0; mu < NDIMS ; mu++) {
          ops[index] -= 0.25 * sub;
          for (nu = mu; nu < NDIMS ; nu++) {
            tr = P[mu][a] * P[nu][b] + P[mu][a] * P[nu][b];
            ops[index] += 0.5 * tr * s->traceBB[a][b];
            index++;
          }
        }
      }
    }
    index = i * len + len - 1;
    ops[index] /= 5.0;    // Konishi normalization
  }

  // Construct and print correlators
  // Use general gathers, at least for now
  for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
    d[XUP] = x_dist;

    // Don't need negative y_dist when x_dist = 0
    if (x_dist > 0)
      y_start = -MAX_X;
    else
      y_start = 0;

    for (y_dist = y_start; y_dist <= MAX_X; y_dist++) {
      d[YUP] = y_dist;

      // Don't need negative z_dist when both x, y non-positive
      if (x_dist > 0 || y_dist > 0)
        z_start = -MAX_X;
      else
        z_start = 0;

      for (z_dist = z_start; z_dist <= MAX_X; z_dist++) {
        d[ZUP] = z_dist;

        // Don't need negative t_dist when x, y and z are all non-positive
        if (x_dist > 0 || y_dist > 0 || z_dist > 0)
          t_start = -MAX_X;
        else
          t_start = 0;

        for (t_dist = t_start; t_dist <= MAX_T; t_dist++) {
          d[TUP] = t_dist;
          mtag = start_general_gather_field(ops, len * sizeof(Real),
                                            d, EVENANDODD, gen_pt[0]);
          wait_general_gather(mtag);

          // Konishi
          corr = 0.0;
          FORALLSITES(i, s) {
            index = i * len + len - 1;
            corr += ops[index] * ((Real *)gen_pt[0][i])[len - 1];
          }
          g_doublesum(&corr);
          node0_printf("CORR_K %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist, corr * normK);

          // SUGRA, averaging over ten components with mu <= nu
          a = 0;
          corr = 0.0;
          for (mu = 0; mu < NDIMS ; mu++) {
            for (nu = mu; nu < NDIMS ; nu++) {
              FORALLSITES(i, s) {
                index = i * len + a;
                corr += ops[index] * ((Real *)gen_pt[0][i])[a];
              }
              a++;
            }
          }
          g_doublesum(&corr);
          node0_printf("CORR_S %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist, corr * normS);
          cleanup_general_gather(mtag);
        } // t_dist
      } // z_dist
    } // y dist
  } // x dist
  free(ops);
}
// -----------------------------------------------------------------
