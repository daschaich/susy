// -----------------------------------------------------------------
// Measure the Konishi and SUGRA scalar field correlation functions
// Use general gathers, but combine Konishi and SUGRA into single vector
// Then only need one general gather per displacement
#include "susy_includes.h"

// Define CHECK_ROT to check rotational invariance
//#define CHECK_ROT
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Map (x, y, z, t) to r on Nx x Ny x Nz x Nt A4* lattice
// Check all possible periodic shifts to find true r
Real A4map(x_in, y_in, z_in, t_in) {
  int x, y, z, t, xSq, ySq, zSq, xy, xpy, xpyz, xpypz;
  Real r = 100.0 * MAX_X, tr;

  for (x = x_in - nx; x <= x_in + nx; x += nx) {
    xSq = x * x;
    for (y = y_in - ny; y <= y_in + ny; y += ny) {
      ySq = y * y;
      xy = x * y;
      xpy = x + y;
      for (z = z_in - nz; z <= z_in + nz; z += nz) {
        zSq = z * z;
        xpyz = xpy * z;
        xpypz = xpy + z;
        for (t = t_in - nt; t <= t_in + nt; t += nt) {
          tr = sqrt((xSq + ySq + zSq + t * t) * 0.8
                    - (xy + xpyz + xpypz * t) * 0.4);
          if (tr < r)
            r = tr;
        }
      }
    }
  }
  return r;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA correlators as functions of (x, y, z, t)
void d_correlator_r() {
  register int i;
  register site *s;
  int a, b, mu, nu, index, x_dist, y_dist, z_dist, t_dist;
  int y_start, z_start, t_start, len = 11;
  int d[NDIMS] = {0, 0, 0, 0}, this_r, total_r = 0;
  int MAX_pts = 8 * MAX_X * MAX_X * MAX_X;  // Should be plenty
  int count[MAX_pts];
  Real MAX_r = 100.0 * MAX_X, tr, sub, OK = 0.0;
  Real lookup[MAX_pts], CK[MAX_pts], CS[MAX_pts];
  Real *ops = malloc(sites_on_node * len * sizeof(Real*));
  msg_tag *mtag;

  // Initialize Konishi and SUGRA operators
  // On each site, components [0, 9] of ops are the SUGRA, last is Konishi
  FORALLSITES(i, s) {
    index = i * len;
    for (a = 0; a < len; a++) {
      ops[index] = 0.0;
      index++;
    }
  }
  for (i = 0; i < MAX_pts; i++) {
    count[i] = 0;
    CK[i] = 0.0;
    CS[i] = 0.0;
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

        // Try subtracting volume average from Konishi
        OK += s->traceBB[a][b];

        // Compute mu--nu trace to be subtracted
        sub = P[0][a] * P[0][b] * s->traceBB[a][b];
        for (mu = 1; mu < NDIMS ; mu++)
          sub += P[mu][a] * P[mu][b] * s->traceBB[a][b];

        // Now SUGRA with mu--nu trace subtraction
        index = i * len;
        for (mu = 0; mu < NDIMS ; mu++) {
          ops[index] -= 0.25 * sub;
          for (nu = mu; nu < NDIMS ; nu++) {
            tr = P[mu][a] * P[nu][b] + P[nu][a] * P[mu][b];
            ops[index] += 0.5 * tr * s->traceBB[a][b];
            index++;
          }
        }
      }
    }
    index = i * len + len - 1;
    ops[index] /= 5.0;    // Konishi normalization
  }

  // Try subtracting volume average from Konishi
  OK /= (5.0  * volume);   // Remove from site loop
  g_doublesum(&OK);
  FORALLSITES(i, s) {
    index = i * len + len - 1;
    ops[index] -= OK;
  }

  // Find smallest scalar distance cut by imposing MAX_X
  for (y_dist = 0; y_dist <= MAX_X; y_dist++) {
    for (z_dist = 0; z_dist <= MAX_X; z_dist++) {
      for (t_dist = 0; t_dist <= MAX_X; t_dist++) {
        tr = A4map(MAX_X + 1, y_dist, z_dist, t_dist);
        if (tr < MAX_r)
          MAX_r = tr;
      }
    }
  }

  // Assume MAX_T >= MAX_X
  node0_printf("d_correlator_r: MAX = %d --> r < %.6g\n", MAX_X, MAX_r);

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

        // Ignore any t_dist > MAX_X even if t_dist <= MAX_T
        for (t_dist = t_start; t_dist <= MAX_X; t_dist++) {
          d[TUP] = t_dist;
          mtag = start_general_gather_field(ops, len * sizeof(Real),
                                            d, EVENANDODD, gen_pt[0]);

          // Figure out scalar distance while general gather runs
          tr = A4map(x_dist, y_dist, z_dist, t_dist);
          wait_general_gather(mtag);

          if (tr > MAX_r - 1.0e-6) {
            cleanup_general_gather(mtag);
            continue;
          }

          // Combine four-vectors with same scalar distance
          this_r = -1;
          for (i = 0; i < total_r; i++) {
            if (fabs(tr - lookup[i]) < 1.0e-6) {
              this_r = i;
              break;
            }
          }
          if (this_r < 0) {   // Add new scalar distance to lookup table
            lookup[total_r] = tr;
            this_r = total_r;
            total_r++;
            if (total_r >= MAX_pts) {
              node0_printf("d_correlator_r: MAX_pts %d too small\n", MAX_pts);
              terminate(1);
            }
          }
          count[this_r]++;

          // Konishi
          tr = 0.0;
          FORALLSITES(i, s) {
            index = i * len + len - 1;
            tr += ops[index] * ((Real *)gen_pt[0][i])[len - 1];
          }
          g_doublesum(&tr);
          CK[this_r] += tr;
#ifdef CHECK_ROT
          // Potentially useful to check rotational invariance
          node0_printf("CORR_K %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist,
                       tr * 1.0 / (Real)volume);
#endif

          // SUGRA, averaging over ten components with mu <= nu
          a = 0;
          tr = 0.0;
          for (mu = 0; mu < NDIMS ; mu++) {
            for (nu = mu; nu < NDIMS ; nu++) {
              FORALLSITES(i, s) {
                index = i * len + a;
                tr += ops[index] * ((Real *)gen_pt[0][i])[a];
              }
              a++;
            }
          }
          g_doublesum(&tr);
          CS[this_r] += tr;
#ifdef CHECK_ROT
          // Potentially useful to check rotational invariance
          node0_printf("CORR_S %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist,
                       tr * 0.1 / (Real)volume);
#endif
          cleanup_general_gather(mtag);
        } // t_dist
      } // z_dist
    } // y dist
  } // x dist

  // Now cycle through unique scalar distances and print results
  // Won't be sorted, but this is easy to do offline
  for (i = 0; i < total_r; i++) {
    tr = 1.0 / (Real)(count[i] * volume);
    node0_printf("CORR_K %d %.6g %.6g\n", i, lookup[i], CK[i] * tr);
  }
  for (i = 0; i < total_r; i++) {
    tr = 0.1 / (Real)(count[i] * volume);
    node0_printf("CORR_S %d %.6g %.6g\n", i, lookup[i], CS[i] * tr);
  }

  free(ops);
}
// -----------------------------------------------------------------
