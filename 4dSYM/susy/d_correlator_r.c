// -----------------------------------------------------------------
// Measure the Konishi and SUGRA correlation functions
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
// Both Konishi and SUGRA correlators as functions of r
// For gathering, all operators live in array of len = 4N_K
// vevK[N_K] are Konishi ensemble averages; assume these vanish for SUGRA
// volK[N_K] and volS[N_K] are Konishi and SUGRA volume averages
// Check to make sure the latter have been set by d_konishi()
void d_correlator_r() {
  register int i;
  register site *s;
  int a, b, j, mu, nu, index, x_dist, y_dist, z_dist, t_dist;
  int y_start, z_start, t_start, len = 4 * N_K, subs;
  int disp[NDIMS] = {0, 0, 0, 0}, this_r, total_r = 0;
  int MAX_pts = 8 * MAX_X * MAX_X * MAX_X, count[MAX_pts];
  Real MAX_r = 100.0 * MAX_X, tr;   // To be overwritten
  Real lookup[MAX_pts], CK[MAX_pts][len * N_K], CS[MAX_pts][len * N_K]; // !!!
  Real *ops = malloc(sites_on_node * len * sizeof(*ops));
  msg_tag *mtag;

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

  // Subtract either ensemble average or volume average (respectively)
  // while initializing operators and correlators
  // First half are Konishi, second half are SUGRA
  // Make sure vevK and vevS have been set up properly
  if (volK[0] < 0.0) {
    node0_printf("ERROR: Improperly initialized vacuum subtraction\n");
    terminate(1);
  }
  FORALLSITES(i, s) {
    index = i * len;
    for (a = 0; a < N_K; a++) {
      ops[index] = -1.0 * vevK[a];
      index++;
    }
    for (a = 0; a < N_K; a++) {
      ops[index] = -1.0 * volK[a];
      index++;
    }
    for (a = 0; a < N_K; a++) {
      ops[index] = 0.0;   // Assume vanishing SUGRA vev
      index++;
    }
    for (a = 0; a < N_K; a++) {
      ops[index] = -1.0 * volS[a];
      index++;
    }
  }
  for (j = 0; j < MAX_pts; j++) {
    count[j] = 0;
    for (a = 0; a < 4 * N_K * N_K; a++) {
      CK[j][a] = 0.0;
      CS[j][a] = 0.0;
    }
  }

  // Compute traces of bilinears of scalar field interpolating ops
  compute_Ba();

  // Construct the operators
  FORALLSITES(i, s) {
    for (a = 0; a < NUMLINK; a++) {
      index = i * len;
      // First Konishi subtracting ensemble average
      //          then subtracting   volume average
      for (subs = 0; subs < 2; subs++) {
        for (j = 0; j < N_K; j++) {
          ops[index] += traceBB[j][a][a][i];
          index++;
        }
      }

      // Now SUGRA summed over six components with mu < nu
      // Same story with subtractions
      for (b = 0; b < NUMLINK; b++) {
        index = i * len + 2 * N_K;
        for (subs = 0; subs < 2; subs++) {
          for (j = 0; j < N_K; j++) {
            for (mu = 0; mu < NDIMS; mu++) {
              for (nu = mu + 1; nu < NDIMS; nu++) {
                tr = P[mu][a] * P[nu][b] + P[nu][a] * P[mu][b];
                ops[index] += 0.5 * tr * traceBB[j][a][b][i];
              }
            }
            index++;
          }
        }
      }
    }

    // Finish averaging SUGRA over six components with mu < nu
    index = i * len + 2 * N_K;
    for (j = 0; j < 2 * N_K; j++) {
      ops[index] /= 6.0;
      index++;
    }
  }

  // Construct and optionally print correlators
  // Use general gathers, at least for now
  for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
    disp[XUP] = x_dist;

    // Don't need negative y_dist when x_dist = 0
    if (x_dist > 0)
      y_start = -MAX_X;
    else
      y_start = 0;

    for (y_dist = y_start; y_dist <= MAX_X; y_dist++) {
      disp[YUP] = y_dist;

      // Don't need negative z_dist when both x, y non-positive
      if (x_dist > 0 || y_dist > 0)
        z_start = -MAX_X;
      else
        z_start = 0;

      for (z_dist = z_start; z_dist <= MAX_X; z_dist++) {
        disp[ZUP] = z_dist;

        // Don't need negative t_dist when x, y and z are all non-positive
        if (x_dist > 0 || y_dist > 0 || z_dist > 0)
          t_start = -MAX_X;
        else
          t_start = 0;

        // Ignore any t_dist > MAX_X even if t_dist <= MAX_T
        for (t_dist = t_start; t_dist <= MAX_X; t_dist++) {
          disp[TUP] = t_dist;
          mtag = start_general_gather_field(ops, len * sizeof(Real),
                                            disp, EVENANDODD, gen_pt[0]);

          // Figure out scalar distance while general gather runs
          tr = A4map(x_dist, y_dist, z_dist, t_dist);
          wait_general_gather(mtag);

          if (tr > MAX_r - 1.0e-6) {
            cleanup_general_gather(mtag);
            continue;
          }

          // Combine four-vectors with same scalar distance
          this_r = -1;
          for (j = 0; j < total_r; j++) {
            if (fabs(tr - lookup[j]) < 1.0e-6) {
              this_r = j;
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

          // Don't mix different kinds of subtractions...
          // So only 2 * N_K^2 Konishi correlators
          // !!! Should each N_K x N_K correlator matrix be symmetric???
#ifdef CHECK_ROT
          // Potentially useful to check rotational invariance
          node0_printf("ROT_K %d %d %d %d", x_dist, y_dist, z_dist, t_dist);
#endif
          for (subs = 0; subs < 2; subs++) {
            for (a = 0; a < N_K; a++) {
              for (b = 0; b < N_K; b++) {
                tr = 0.0;
                FORALLSITES(i, s) {
                  index = i * len + subs * N_K + a;
                  tr += ops[index] * ((Real *)gen_pt[0][i])[subs * N_K + b];
                }
                g_doublesum(&tr);
                index = subs * N_K * N_K + a * N_K + b;
                CK[this_r][index] += tr;
#ifdef CHECK_ROT
                // Potentially useful to check rotational invariance
                node0_printf(" %.6g", tr / (Real)volume);
                if (index == N_K * N_K - 1)
                  node0_printf("\n");
#endif
              }
            }
          }

          // 2 * N_K^2 SUGRA correlators
#ifdef CHECK_ROT
          // Potentially useful to check rotational invariance
          node0_printf("ROT_S %d %d %d %d", x_dist, y_dist, z_dist, t_dist);
#endif
          for (subs = 0; subs < 2; subs++) {
            for (a = 2 * N_K; a < 3 * N_K; a++) {
              for (b = 2 * N_K; b < 3 * N_K; b++) {
                tr = 0.0;
                FORALLSITES(i, s) {
                  index = i * len + subs * N_K + a;
                  tr += ops[index] * ((Real *)gen_pt[0][i])[subs * N_K + b];
                }
                g_doublesum(&tr);
                index = subs * N_K * N_K + (a - 2 * N_K) * N_K + (b - 2 * N_K);
                CS[this_r][index] += tr;
#ifdef CHECK_ROT
                // Potentially useful to check rotational invariance
                node0_printf(" %.6g", tr / (Real)volume);
                if (index == N_K * N_K - 1)
                  node0_printf("\n");
#endif
              }
            }
          }
          cleanup_general_gather(mtag);
        } // t_dist
      } // z_dist
    } // y dist
  } // x dist

  // Now cycle through unique scalar distances and print results
  // Won't be sorted, but this is easy to do offline
  for (j = 0; j < total_r; j++) {
    tr = 1.0 / (Real)(count[j] * volume);
    node0_printf("CORR_K %d %.6g", j, lookup[j]);
    for (subs = 0; subs < 2; subs++) {
      for (a = 0; a < N_K; a++) {
        for (b = 0; b < N_K; b++) {
          index = subs * N_K * N_K + a * N_K + b;
          node0_printf(" %.8g", CK[j][index] * tr);
        }
      }
    }
    node0_printf("\n");
  }
  for (j = 0; j < total_r; j++) {
    tr = 1.0 / (Real)(count[j] * volume);
    node0_printf("CORR_S %d %.6g", j, lookup[j]);
    for (subs = 0; subs < 2; subs++) {
      for (a = 0; a < N_K; a++) {
        for (b = 0; b < N_K; b++) {
          index = subs * N_K * N_K + a * N_K + b;
          node0_printf(" %.8g", CS[j][index] * tr);
        }
      }
    }
    node0_printf("\n");
  }
  free(ops);
}
// -----------------------------------------------------------------
