// -----------------------------------------------------------------
// Measure the Konishi and SUGRA correlation functions
// Combine Konishi and SUGRA into single structs
// Then only need one gather per shift
#include "susy_includes.h"

// Define CHECK_ROT to check rotational invariance
//#define CHECK_ROT
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Map (x, y, z, t) to r on Nx x Ny x Nz x Nt A4* lattice
// Check all possible periodic shifts to find true r
Real A4map(int x_in, int y_in, int z_in, int t_in) {
  int x, y, z, t, xSq, ySq, zSq, xy, xpy, xpyz, xpypz;
  Real r = 100.0 * MAX_X, tr;   // r to be overwritten

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
// Count and return total number of unique scalar distances
int count_points(Real MAX_r, Real *save, int *count, int MAX_pts) {
  int total = 0, this, j, x_dist, y_dist, z_dist, t_dist;
  int y_start, z_start, t_start;
  Real tr;

  for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
    // Don't need negative y_dist when x_dist = 0
    if (x_dist > 0)
      y_start = -MAX_X;
    else
      y_start = 0;

    for (y_dist = y_start; y_dist <= MAX_X; y_dist++) {
      // Don't need negative z_dist when both x, y non-positive
      if (x_dist > 0 || y_dist > 0)
        z_start = -MAX_X;
      else
        z_start = 0;

      for (z_dist = z_start; z_dist <= MAX_X; z_dist++) {
        // Don't need negative t_dist when x, y and z are all non-positive
        if (x_dist > 0 || y_dist > 0 || z_dist > 0)
          t_start = -MAX_X;
        else
          t_start = 0;

        // Ignore any t_dist > MAX_X even if t_dist <= MAX_T
        for (t_dist = t_start; t_dist <= MAX_X; t_dist++) {
          // Figure out the scalar distance, see if it's smaller than MAX_r
          tr = A4map(x_dist, y_dist, z_dist, t_dist);
          if (tr > MAX_r - 1.0e-6)
            continue;

          // See if this scalar distance is already accounted for
          this = -1;
          for (j = 0; j < total; j++) {
            if (fabs(tr - save[j]) < 1.0e-6) {
              this = j;
              count[this]++;
              break;
            }
          }
          // If this scalar distance is new, add it to the list
          if (this < 0) {
            save[total] = tr;
            this = total;
            total++;
            count[this] = 1;    // Initialize
            if (total >= MAX_pts) {
              node0_printf("count_points: MAX_pts %d too small\n", MAX_pts);
              terminate(1);
            }
          }
        }
      }
    }
  }
  return total;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Simple helper function to copy operators
void copy_ops(Kops *src, Kops *dest) {
  int j;
  for (j = 0; j < N_K; j++) {
    dest->OK[j] = src->OK[j];
    dest->OS[j] = src->OS[j];
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Shift all operators without parallel transport
// The dir should come from goffset
void shift_ops(Kops *dat, Kops *temp, int dir) {
  register int i;
  register site *s;
  msg_tag *mtag;

  mtag = start_gather_field(dat, sizeof(Kops), dir, EVENANDODD, gen_pt[0]);
  wait_gather(mtag);
  FORALLSITES(i, s)
    copy_ops((Kops *)gen_pt[0][i], &(temp[i]));
  cleanup_gather(mtag);
  FORALLSITES(i, s)
    copy_ops(&(temp[i]), &(dat[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA correlators as functions of r
// For gathering, all operators live in Kops structs
void correlator_r() {
  register int i;
  register site *s;
  int a, b, j, x_dist, y_dist, z_dist, t_dist;
  int y_start, z_start, t_start, this_r, total_r = 0;
  int *count, *temp, MAX_pts = 8 * MAX_X * MAX_X * MAX_X;
  Real MAX_r = 100.0 * MAX_X, tr;   // MAX_r to be overwritten
  Real *lookup, *sav;
  Kops *ops = malloc(sites_on_node * sizeof(*ops));
  Kcorrs *CK, *CS;

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
  node0_printf("correlator_r: MAX = %d --> r < %.6g\n", MAX_X, MAX_r);

  // Count number of unique scalar distances r < MAX_r
  // Copy sav and temp into lookup and count, respectively
  sav = malloc(MAX_pts * sizeof(*sav));
  temp = malloc(MAX_pts * sizeof(*temp));
  total_r = count_points(MAX_r, sav, temp, MAX_pts);
  lookup = malloc(total_r * sizeof(*lookup));
  count = malloc(total_r * sizeof(*count));
  for (j = 0; j < total_r; j++) {
    lookup[j] = sav[j];
    count[j] = temp[j];
  }
  free(sav);
  free(temp);

  // Allocate correlator arrays now that we know how big they should be
  CK = malloc(total_r * sizeof(*CK));
  CS = malloc(total_r * sizeof(*CS));

  // Initialize operators and correlators
  FORALLSITES(i, s) {
    for (j = 0; j < N_K; j++) {
      ops[i].OK[j] = 0.0;
      ops[i].OS[j] = 0.0;
    }
  }
  for (j = 0; j < total_r; j++) {
    for (a = 0; a < N_K; a++) {
      for (b = 0; b < N_K; b++) {
        CK[j].C[a][b] = 0.0;
        CS[j].C[a][b] = 0.0;
      }
    }
  }

  // Compute traces of bilinears of scalar field interpolating ops
  compute_Ba();

  // Construct the operators for all definitions
  FORALLSITES(i, s) {
    FORALLDIR(a) {
      for (j = 0; j < N_K; j++) {
        // First Konishi
        ops[i].OK[j] += traceBB[j][a][a][i];

        // Now SUGRA, averaged over 20 off-diagonal components
        FORALLDIR(b) {
          if (a == b)
            continue;
          ops[i].OS[j] += 0.05 * traceBB[j][a][b][i];
        }
      }
    }
  }

  // Construct and optionally print correlators
  for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
    // Don't need negative y_dist when x_dist = 0
    if (x_dist > 0)
      y_start = -MAX_X;
    else
      y_start = 0;

    for (y_dist = y_start; y_dist <= MAX_X; y_dist++) {
      // Don't need negative z_dist when both x, y non-positive
      if (x_dist > 0 || y_dist > 0)
        z_start = -MAX_X;
      else
        z_start = 0;

      for (z_dist = z_start; z_dist <= MAX_X; z_dist++) {
        // Gather ops to tempops along spatial offset, using tempops2
        FORALLSITES(i, s)
          copy_ops(&(ops[i]), &(tempops[i]));
        for (j = 0; j < x_dist; j++)
          shift_ops(tempops, tempops2, goffset[XUP]);
        for (j = 0; j < y_dist; j++)
          shift_ops(tempops, tempops2, goffset[YUP]);
        for (j = y_dist; j < 0; j++)
          shift_ops(tempops, tempops2, goffset[YUP] + 1);
        for (j = 0; j < z_dist; j++)
          shift_ops(tempops, tempops2, goffset[ZUP]);
        for (j = z_dist; j < 0; j++)
          shift_ops(tempops, tempops2, goffset[ZUP] + 1);

        // Don't need negative t_dist when x, y and z are all non-positive
        // Otherwise we need to start with MAX_X shifts in the -t direction
        if (x_dist > 0 || y_dist > 0 || z_dist > 0)
          t_start = -MAX_X;
        else
          t_start = 0;
        for (j = t_start; j < 0; j++)
          shift_ops(tempops, tempops2, goffset[TUP] + 1);

        // Ignore any t_dist > MAX_X even if t_dist <= MAX_T
        for (t_dist = t_start; t_dist <= MAX_X; t_dist++) {
          // Figure out scalar distance
          tr = A4map(x_dist, y_dist, z_dist, t_dist);
          if (tr > MAX_r - 1.0e-6) {
            // Only increment t, but still need to shift in t direction
            shift_ops(tempops, tempops2, goffset[TUP]);
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
          if (this_r < 0) {
            node0_printf("correlator_r: bad scalar distance %.4g ", tr);
            node0_printf("from displacement %d %d %d %d\n",
                         x_dist, y_dist, z_dist, t_dist);
            terminate(1);
          }

          // N_K^2 Konishi correlators
#ifdef CHECK_ROT
          // Potentially useful to check rotational invariance
          node0_printf("ROT_K %d %d %d %d", x_dist, y_dist, z_dist, t_dist);
#endif
          for (a = 0; a < N_K; a++) {
            for (b = 0; b < N_K; b++) {
              tr = 0.0;
              FORALLSITES(i, s)
                tr += ops[i].OK[a] * tempops[i].OK[b];
              g_doublesum(&tr);
              CK[this_r].C[a][b] += tr;
#ifdef CHECK_ROT
              node0_printf(" %.6g", tr / (Real)volume);
#endif
            }
          }

          // N_K^2 SUGRA correlators
#ifdef CHECK_ROT
          // Potentially useful to check rotational invariance
          node0_printf("\nROT_S %d %d %d %d", x_dist, y_dist, z_dist, t_dist);
#endif
          for (a = 0; a < N_K; a++) {
            for (b = 0; b < N_K; b++) {
              tr = 0.0;
              FORALLSITES(i, s)
                tr += ops[i].OS[a] * tempops[i].OS[b];
              g_doublesum(&tr);
              CS[this_r].C[a][b] += tr;
#ifdef CHECK_ROT
              node0_printf(" %.6g", tr / (Real)volume);
#endif
            }
          }
#ifdef CHECK_ROT
          node0_printf("\n");
#endif
          // As we increment t, shift in t direction
          shift_ops(tempops, tempops2, goffset[TUP]);
        } // t_dist
      } // z_dist
    } // y dist
  } // x dist
  free(ops);

  // Now cycle through unique scalar distances and print results
  // Distances won't be sorted, but this is easy to do offline
  // Format: CORR_?  tag  r  a  b  corr[a][b]
  for (j = 0; j < total_r; j++) {
    tr = 1.0 / (Real)(count[j] * volume);
    for (a = 0; a < N_K; a++) {
      for (b = 0; b < N_K; b++) {
        node0_printf("CORR_K %d %.6g %d %d %.8g\n", j, lookup[j], a, b,
                     CK[j].C[a][b] * tr);
      }
    }
  }
  for (j = 0; j < total_r; j++) {
    tr = 1.0 / (Real)(count[j] * volume);
    for (a = 0; a < N_K; a++) {
      for (b = 0; b < N_K; b++) {
        node0_printf("CORR_S %d %.6g %d %d %.8g\n", j, lookup[j], a, b,
                     CS[j].C[a][b] * tr);
      }
    }
  }
  free(lookup);
  free(count);
  free(CK);
  free(CS);
}
// -----------------------------------------------------------------
