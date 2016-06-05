// -----------------------------------------------------------------
// Measure the Konishi and SUGRA correlation functions
// Use general gathers, but combine Konishi and SUGRA into single vector
// Then only need one general gather per displacement
#include "susy_includes.h"

// Define CHECK_ROT to check rotational invariance
//#define CHECK_ROT
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Map (x, t) to r on Nx x Nt square lattice
// Check all possible periodic shifts to find true r
Real r_map(x_in, t_in) {
  int x, t, xSq;
  Real r = 100.0 * MAX_X, tr;   // r to be overwritten

  for (x = x_in - nx; x <= x_in + nx; x += nx) {
    xSq = x * x;
    for (t = t_in - nt; t <= t_in + nt; t += nt) {
      tr = sqrt((xSq + t * t));
      if (tr < r)
        r = tr;
    }
  }
  return r;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Count and return total number of unique scalar distances
int count_points(Real MAX_r, Real *save, int *count, int MAX_pts) {
  int total = 0, this, j, x_dist, t_dist, t_start;
  Real tr;

  for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
    // Don't need negative t_dist when x, y and z are all non-positive
    if (x_dist > 0)
      t_start = -MAX_X;
    else
      t_start = 0;

    // Ignore any t_dist > MAX_X even if t_dist <= MAX_T
    for (t_dist = t_start; t_dist <= MAX_X; t_dist++) {
      // Figure out the scalar distance, see if it's smaller than MAX_r
      tr = r_map(x_dist, t_dist);
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
  return total;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Simple helper function to copy operators
void copy_ops(Kops *src, Kops *dest) {
  int j, sub;
  for (j = 0; j < N_K; j++) {
    for (sub = 0; sub < 2; sub++) {
      dest->OK[j][sub] = src->OK[j][sub];
      dest->OS[j][sub] = src->OS[j][sub];
    }
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
// For gathering, all operators live in array of len = 4N_K
// vevK[N_K] and vevS[N_K] are Konishi and SUGRA ensemble averages
// volK[N_K] and volS[N_K] are Konishi and SUGRA volume averages
// Check to make sure the latter have been set by konishi()
void correlator_r() {
  register int i;
  register site *s;
  int a, b, j, x_dist, t_dist, t_start, sub, this_r, total_r = 0;
  int *count, *temp, MAX_pts = 8 * MAX_X;
  Real MAX_r = 100.0 * MAX_X, tr;   // MAX_r to be overwritten
  Real *lookup, *sav;
  Kops *ops = malloc(sites_on_node * sizeof(*ops));
  Kcorrs *CK, *CS;

  // Find smallest scalar distance cut by imposing MAX_X
  for (t_dist = 0; t_dist <= MAX_X; t_dist++) {
    tr = r_map(MAX_X + 1, t_dist);
    if (tr < MAX_r)
      MAX_r = tr;
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

  // Subtract either ensemble average or volume average (respectively)
  // while initializing operators and correlators
  // Make sure vevK and vevS have been set up properly
  if (volK[0] < 0.0) {
    node0_printf("ERROR: Improperly initialized vacuum subtraction\n");
    terminate(1);
  }
  FORALLSITES(i, s) {
    for (j = 0; j < N_K; j++) {
      ops[i].OK[j][0] = -1.0 * vevK[j];
      ops[i].OK[j][1] = -1.0 * volK[j];
      ops[i].OS[j][0] = -1.0 * vevS[j];
      ops[i].OS[j][1] = -1.0 * volS[j];
    }
  }
  for (j = 0; j < total_r; j++) {
    for (a = 0; a < N_K; a++) {
      for (b = 0; b < N_K; b++) {
        for (sub = 0; sub < 2; sub++) {
          CK[j].C[a][b][sub] = 0.0;
          CS[j].C[a][b][sub] = 0.0;
        }
      }
    }
  }

  // Compute traces of bilinears of scalar field interpolating ops
  compute_Ba();

  // Construct the operators for all definitions and subtractions
  FORALLSITES(i, s) {
    FORALLDIR(a) {
      for (j = 0; j < N_K; j++) {
        for (sub = 0; sub < 2; sub++) {
          // First Konishi
          ops[i].OK[j][sub] += traceBB[j][a][a][i];

          FORALLDIR(b) {
            // Add Tr phi_6^2 = (1 / 5) Tr \sum_{a, b} A_a A_b
            // to cancel SUGRA mixing in log-polar Konishi
            if (j == 0)
              ops[i].OK[j][sub] += 0.2 * traceAA[a][b][i];

            // Now SUGRA, averaged over 20 off-diagonal components
            if (a == b)
              continue;
            ops[i].OS[j][sub] += 0.05 * traceBB[j][a][b][i];
          }
        }
      }
    }
  }

  // Construct and optionally print correlators
  // Use general gathers, at least for now
  for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
    // Gather ops to tempops along spatial offset, using tempops2
    FORALLSITES(i, s)
      copy_ops(&(ops[i]), &(tempops[i]));
    for (j = 0; j < x_dist; j++)
      shift_ops(tempops, tempops2, goffset[XUP]);

    // Don't need negative t_dist when x_dist = 0
    // Otherwise we need to start with MAX_X shifts in the -t direction
    if (x_dist > 0)
      t_start = -MAX_X;
    else
      t_start = 0;
    for (j = t_start; j < 0; j++)
      shift_ops(tempops, tempops2, goffset[TUP] + 1);

    // Ignore any t_dist > MAX_X even if t_dist <= MAX_T
    for (t_dist = t_start; t_dist <= MAX_X; t_dist++) {
      // Figure out scalar distance
      tr = r_map(x_dist, t_dist);
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
        node0_printf("from displacement %d %d\n", x_dist, t_dist);
        terminate(1);
      }

      // 2 * N_K^2 Konishi correlators
#ifdef CHECK_ROT
      // Potentially useful to check rotational invariance
      node0_printf("ROT_K %d %d", x_dist, t_dist);
#endif
      for (sub = 0; sub < 2; sub++) {
        for (a = 0; a < N_K; a++) {
          for (b = 0; b < N_K; b++) {
            tr = 0.0;
            FORALLSITES(i, s)
              tr += ops[i].OK[a][sub] * tempops[i].OK[b][sub];
            g_doublesum(&tr);
            CK[this_r].C[a][b][sub] += tr;
#ifdef CHECK_ROT
            node0_printf(" %.6g", tr / (Real)volume);
#endif
          }
        }
      }

      // 2 * N_K^2 SUGRA correlators
#ifdef CHECK_ROT
      // Potentially useful to check rotational invariance
      node0_printf("\nROT_S %d %d", x_dist, t_dist);
#endif
      for (sub = 0; sub < 2; sub++) {
        for (a = 0; a < N_K; a++) {
          for (b = 0; b < N_K; b++) {
            tr = 0.0;
            FORALLSITES(i, s)
              tr += ops[i].OS[a][sub] * tempops[i].OS[b][sub];
            g_doublesum(&tr);
            CS[this_r].C[a][b][sub] += tr;
#ifdef CHECK_ROT
            node0_printf(" %.6g", tr / (Real)volume);
#endif
          }
        }
      }
#ifdef CHECK_ROT
      node0_printf("\n");
#endif
      // As we increment t, shift in t direction
      shift_ops(tempops, tempops2, goffset[TUP]);
    } // t_dist
  } // x dist
  free(ops);

  // Now cycle through unique scalar distances and print results
  // Distances won't be sorted, but this is easy to do offline
  // Format: CORR_?  tag  r  a  b  corr[a][b]  subtracted[a][b]
  for (j = 0; j < total_r; j++) {
    tr = 1.0 / (Real)(count[j] * volume);
    for (a = 0; a < N_K; a++) {
      for (b = 0; b < N_K; b++) {
        node0_printf("CORR_K %d %.6g %d %d %.8g %.8g\n", j, lookup[j], a, b,
                     CK[j].C[a][b][0] * tr, CK[j].C[a][b][1] * tr);
      }
    }
  }
  for (j = 0; j < total_r; j++) {
    tr = 1.0 / (Real)(count[j] * volume);
    for (a = 0; a < N_K; a++) {
      for (b = 0; b < N_K; b++) {
        node0_printf("CORR_S %d %.6g %d %d %.8g %.8g\n", j, lookup[j], a, b,
                     CS[j].C[a][b][0] * tr, CS[j].C[a][b][1] * tr);
      }
    }
  }
  free(lookup);
  free(count);
  free(CK);
  free(CS);
}
// -----------------------------------------------------------------
