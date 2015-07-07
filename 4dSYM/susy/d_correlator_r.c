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
// vevK[numK] are Konishi vacuum subtractions; assume these vanish for SUGRA
void d_correlator_r() {
  register int i;
  register site *s;
  int a, b, j, index, x_dist, y_dist, z_dist, t_dist;
  int y_start, z_start, t_start;
  int disp[NDIMS] = {0, 0, 0, 0}, this_r, total_r = 0;
  int MAX_pts = 8 * MAX_X * MAX_X * MAX_X, count[MAX_pts];
  Real MAX_r = 100.0 * MAX_X, tr;
  Real lookup[MAX_pts], CK[MAX_pts][numK * numK], CS[MAX_pts][numK * numK];
  Real *OK[numK], *OS[numK];
  msg_tag *Ktag[numK], *Stag[numK];

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

  // Allocate operators and initialize correlators
  for (j = 0; j < numK; j++) {
    OK[j] = malloc(sites_on_node * sizeof(Real));
    OS[j] = malloc(sites_on_node * sizeof(Real));
  }
  for (j = 0; j < MAX_pts; j++) {
    count[j] = 0;
    for (a = 0; a < numK * numK; a++) {
      CK[j][a] = 0.0;
      CS[j][a] = 0.0;
    }
  }

  // Compute traces of bilinears of scalar field interpolating ops
  compute_Ba();

  // Construct the operators
  node0_printf("(Simple SUGRA)\n");
  for (j = 0; j < numK; j++) {
    FORALLSITES(i, s) {
      OK[j][i] = -1.0 * vevK[j];
      OS[j][i] = 0.0;
      for (a = 0; a < NUMLINK; a++) {
        // First Konishi, summing over diagonal less vacuum
        OK[j][i] += traceBB[j][a][a][i];

        // Now SUGRA, averaging over ten off-diagonal a < b
        for (b = a + 1; b < NUMLINK; b++)
          OS[j][i] += traceBB[j][a][b][i];
      }

      // Finish averaging SUGRA over ten components with a < b
      OS[j][i] *= 0.1;
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
          // Check scalar distance before starting general gathers
          tr = A4map(x_dist, y_dist, z_dist, t_dist);
          if (tr > MAX_r - 1.0e-6)
            continue;

          // Konishi general gathers
          disp[TUP] = t_dist;
          Ktag[0] = start_general_gather_field(OK[0], sizeof(Real), disp,
                                               EVENANDODD, gen_pt[0]);
          for (j = 1; j < numK; j++) {
            wait_general_gather(Ktag[j - 1]);
            Ktag[j] = start_general_gather_field(OK[j], sizeof(Real), disp,
                                                 EVENANDODD, gen_pt[j]);
          }

          // Combine four-vectors with same scalar distance
          // Overlap with last Konishi general gather
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

          // SUGRA general gathers
          wait_general_gather(Ktag[numK - 1]);
          Stag[0] = start_general_gather_field(OS[0], sizeof(Real), disp,
                                               EVENANDODD, gen_pt[numK]);
          for (j = 1; j < numK; j++) {
            wait_general_gather(Stag[j - 1]);
            Stag[j] = start_general_gather_field(OS[j], sizeof(Real), disp,
                                                 EVENANDODD, gen_pt[numK + j]);
          }

          // numK^2 Konishi correlators
          // Overlap with last SUGRA general gather
#ifdef CHECK_ROT
          // Potentially useful to check rotational invariance
          node0_printf("ROT_K %d %d %d %d", x_dist, y_dist, z_dist, t_dist);
#endif
          for (a = 0; a < numK; a++) {
            for (b = 0; b < numK; b++) {
              tr = 0.0;
              FORALLSITES(i, s)
                tr += *((Real *)(gen_pt[b][i])) * OK[a][i];
              g_doublesum(&tr);
              CK[this_r][a * numK + b] += tr;
#ifdef CHECK_ROT
              // Potentially useful to check rotational invariance
              node0_printf(" %.6g", tr / (Real)volume);
              if (a == numK - 1 && b == numK - 1)
                node0_printf("\n");
#endif
            }
          }

          // numK^2 SUGRA correlators
          wait_general_gather(Stag[numK - 1]);
#ifdef CHECK_ROT
          // Potentially useful to check rotational invariance
          node0_printf("ROT_S %d %d %d %d", x_dist, y_dist, z_dist, t_dist);
#endif
          for (a = 0; a < numK; a++) {
            for (b = 0; b < numK; b++) {
              tr = 0.0;
              FORALLSITES(i, s)
                tr += *((Real *)(gen_pt[numK + b][i])) * OS[a][i];
              g_doublesum(&tr);
              CS[this_r][a * numK + b] += tr;
#ifdef CHECK_ROT
              // Potentially useful to check rotational invariance
              node0_printf(" %.6g", tr / (Real)volume);
              if (a == numK - 1 && b == numK - 1)
                node0_printf("\n");
#endif
            }
          }
          for (j = 0; j < numK; j++) {
            cleanup_general_gather(Ktag[j]);
            cleanup_general_gather(Stag[j]);
          }
        } // t_dist
      } // z_dist
    } // y dist
  } // x dist

  // Now cycle through unique scalar distances and print results
  // Won't be sorted, but this is easy to do offline
  for (j = 0; j < total_r; j++) {
    tr = 1.0 / (Real)(count[j] * volume);
    node0_printf("CORR_K %d %.6g", j, lookup[j]);
    for (a = 0; a < numK; a++) {
      for (b = 0; b < numK; b++)
        node0_printf(" %.8g", CK[j][a * numK + b] * tr);
    }
    node0_printf("\n");
  }
  for (j = 0; j < total_r; j++) {
    tr = 1.0 / (Real)(count[j] * volume);
    node0_printf("CORR_S %d %.6g", j, lookup[j]);
    for (a = 0; a < numK; a++) {
      for (b = 0; b < numK; b++)
        node0_printf(" %.8g", CS[j][a * numK + b] * tr);
    }
    node0_printf("\n");
  }
  for (j = 0; j < numK; j++) {
    free(OK[j]);
    free(OS[j]);
  }
}
// -----------------------------------------------------------------
