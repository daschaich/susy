// -----------------------------------------------------------------
// Static potential for all displacements
// Evaluate in different spatial dirs to check rotational invariance
// Must gauge fix to Coulomb gauge before calling
// This version computes spatial correlators of temporal products
#include "susy_includes.h"

// Define CHECK_ROT to check rotational invariance
//#define CHECK_ROT
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Map (x, y, z) to r on Nx x Ny x Nz A4* time slice
// Check all possible periodic shifts to find true (smallest) r
Real A4map_slice(x_in, y_in, z_in) {
  int x, y, z, xSq, ySq, xy, xpy;
  Real r = 100.0 * MAX_X, tr;
#ifdef DEBUG_CHECK
  Real invSq2 = 1.0 / sqrt(2.0), invSq6  = 1.0 / sqrt(6.0);
  Real invSq12 = 1.0 / sqrt(12.0), x_a4, y_a4, z_a4, check;
#endif

  for (x = x_in - nx; x <= x_in + nx; x += nx) {
    xSq = x * x;
    for (y = y_in - ny; y <= y_in + ny; y += ny) {
      ySq = y * y;
      xy = x * y;
      xpy = x + y;
      for (z = z_in - nz; z <= z_in + nz; z += nz) {
        tr = sqrt((xSq + ySq + z * z) * 0.75 - (xy + xpy * z) * 0.5);
        if (tr < r)
          r = tr;

#ifdef DEBUG_CHECK
        x_a4 = (x - y) * invSq2;
        y_a4 = (x + y - 2.0 * z) * invSq6;
        z_a4 = (x + y + z) * invSq12;
        check = sqrt(x_a4 * x_a4 + y_a4 * y_a4 + z_a4 * z_a4);
        if (fabs(tr - check) > IMAG_TOL) {
          node0_printf("ERROR: %4g isn't %.4g for (%d, %d, %d)\n",
                       tr, check, x, y, z);
          terminate(1);
        }
#endif
      }
    }
  }
  return r;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use tempmat, tempmat2 and staple for temporary storage
void hvy_pot(int do_det) {
  register int i;
  register site *s;
  int j, t_dist, x_dist, y_dist, z_dist, y_start, z_start;
  int MAX_pts = 8 * MAX_X * MAX_X;            // Should be plenty
  int count[MAX_pts], this_r, total_r = 0;
  Real MAX_r = 100.0 * MAX_X, tr, lookup[MAX_pts], W[MAX_pts];
  double wloop;
  complex tc;
  matrix_f tmat, tmat2, *mat;
  msg_tag *mtag = NULL;

  // Initialize Wilson loop accumulators
  for (j = 0; j < MAX_pts; j++) {
    count[j] = 0;
    W[j] = 0.0;
  }

  // Find smallest scalar distance cut by imposing MAX_X
  for (y_dist = 0; y_dist <= MAX_X; y_dist++) {
    for (z_dist = 0; z_dist <= MAX_X; z_dist++) {
      tr = A4map_slice(MAX_X + 1, y_dist, z_dist);
      if (tr < MAX_r)
        MAX_r = tr;
    }
  }
  node0_printf("hvy_pot: MAX_T = %d, MAX_X = %d --> r < %.6g\n",
               MAX_T, MAX_X, MAX_r);

  // Use staple to hold product of t_dist links at each (x, y, z)
  for (t_dist = 1; t_dist <= MAX_T; t_dist++) {
    if (t_dist == 1) {
      FORALLSITES(i, s)
        mat_copy_f(&(s->linkf[TUP]), &(staple[i]));
    }
    else {
      mtag = start_gather_field(staple, sizeof(matrix_f),
                                goffset[TUP], EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s) {
        mat = (matrix_f *)gen_pt[0][i];
        mult_nn_f(&(s->linkf[TUP]), mat, &(tempmat2[i]));
      }
      cleanup_gather(mtag);
      FORALLSITES(i, s)
        mat_copy_f(&(tempmat2[i]), &(staple[i]));
    }

    for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
      if (x_dist > 0)
        y_start = -MAX_X;
      else
        y_start = 0;    // Don't need negative y_dist when x_dist = 0
      for (y_dist = y_start; y_dist <= MAX_X; y_dist++) {
        // Gather staple to tempmat along spatial offset, using tempmat2
        FORALLSITES(i, s)
          mat_copy_f(&(staple[i]), &(tempmat[i]));
        for (j = 0; j < x_dist; j++)
          shiftmat(tempmat, tempmat2, goffset[XUP]);
        for (j = 0; j < y_dist; j++)
          shiftmat(tempmat, tempmat2, goffset[YUP]);
        for (j = y_dist; j < 0; j++)
          shiftmat(tempmat, tempmat2, goffset[YUP] + 1);

        // If either x_dist or y_dist are positive,
        // we need to start with MAX_X shifts in the -z direction
        if (x_dist > 0 || y_dist > 0)
          z_start = -MAX_X;
        else
          z_start = 0;
        for (j = z_start; j < 0; j++)
          shiftmat(tempmat, tempmat2, goffset[ZUP] + 1);
        for (z_dist = z_start; z_dist <= MAX_X; z_dist++) {
          // Figure out scalar distance
          tr = A4map_slice(x_dist, y_dist, z_dist);
          if (tr > MAX_r - 1.0e-6) {
#ifdef CHECK_ROT
#ifdef DEBUG_CHECK
            // Potentially useful to check against old output
            if (do_det == 1) {  // Braces fix compiler error
              node0_printf("D_LOOP   ");
            }
            else
              node0_printf("POT_LOOP ");
            node0_printf("%d %d %d %d SKIP\n", x_dist, y_dist, z_dist, t_dist);
#endif
#endif
            // Need to increment z and shift in z direction
            shiftmat(tempmat, tempmat2, goffset[ZUP]);
            continue;
          }

          // Combine three-vectors with same scalar distance
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
              node0_printf("hvy_pot: MAX_pts %d too small\n", MAX_pts);
              terminate(1);
            }
          }
          count[this_r]++;

          // Evaluate potential at this separation
          wloop = 0.0;
          FORALLSITES(i, s) {
            // Compute the actual Coulomb gauge Wilson loop product
            mult_na_f(&(staple[i]), &(tempmat[i]), &tmat);

            if (do_det == 1)
              det_project(&tmat, &tmat2);
            else
              mat_copy_f(&tmat, &tmat2);

            tc = trace_f(&tmat2);
            wloop += tc.real;
          }
          g_doublesum(&wloop);
          W[this_r] += wloop;
#ifdef CHECK_ROT
          // Potentially useful to check rotational invariance
          if (do_det == 1) {  // Braces fix compiler error
            node0_printf("D_LOOP   ");
          }
          else
            node0_printf("POT_LOOP ");
          node0_printf("%d %d %d %d %.6g\n", x_dist, y_dist, z_dist, t_dist,
                       wloop / volume);
#endif

          // As we increment z, shift in z direction
          shiftmat(tempmat, tempmat2, goffset[ZUP]);
        } // z_dist
      } // y_dist
    } // x_dist

    // Now cycle through unique scalar distances and print results
    // Won't be sorted, but this is easy to do offline
    // Also reset for next t_dist
    for (j = 0; j < total_r; j++) {
      tr = 1.0 / (Real)(count[j] * volume);
      if (do_det == 1) {            // Braces suppress compiler complaint
        node0_printf("D_LOOP   ");
      }
      else
        node0_printf("POT_LOOP ");
      node0_printf("%d %.6g %d %.6g\n", j, lookup[j], t_dist, W[j] * tr);
      count[j] = 0;
      W[j] = 0.0;
    }
  } // t_dist
}
// -----------------------------------------------------------------
