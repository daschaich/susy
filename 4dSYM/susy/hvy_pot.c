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
// Check all possible periodic shifts to find true r
Real A4map_slice(x_in, y_in, z_in) {
  int x, y, z, xSq, ySq, xy, xpy;
  Real r = 100.0 * MAX_X, tr;

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
      }
    }
  }
  return r;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use tempmat1, tempmat2 and staple for temporary storage
void hvy_pot(int do_det) {
  register int i;
  register site *s;
  int t_dist, x_dist, y_dist, z_dist, y_start, z_start;
  int MAX_pts = 8 * MAX_X * MAX_X;            // Should be plenty
  int count[MAX_pts], this_r, total_r = 0;
  Real MAX_r = 100.0 * MAX_X, tr;
  Real lookup[MAX_pts], W[MAX_pts];
  double wloop;
  complex tc;
  su3_matrix_f tmat, tmat2, *mat;
  msg_tag *mtag = NULL;

  // Initialize Wilson loop accumulators
  for (i = 0; i < MAX_pts; i++) {
    count[i] = 0;
    W[i] = 0.0;
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
        su3mat_copy_f(&(s->linkf[TUP]), &(staple[i]));
    }
    else {
      mtag = start_gather_field(staple, sizeof(su3_matrix_f),
                                goffset[TUP], EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s) {
        mat = (su3_matrix_f *)gen_pt[0][i];
        mult_su3_nn_f(&(s->linkf[TUP]), mat, &(tempmat2[i]));
      }
      cleanup_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy_f(&(tempmat2[i]), &(staple[i]));
    }

    for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
      if (x_dist > 0)
        y_start = -MAX_X;
      else
        y_start = 0;    // Don't need negative y_dist when x_dist = 0
      for (y_dist = y_start; y_dist <= MAX_X; y_dist++) {
        // Gather staple to tempmat1 along spatial offset, using tempmat2
        FORALLSITES(i, s)
          su3mat_copy_f(&(staple[i]), &(tempmat1[i]));
        for (i = 0; i < x_dist; i++) {
          shiftmat(tempmat1, tempmat2, goffset[XUP]);
          FORALLSITES(i, s)
            su3mat_copy_f(&(tempmat2[i]), &(tempmat1[i]));
        }
        if (y_dist > 0) {
          for (i = 0; i < y_dist; i++) {
            shiftmat(tempmat1, tempmat2, goffset[YUP]);
            FORALLSITES(i, s)
              su3mat_copy_f(&(tempmat2[i]), &(tempmat1[i]));
          }
        }
        else if (y_dist < 0) {
          for (i = y_dist; i < 0; i++) {
            shiftmat(tempmat1, tempmat2, goffset[YUP] + 1);
            FORALLSITES(i, s)
              su3mat_copy_f(&(tempmat2[i]), &(tempmat1[i]));
          }
        }

        // If either x_dist or y_dist are positive,
        // we need to start with MAX_X shifts in the -z direction
        if (x_dist > 0 || y_dist > 0)
          z_start = -MAX_X;
        else
          z_start = 0;
        for (i = z_start; i < 0; i++) {
          shiftmat(tempmat1, tempmat2, goffset[ZUP] + 1);
          FORALLSITES(i, s)
            su3mat_copy_f(&(tempmat2[i]), &(tempmat1[i]));
        }
        for (z_dist = z_start; z_dist <= MAX_X; z_dist++) {
          // Figure out scalar distance
          tr = A4map_slice(x_dist, y_dist, z_dist);
          if (tr > MAX_r - 1.0e-6) {
            // Need to move on to next t_dist
            shiftmat(tempmat1, tempmat2, goffset[ZUP]);
            FORALLSITES(i, s)
              su3mat_copy_f(&(tempmat2[i]), &(tempmat1[i]));
            continue;
          }

          // Combine three-vectors with same scalar distance
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
              node0_printf("hvy_pot: MAX_pts %d too small\n", MAX_pts);
              terminate(1);
            }
          }
          count[this_r]++;

          // Evaluate potential at this separation
          wloop = 0.0;
          FORALLSITES(i, s) {
            // Compute the actual Coulomb gauge Wilson loop product
            mult_su3_na_f(&(staple[i]), &(tempmat1[i]), &tmat);

            if (do_det == 1)
              det_project(&tmat, &tmat2);
            else
              su3mat_copy_f(&tmat, &tmat2);

            tc = trace_su3_f(&tmat2);
            wloop += tc.real;
          }
          g_doublesum(&wloop);
          W[this_r] += wloop;
#ifdef CHECK_ROT
          // Potentially useful to check rotational invariance
          if (do_det == 1)
            node0_printf("D_LOOP   ");
          else
            node0_printf("POT_LOOP ");
          node0_printf("%d %d %d %d %.6g\n", x_dist, y_dist, z_dist, t_dist,
                       wloop / volume);
#endif

          // As we increment z, shift in z direction
          shiftmat(tempmat1, tempmat2, goffset[ZUP]);
          FORALLSITES(i, s)
            su3mat_copy_f(&(tempmat2[i]), &(tempmat1[i]));
        } // z_dist
      } // y_dist
    } // x_dist

    // Now cycle through unique scalar distances and print results
    // Won't be sorted, but this is easy to do offline
    // Also reset for next t_dist
    for (i = 0; i < total_r; i++) {
      tr = 1.0 / (Real)(count[i] * volume);
      if (do_det == 1) {            // Braces suppress compiler complaint
        node0_printf("D_LOOP   ");
      }
      else
        node0_printf("POT_LOOP ");
      node0_printf("%d %.6g %d %.6g\n", i, lookup[i], t_dist, W[i] * tr);
      count[i] = 0;
      W[i] = 0.0;
    }
  } // t_dist
}
// -----------------------------------------------------------------
