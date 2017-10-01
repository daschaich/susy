// -----------------------------------------------------------------
// Helper functions for the correlators
#include "corr_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Map (x, y, z, t) to r on Nx x Ny x Nz x Nt A4* lattice
// Check all possible periodic shifts to find true r
Real A4map(x_in, y_in, z_in, t_in) {
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
int count_points(Real MAX_r) {
  int total = 0, this, j, x_dist, y_dist, z_dist, t_dist;
  int y_start, z_start, t_start;
  int *count = malloc(MAX_pts * sizeof(*count));
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
            if (fabs(tr - lookup[j]) < 1.0e-6) {
              this = j;
              count[this]++;
              break;
            }
          }
          // If this scalar distance is new, add it to the list
          if (this < 0) {
            lookup[total] = tr;
            this = total;
            total++;
            count[this] = 1;    // Initialize
            if (total > MAX_pts) {
              node0_printf("count_points: MAX_pts %d too small\n", MAX_pts);
              terminate(1);
            }
          }
        }
      }
    }
  }

  // Convert count into normalization
  for (j = 0; j < total; j++)
    norm[j] = 1.0 / (Real)(count[j] * volume);
  free(count);

  return total;
}
// -----------------------------------------------------------------
