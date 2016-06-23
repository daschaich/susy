// -----------------------------------------------------------------
// Static potential for all displacements
// Print all permutations to enable tree-level improvement
// Must gauge fix to Coulomb gauge before calling
// Gauge-fixed links unitarized via polar projection -- overwritten!!!
// This version computes spatial correlators of temporal products
// Use tempmat, tempmat2 and staple for temporary storage
#include "susy_includes.h"

void hvy_pot_polar() {
  register int i;
  register site *s;
  int j, t_dist, x_dist, y_dist, z_dist, y_start, z_start;
  double wloop;
  matrix tmat, tmat2, *mat;
  msg_tag *mtag = NULL;

  node0_printf("hvy_pot_polar: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  FORALLSITES(i, s) {
   // Polar projection of gauge-fixed links
   // To be multiplied together after projecting
   // !!! Overwrites links
   polar(&(s->link[TUP]), &tmat, &tmat2);
   mat_copy(&tmat, &(s->link[TUP]));
  }

  // Use staple to hold product of t_dist links at each (x, y, z)
  for (t_dist = 1; t_dist <= MAX_T; t_dist++) {
    if (t_dist == 1) {
      FORALLSITES(i, s)
        mat_copy(&(s->link[TUP]), &(staple[i]));
    }
    else {
      mtag = start_gather_field(staple, sizeof(matrix),
                                goffset[TUP], EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s) {
        mat = (matrix *)gen_pt[0][i];
        mult_nn(&(s->link[TUP]), mat, &(tempmat2[i]));
      }
      cleanup_gather(mtag);
      FORALLSITES(i, s)
        mat_copy(&(tempmat2[i]), &(staple[i]));
    }

    for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
      if (x_dist > 0)
        y_start = -MAX_X;
      else
        y_start = 0;    // Don't need negative y_dist when x_dist = 0
      for (y_dist = y_start; y_dist <= MAX_X; y_dist++) {
        // Gather staple to tempmat along spatial offset, using tempmat2
        FORALLSITES(i, s)
          mat_copy(&(staple[i]), &(tempmat[i]));
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
          // Evaluate potential at this separation
          // Print even if scalar distance might wrap around lattice
          // Tree-level improved result might fit...
          // Need to check in offline analyses
          wloop = 0.0;
          FORALLSITES(i, s)
            wloop += (double)realtrace(&(staple[i]), &(tempmat[i]));
          g_doublesum(&wloop);

          node0_printf("POLAR_LOOP %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist, wloop / volume);

          // As we increment z, shift in z direction
          shiftmat(tempmat, tempmat2, goffset[ZUP]);
        } // z_dist
      } // y_dist
    } // x_dist
  } // t_dist
}
// -----------------------------------------------------------------
