// -----------------------------------------------------------------
// Static potential for all displacements
// Must gauge fix to Coulomb gauge before calling
// Gauge-fixed links unitarized via polar projection -- overwritten!!!
// This version computes spatial correlators of temporal products
// Use tempmat, tempmat2 and staple for temporary storage
#include "susy_includes.h"

void hvy_pot_polar() {
  register int i;
  register site *s;
  int j, t_dist, x_dist;
  double wloop;
  matrix tmat, tmat2;
  msg_tag *mtag = NULL;

  node0_printf("hvy_pot_polar: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  FORALLSITES(i, s) {
   // Polar projection of gauge-fixed links
   // To be multiplied together after projecting
   // !!! Overwrites links
   polar(&(s->link[TUP]), &tmat, &tmat2);
   mat_copy(&tmat, &(s->link[TUP]));
  }

  // Use staple to hold product of t_dist links at each point
  for (t_dist = 1; t_dist <= MAX_T; t_dist++) {
    if (t_dist == 1) {
      FORALLSITES(i, s)
        mat_copy(&(s->link[TUP]), &(staple[i]));
    }
    else {
      mtag = start_gather_field(staple, sizeof(matrix),
                                goffset[TUP], EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s)
        mult_nn(&(s->link[TUP]), (matrix *)gen_pt[0][i], &(tempmat2[i]));
      cleanup_gather(mtag);
      FORALLSITES(i, s)
        mat_copy(&(tempmat2[i]), &(staple[i]));
    }

    for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
      // Gather staple to tempmat along spatial offset, using tempmat2
      FORALLSITES(i, s)
        mat_copy(&(staple[i]), &(tempmat[i]));
      for (j = 0; j < x_dist; j++)
        shiftmat(tempmat, tempmat2, goffset[XUP]);

      // Evaluate potential at this separation
      wloop = 0.0;
      FORALLSITES(i, s)
        wloop += (double)realtrace(&(staple[i]), &(tempmat[i]));
      g_doublesum(&wloop);

      node0_printf("POLAR_LOOP %d %d %.6g\n", x_dist, t_dist, wloop / volume);

      // As we increment x, shift in x direction
      shiftmat(tempmat, tempmat2, goffset[XUP]);
    } // x_dist
  } // t_dist
}
// -----------------------------------------------------------------
