// -----------------------------------------------------------------
// Static potential for all displacements
// Must gauge fix to Coulomb gauge before calling
// This version computes spatial correlators of temporal products
// Use tempmat, tempmat2 and staple for temporary storage
#include "susy_includes.h"

void hvy_pot(int do_det) {
  register int i;
  register site *s;
  int j, t_dist, x_dist;
  double wloop;
  complex tc;
  matrix tmat, tmat2;
  msg_tag *mtag = NULL;

  node0_printf("hvy_pot: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

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
      FORALLSITES(i, s) {
        // Compute the actual Coulomb gauge Wilson loop product
        mult_na(&(staple[i]), &(tempmat[i]), &tmat);

        if (do_det == 1)
          det_project(&tmat, &tmat2);
        else
          mat_copy(&tmat, &tmat2);

        tc = trace(&tmat2);
        wloop += tc.real;
      }
      g_doublesum(&wloop);

      if (do_det == 1) {  // Braces fix compiler error
        node0_printf("D_LOOP   ");
      }
      else
        node0_printf("POT_LOOP ");
      node0_printf("%d %d %.6g\n", x_dist, t_dist, wloop / volume);

      // As we increment x, shift in x direction
      shiftmat(tempmat, tempmat2, goffset[XUP]);
    } // x_dist
  } // t_dist
}
// -----------------------------------------------------------------
