// -----------------------------------------------------------------
// Static potential for all displacements
// Evaluate in different spatial dirs to check rotational invariance
// Must gauge fix to Coulomb gauge before calling
// Gauge-fixed links unitarized via polar projection -- overwritten!!!
// This version computes spatial correlators of temporal products
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void hvy_pot_polar() {
  register int i;
  register site *s;
  int t_dist, x_dist;
  double wloop;
  su3_matrix_f tmat;
  msg_tag *tag;
  field_offset oldmat, newmat, tt;

  node0_printf("hvy_pot_polar: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  FORALLSITES(i, s) {
   // Polar projection of gauge-fixed links
   // To be multiplied together after projecting
   // !!! Overwrites links
   polar(&(s->linkf[TUP]), &tmat);
   su3mat_copy_f(&tmat, &(s->linkf[TUP]));
  }

  // Use tempmat1 to construct linear product from each point
  for (t_dist = 1; t_dist <= MAX_T; t_dist++) {
    if (t_dist == 1) {
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->linkf[TUP]), &(s->tempmat1));
    }
    else {
      tag = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix_f),
                              goffset[TUP], EVENANDODD, gen_pt[0]);
      wait_gather(tag);
      FORALLSITES(i, s) {
        mult_su3_nn_f(&(s->linkf[TUP]), (su3_matrix_f *)gen_pt[0][i],
                      &(s->staple));
      }
      cleanup_gather(tag);
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->staple), &(s->tempmat1));
    }
    // Now tempmat1 is product of t_dist links at each x
    oldmat = F_OFFSET(tempmat2);
    newmat = F_OFFSET(staple);    // Will switch these two

    for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
      // Evaluate potential at this x
      wloop = 0.0;
      FORALLSITES(i, s) {
        wloop += (double)realtrace_su3_f(&(s->tempmat1),
                          (su3_matrix_f *)F_PT(s, oldmat));
      }
      g_doublesum(&wloop);
      node0_printf("POLAR_LOOP %d %d %.6g\n",
                   x_dist, t_dist, wloop / volume);

      // As we increment x, shift in x direction
      shiftmat(oldmat, newmat, XUP);
      tt = oldmat;
      oldmat = newmat;
      newmat = tt;
    } // x dist
  } // t_dist
}
// -----------------------------------------------------------------
