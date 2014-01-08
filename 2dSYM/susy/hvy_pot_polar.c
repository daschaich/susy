// -----------------------------------------------------------------
// Heavy quark potential modified to do all displacements
// Evaluate in different spatial dirs to check rotational invariance
// Must gauge fix to Coulomb gauge before calling to do spatial segments
// Gauge-fixed links unitarized via polar projection
// !!! Links overwritten -- only run last on saved configurations
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Argument tells whether to use ordinary links or fat links
void hvy_pot_polar(field_offset links) {
  register int i;
  register site *s;
  int t_dist, x_dist, y_dist, z_dist;
  double wloop;
  su3_matrix_f tmat1, tmat2;
  msg_tag *mtag0;
  field_offset oldmat, newmat, tt;

  node0_printf("hvy_pot_polar(): MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  FORALLSITES(i, s) {
   // Polar projection of gauge-fixed links
   // To be multiplied together after projecting
   // !!! Links overwritten
   polar(&(((su3_matrix_f *)(F_PT(s, links)))[TUP]), &tmat1, &tmat2);
   su3mat_copy_f(&tmat1, &(((su3_matrix_f *)(F_PT(s, links)))[TUP]));
  }

  // Use tempmat1 to construct t-direction path from each point
  for (t_dist = 1; t_dist <= MAX_T; t_dist ++) {
    if (t_dist == 1) {
      FORALLSITES(i, s) {
        su3mat_copy_f(&(((su3_matrix_f *)(F_PT(s, links)))[TUP]),
                      &(s->tempmat1));
      }
    }
    else {
      mtag0 = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix_f),
                                TUP, EVENANDODD, gen_pt[0]);
      wait_gather(mtag0);
      FORALLSITES(i, s) {
        mult_su3_nn_f(&(((su3_matrix_f *)(F_PT(s, links)))[TUP]),
                      (su3_matrix_f *)gen_pt[0][i], &(s->staple));
      }
      cleanup_gather(mtag0);
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->staple), &(s->tempmat1));
    }
    // Now tempmat1 is path of length t_dist in TUP direction
    oldmat = F_OFFSET(tempmat2);
    newmat = F_OFFSET(staple);    // Will switch these two

    for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
      for (y_dist = 0; y_dist <= MAX_X; y_dist++) {
        // Gather from spatial dirs, compute products of paths
        FORALLSITES(i, s)
          su3mat_copy_f(&(s->tempmat1), (su3_matrix_f *)F_PT(s, oldmat));
        for (i = 0; i < x_dist; i++) {
          shiftmat(oldmat, newmat, XUP);
          tt = oldmat;
          oldmat = newmat;
          newmat = tt;
        }
        for (i = 0; i < y_dist; i++) {
          shiftmat(oldmat, newmat, YUP);
          tt = oldmat;
          oldmat = newmat;
          newmat = tt;
        }
        for (z_dist = 0; z_dist <= MAX_X; z_dist++) {
          // Evaluate potential at this separation
          wloop = 0.0;
          FORALLSITES(i, s) {
            wloop += (double)realtrace_su3_f(&(s->tempmat1),
                             (su3_matrix_f *)F_PT(s, oldmat));
          }
          g_doublesum(&wloop);
          node0_printf("POLAR_LOOP %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist, wloop / volume);

          // As we increment z, shift in z direction
          shiftmat(oldmat, newmat, ZUP);
          tt = oldmat;
          oldmat = newmat;
          newmat = tt;
        } // z_dist
      } // y dist
    } // x dist
  } // t_dist
}
// -----------------------------------------------------------------
