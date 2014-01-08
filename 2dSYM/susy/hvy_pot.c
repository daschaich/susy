// -----------------------------------------------------------------
// Heavy quark potential modified to do all displacements
// Evaluate in different spatial dirs to check rotational invariance
// Must gauge fix to Coulomb gauge before calling to do spatial segments
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Argument tells whether to use ordinary links or fat links
void hvy_pot(field_offset links) {
  register int i;
  register site *s;
  int t_dist, x_dist, y_dist, z_dist;
  Real frac = -1.0 / (Real)NCOL;
  double wloop, detloop;
  complex det_wloop, c_loop, c1, c2, mult;
  su3_matrix_f tmat;
  msg_tag *mtag0;
  field_offset oldmat, newmat, tt;

  node0_printf("hvy_pot(): MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

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
          detloop = 0.0;
          FORALLSITES(i, s) {
            wloop += (double)realtrace_su3_f(&(s->tempmat1),
                             (su3_matrix_f *)F_PT(s, oldmat));

            // Compute the actual Coulomb gauge Wilson loop product
            mult_su3_na_f(&(s->tempmat1),
                          (su3_matrix_f *)F_PT(s, oldmat), &tmat);

            // Divide out the det raised to fractional power as
            // x^y = exp[y log x]
            c_loop = trace_su3_f(&tmat);
            det_wloop = find_det(&tmat);
            c1 = clog(&det_wloop);
            CMULREAL(c1, frac, c2);
            mult = cexp(&c2);
            CMUL(c_loop, mult, c1);
            detloop += c1.real;
          }
          g_doublesum(&wloop);
          g_doublesum(&detloop);
          node0_printf("POT_LOOP %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist, wloop / volume);
          node0_printf("D_LOOP   %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist, detloop / volume);

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
