// -----------------------------------------------------------------
// Static potential for all displacements
// Evaluate in different spatial dirs to check rotational invariance
// Must gauge fix to Coulomb gauge before calling
// This version computes spatial correlators of temporal products
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void hvy_pot() {
  register int i;
  register site *s;
  int t_dist, x_dist;
  Real frac = -1.0 / (Real)NCOL;
  double wloop, detloop;
  complex det_wloop, c_loop, c1, c2, mult;
  su3_matrix_f tmat;
  msg_tag *mtag = NULL;
  field_offset oldmat, newmat, tt;

  node0_printf("hvy_pot: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  // Use tempmat1 to construct linear product from each point
  for (t_dist = 1; t_dist <= MAX_T; t_dist++) {
    if (t_dist == 1) {
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->linkf[TUP]), &(s->tempmat1));
    }
    else {
      mtag = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix_f),
                               goffset[TUP], EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s) {
        mult_su3_nn_f(&(s->linkf[TUP]), (su3_matrix_f *)gen_pt[0][i],
                      &(s->staple));
      }
      cleanup_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->staple), &(s->tempmat1));
    }
    // Now tempmat1 is product of t_dist links at each x
    oldmat = F_OFFSET(tempmat2);
    newmat = F_OFFSET(staple);    // Will switch these two

    for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
      // Evaluate potential at this separation
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->tempmat1), (su3_matrix_f *)F_PT(s, oldmat));
      // Evaluate potential at this separation
      wloop = 0.0;
      detloop = 0.0;
      FORALLSITES(i, s) {
        // Compute the actual Coulomb gauge Wilson loop product
        mult_su3_na_f(&(s->tempmat1),
                      (su3_matrix_f *)F_PT(s, oldmat), &tmat);
        c_loop = trace_su3_f(&tmat);
        wloop += c_loop.real;

        // Divide out the det raised to fractional power as
        // x^y = exp[y log x]
        det_wloop = find_det(&tmat);
        c1 = clog(&det_wloop);
        CMULREAL(c1, frac, c2);
        mult = cexp(&c2);
        CMUL(c_loop, mult, c1);
        detloop += c1.real;
      }
      g_doublesum(&wloop);
      g_doublesum(&detloop);
      node0_printf("POT_LOOP %d %d %.6g\n",
                   x_dist, t_dist, wloop / volume);
      node0_printf("D_LOOP   %d %d %.6g\n",
                   x_dist, t_dist, detloop / volume);

      // As we increment z, shift in z direction
      shiftmat(oldmat, newmat, goffset[XUP]);
      tt = oldmat;
      oldmat = newmat;
      newmat = tt;
    } // x_dist
  } // t_dist
}
// -----------------------------------------------------------------
