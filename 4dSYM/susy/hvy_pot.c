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
  int t_dist, x_dist, y_dist, z_dist;
  int d[4] = {0, 0, 0, 0};
  Real frac = -1.0 / (Real)NCOL;
  double wloop, detloop;
  complex det_wloop, c_loop, c1, c2, mult;
  su3_matrix_f tmat;
  msg_tag *mtag = NULL;

  node0_printf("hvy_pot: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  // Use tempmat1 to construct linear product from each point
  for (t_dist = 1; t_dist <= MAX_T; t_dist++) {
    if (t_dist == 1) {
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->linkf[DIR_5]), &(s->tempmat1));
    }
    else {
      mtag = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix_f),
                               goffset[DIR_5], EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s) {
        mult_su3_nn_f(&(s->linkf[DIR_5]), (su3_matrix_f *)gen_pt[0][i],
                      &(s->staple));
      }
      cleanup_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy_f(&(s->staple), &(s->tempmat1));
    }

    // Now tempmat1 is product of t_dist links at each (x, y, z)
    // Construct and print correlator, using general gathers for now
    for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
      d[XUP] = x_dist;
      for (y_dist = 0; y_dist <= MAX_X; y_dist++) {
        d[YUP] = y_dist;
        for (z_dist = 0; z_dist <= MAX_X; z_dist++) {
          d[ZUP] = z_dist;
          mtag = start_general_gather_site(F_OFFSET(tempmat1),
                                           sizeof(su3_matrix_f), d,
                                           EVENANDODD, gen_pt[0]);
          wait_general_gather(mtag);

          // Evaluate potential at this separation
          wloop = 0.0;
          detloop = 0.0;
          FORALLSITES(i, s) {
            // Compute the actual Coulomb gauge Wilson loop product
            mult_su3_na_f(&(s->tempmat1), (su3_matrix_f *)gen_pt[0][i], &tmat);
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
          node0_printf("POT_LOOP %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist, wloop / volume);
          node0_printf("D_LOOP   %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist, detloop / volume);
          cleanup_general_gather(mtag);
        } // z_dist
      } // y_dist
    } // x_dist
  } // t_dist
}
// -----------------------------------------------------------------
