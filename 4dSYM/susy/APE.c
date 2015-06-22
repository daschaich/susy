// -----------------------------------------------------------------
// Construct APE-smeared links without unitarization
// Overwrite s->linkf and save original values in thin_link field
// The application program must define and allocate staples stp
#include "susy_includes.h"

void APE_smear(int Nsmear, double alpha, int do_det) {
  register int i, n, dir, dir2;
  register site *s;
  Real tr, tr2;
  su3_matrix_f tmat, tmat2;

  tr = alpha / (6.0 * (1.0 - alpha));
  tr2 = 1.0 - alpha;

  for (dir = XUP; dir < NUMLINK; dir++) {
    FORALLSITES(i, s)
      su3mat_copy_f(&(s->linkf[dir]), &(thin_link[dir][i]));
  }

  for (n = 0; n < Nsmear; n++) {
    for (dir = XUP; dir < NUMLINK; dir++) {
      FORALLSITES(i, s)
        clear_su3mat_f(&(stp[dir][i]));    // Initialize staple sum

      // Compute staple sums
      for (dir2 = XUP; dir2 < NUMLINK; dir2++) {
        if (dir != dir2)                 // Accumulate staples
          directional_staple(dir, dir2, F_OFFSET(linkf[dir]),
                             F_OFFSET(linkf[dir2]), stp[dir]);
      }

      // Combine (1 - alpha).link + (alpha / 6).staple
      FORALLSITES(i, s) {
        if (do_det == 1)
          det_project(&(stp[dir][i]), &tmat);
        else
          su3mat_copy_f(&(stp[dir][i]), &tmat);

        scalar_mult_add_su3_matrix_f(&(s->linkf[dir]), &tmat, tr, &tmat);
        scalar_mult_su3_matrix_f(&tmat2, tr2, &(smeared_link[dir][i]));
      }
    }

    for (dir = XUP; dir < NUMLINK; dir++) {
      FORALLSITES(i, s)
        su3mat_copy_f(&(smeared_link[dir][i]), &(s->linkf[dir]));
    }
  }
}
// -----------------------------------------------------------------
