// -----------------------------------------------------------------
// Construct APE-smeared links without unitarization
// Overwrite s->linkf
// Use staple and s->mom for temporary storage
#include "susy_includes.h"

void APE_smear(int Nsmear, double alpha, int project) {
  register int i, n, dir, dir2;
  register site *s;
  Real tr, tr2;
  su3_matrix_f tmat, tmat2;

  tr = alpha / (6.0 * (1.0 - alpha));
  tr2 = 1.0 - alpha;

  for (n = 0; n < Nsmear; n++) {
    FORALLSITES(i, s) {
      for (dir = XUP; dir < NUMLINK; dir++) {
        // Decide what to do with links before smearing
        // Polar project, divide out determinant, or nothing
        if (project == 1)
          polar(&(s->linkf[dir]), &(s->mom[dir]));
        else
          su3mat_copy_f(&(s->linkf[dir]), &(s->mom[dir]));
      }
    }

    for (dir = XUP; dir < NUMLINK; dir++) {
      FORALLSITES(i, s)
        clear_su3mat_f(&(staple[i]));     // Initialize staple sum

      // Accumulate staple sum in staple
      for (dir2 = XUP; dir2 < NUMLINK; dir2++) {
        if (dir != dir2)
          directional_staple(dir, dir2);
      }

      // Combine (1 - alpha).link + (alpha / 6).staple
      // Don't do anything to this link!
      FORALLSITES(i, s) {
//        if (project == 1)
//          polar(&(staple[i]), &tmat2);
//        else
          su3mat_copy_f(&(staple[i]), &tmat2);

        scalar_mult_add_su3_matrix_f(&(s->linkf[dir]), &tmat2, tr, &tmat);
        scalar_mult_su3_matrix_f(&tmat, tr2, &(s->linkf[dir]));
      }
    }
  }
}
// -----------------------------------------------------------------
