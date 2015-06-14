// -----------------------------------------------------------------
// Construct APE-smeared links
// The application program must define and allocate staples stp
#include "generic_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Do APE smearing
// Overwrite s->linkf and save original values in thin_link field
void APE_smear(int Nsmear, double alpha) {
  register int i, n, dir, dir2;
  register site *s;
  Real tr, tr2;
  su3_matrix_f tmat;

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
        scalar_mult_add_su3_matrix_f(&(s->linkf[dir]), &(stp[dir][i]),
                                     tr, &tmat);
        scalar_mult_su3_matrix_f(&tmat, tr2, &(smeared_link[dir][i]));
      }
    }

    for (dir = XUP; dir < NUMLINK; dir++) {
      FORALLSITES(i, s)
        su3mat_copy_f(&(smeared_link[dir][i]), &(s->linkf[dir]));
    }
  }
}
// -----------------------------------------------------------------
