// -----------------------------------------------------------------
// Either stout smearing or APE-like smearing after RG blocking
// No SU(N) projection for APE-like smeared links
// For stout smearing see Morningstar and Peardon, hep-lat/0311018
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Accumulate strided staple sum in staple (add instead of overwrite)
// Must copy or project desired links to s->mom before calling
// dir is the direction of the original link
// dir2 is the other direction that defines the staple
// Use gather offsets to handle all five links!
// Use tempmat, tempmat2 and UpsiU[01] for temporary storage
void blocked_staple(int stride, int dir, int dir2) {
  register int i;
  register site *s;
  int j;
  matrix tmat, tmat2;

  // Copy links to tempmat and tempmat2 to be shifted
  FORALLSITES(i, s) {
    mat_copy_f(&(s->mom[dir]), &(tempmat[i]));
    mat_copy_f(&(s->mom[dir2]), &(tempmat2[i]));
  }

  // Get mom[dir2] from dir and mom[dir] from dir2, both with stride
  // This order may be easier on cache
  for (j = 0; j < stride; j++)
    shiftmat(tempmat2, UpsiU[0], goffset[dir]);
  for (j = 0; j < stride; j++)
    shiftmat(tempmat, UpsiU[0], goffset[dir2]);

  // Compute the lower staple at x - stride * dir2, store it in UpsiU[0],
  // and gather it to x, again with stride
  // then gathered to x
  FORALLSITES(i, s) {
    mult_an(&(s->mom[dir2]), &(s->mom[dir]), &tmat);
    mult_nn(&tmat, &(tempmat2[i]), &(UpsiU[0][i]));
  }
  for (j = 0; j < stride; j++)
    shiftmat(UpsiU[0], UpsiU[1], goffset[dir2] + 1);

  // Calculate upper staple, add both staples to sum
  FORALLSITES(i, s) {
    mult_nn(&(s->mom[dir2]), &(tempmat[i]), &tmat);
    mult_na_sum(&tmat, &(tempmat2[i]), &(staple[i]));
    sum_matrix(&(UpsiU[0][i]), &(staple[i]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Nsmear steps of stout smearing, overwriting s->link
// Use staple, s->mom and Q for temporary storage
// bl <= 0 (stride = 1) allows sanity check of reproducing stout_smear()
void blocked_stout(int Nsmear, double alpha, int bl) {
  register int i, n, dir, dir2;
  register site *s;
  int j, stride = 1;
  matrix tmat;

  // Set number of links to stride, 2^bl
  for (j = 0; j < bl; j++)
    stride *= 2;

  for (n = 0; n < Nsmear; n++) {
    FORALLSITES(i, s) {
      // Unmodified links -- no projection or determinant division
      for (dir = XUP; dir < NUMLINK; dir++)
        mat_copy_f(&(s->link[dir]), &(s->mom[dir]));
    }

    for (dir = XUP; dir < NUMLINK; dir++) {
      FORALLSITES(i, s)
        clear_mat(&(staple[i]));     // Initialize staple sum

      // Accumulate staple sum in staple
      for (dir2 = XUP; dir2 < NUMLINK; dir2++) {
        if (dir != dir2)
          blocked_staple(stride, dir, dir2);
      }

      // Multiply by alpha Udag, take traceless anti-hermitian part
      FORALLSITES(i, s) {
        mult_na(&(staple[i]), &(s->link[dir]), &tmat);  // C.Udag
        scalar_mult_matrix(&tmat, alpha, &(staple[i]));
        make_anti_hermitian(&(staple[i]), &(Q[dir][i]));
      }
    }

    // Do all exponentiations at once to reuse divisions
    // Overwrites s->link
    exp_mult();
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Nsmear steps of APE smearing without unitarization, overwriting s->link
// Optionally project unsmeared links and / or staple sums
// bl <= 0 (stride = 1) allows sanity check of reproducing APE_smear()
void blocked_APE(int Nsmear, double alpha, int project, int bl) {
  register int i, n, dir, dir2;
  register site *s;
  int j, stride = 1;
  Real tr, tr2;
  matrix tmat, tmat2;

  tr = alpha / (8.0 * (1.0 - alpha));
  tr2 = 1.0 - alpha;

  // Set number of links to stride, 2^bl
  for (j = 0; j < bl; j++)
    stride *= 2;

  for (n = 0; n < Nsmear; n++) {
    FORALLSITES(i, s) {
      for (dir = XUP; dir < NUMLINK; dir++) {
        // Decide what to do with links before smearing
        // Polar project, divide out determinant, or nothing
        if (project == 1) {
//          det_project(&(s->link[dir]), &(s->mom[dir]));
          polar(&(s->link[dir]), &(s->mom[dir]), &tmat);
        }
        else
          mat_copy_f(&(s->link[dir]), &(s->mom[dir]));
      }
    }

    for (dir = XUP; dir < NUMLINK; dir++) {
      FORALLSITES(i, s)
        clear_mat(&(staple[i]));     // Initialize staple sum

      // Accumulate staple sum in staple
      for (dir2 = XUP; dir2 < NUMLINK; dir2++) {
        if (dir != dir2)
          blocked_staple(stride, dir, dir2);
      }

      // Combine (1 - alpha).link + (alpha / 8).staple
      // optionally projecting the staple, but not the end result
      FORALLSITES(i, s) {
//        if (project == 1)
//          det_project(&(staple[i]), &tmat2);
//        else
          mat_copy_f(&(staple[i]), &tmat2);

        scalar_mult_add_matrix(&(s->link[dir]), &tmat2, tr, &tmat);
        scalar_mult_matrix(&tmat, tr2, &(s->link[dir]));
      }
    }
  }
}
// -----------------------------------------------------------------
