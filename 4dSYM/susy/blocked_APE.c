// -----------------------------------------------------------------
// Construct APE-smeared links without unitarization
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Accumulate bl-strided staple sum in staple (add instead of overwrite)
// Must copy or project desired links to s->mom before calling
// dir is the direction of the original link
// dir2 is the other direction that defines the staple
// Use gather offsets to handle all five links!
// Use tempmat1 for temporary storage
void blocked_staple(int bl, int dir, int dir2) {
  register int i;
  register site *s;
  int j, d1[4] = {0, 0, 0, 0}, d2[4] = {0, 0, 0, 0};
  msg_tag *tag0, *tag1, *tag2;
  su3_matrix_f tmat, tmat2, *mat0, *mat1;

  for (j = 0; j < NDIMS; j++) {
    d1[j] = bl * offset[dir][j];
    d2[j] = bl * offset[dir2][j];
  }

  // Get mom[dir2] from displacement d1
  // Can only have one general gather at once...
  tag0 = start_general_gather_site(F_OFFSET(mom[dir2]), sizeof(su3_matrix_f),
                                   d1, EVENANDODD, gen_pt[0]);
  wait_general_gather(tag0);

  //  blocked_link[dir] from displacement d2
  tag1 = start_general_gather_site(F_OFFSET(mom[dir]), sizeof(su3_matrix_f),
                                   d2, EVENANDODD, gen_pt[1]);

  // Compute the lower staple while we wait for the gather
  // The lower staple is prepared at x-dir2 and stored in staple,
  // then gathered to x
  FORALLSITES(i, s) {
    mat0 = (su3_matrix_f *)gen_pt[0][i];
    mult_su3_an_f(&(s->mom[dir2]), &(s->mom[dir]), &(tempmat1[i]));
    mult_su3_nn_f(&(tempmat1[i]), mat0, &tmat);
    su3mat_copy_f(&tmat, &(tempmat1[i]));
  }

  // Gather staple from direction -dir2 to "home" site
  wait_general_gather(tag1);
  for (j = 0; j < NDIMS; j++)
    d2[j] *= -1;
  tag2 = start_general_gather_field(tempmat1, sizeof(su3_matrix_f), d2,
                                    EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mat0 = (su3_matrix_f *)gen_pt[0][i];
    mat1 = (su3_matrix_f *)gen_pt[1][i];
    mult_su3_nn_f(&(s->mom[dir2]), mat1, &tmat);
    mult_su3_na_f(&tmat, mat0, &tmat2);
    add_su3_matrix_f(&(staple[i]), &tmat, &(staple[i]));
  }

  // Finally add the lower staple
  wait_general_gather(tag2);
  FORALLSITES(i, s)
    add_su3_matrix_f(&(staple[i]), (su3_matrix_f *)gen_pt[2][i], &(staple[i]));

  cleanup_general_gather(tag0);
  cleanup_general_gather(tag1);
  cleanup_general_gather(tag2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Overwrite s->linkf
void blocked_APE(int Nsmear, double alpha, int project, int block) {
  register int i, n, dir, dir2;
  register site *s;
  int j, bl = 2;
  Real tr, tr2;
  su3_matrix_f tmat, tmat2;

  tr = alpha / (6.0 * (1.0 - alpha));
  tr2 = 1.0 - alpha;

  // Set number of links to stride, bl = 2^block
  // Allow sanity check of reproducing ploop() with this routine
  for (j = 1; j < block; j++)
    bl *= 2;
  if (block <= 0)
    bl = 1;

  for (n = 0; n < Nsmear; n++) {
    FORALLSITES(i, s) {
      for (dir = XUP; dir < NUMLINK; dir++) {
        // Decide what to do with links before smearing
        // Polar project, divide out determinant, or nothing
        if (project == 1) {
//          det_project(&(s->linkf[dir]), &tmat);
          polar(&(s->linkf[dir]), &tmat);
        }
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
          blocked_staple(bl, dir, dir2);
      }

      // Combine (1 - alpha).link + (alpha / 6).staple
      FORALLSITES(i, s) {
//        if (project == 1)
//          det_project(&(staple[i]), &tmat2);
//        else
          su3mat_copy_f(&(staple[i]), &tmat2);

        scalar_mult_add_su3_matrix_f(&(s->linkf[dir]), &tmat2, tr, &tmat);
        scalar_mult_su3_matrix_f(&tmat, tr2, &(s->linkf[dir]));
      }
    }
  }
}
// -----------------------------------------------------------------
