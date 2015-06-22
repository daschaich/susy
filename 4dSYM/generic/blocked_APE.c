// -----------------------------------------------------------------
// Construct APE-smeared links without unitarization
// Overwrite s->linkf and save original values in thin_link field
// The application program must define and allocate staples stp
#include "generic_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Accumulate bl-strided staple sum in stp
// dir is the direction of the original link
// dir2 is the other direction that defines the staple
// Use gather offsets to handle all five links!
// Use tempmat1, tempmat2 and staple for temporary storage
void blocked_staple(int bl, int dir, int dir2, field_offset lnk1,
                    field_offset lnk2, su3_matrix_f *stp) {

  register int i;
  register site *s;
  int j, d1[4] = {0, 0, 0, 0}, d2[4] = {0, 0, 0, 0};
  msg_tag *tag;
  su3_matrix_f tmat1, tmat2;

  for (j = 0; j < NDIMS; j++) {
    d1[j] = bl * offset[dir][j];
    d2[j] = bl * offset[dir2][j];
  }

  // Get blocked_link[dir2] from displacement d1, save in tempmat2
  // Can only have one general gather at once...
  tag = start_general_gather_site(lnk2, sizeof(su3_matrix_f), d1,
                                  EVENANDODD, gen_pt[0]);
  wait_general_gather(tag);
  FORALLSITES(i, s)
    su3mat_copy_f((su3_matrix_f *)gen_pt[0][i], &(tempmat2[i]));
  cleanup_general_gather(tag);

  // Get blocked_link[dir] from displacement d2, save in tempmat1
  tag = start_general_gather_site(lnk1, sizeof(su3_matrix_f), d2,
                                  EVENANDODD, gen_pt[0]);
  wait_general_gather(tag);
  FORALLSITES(i, s)
    su3mat_copy_f((su3_matrix_f *)gen_pt[0][i], &(tempmat1[i]));
  cleanup_general_gather(tag);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in staple,
  // then gathered to x
  FORALLSITES(i, s)
    mult_su3_an_f((su3_matrix_f *)F_PT(s, lnk2), (su3_matrix_f *)F_PT(s, lnk1),
                  &(staple[i]));

  // Finish lower staple
  FORALLSITES(i, s) {
    mult_su3_nn_f(&(staple[i]), &(tempmat2[i]), &tmat1);
    su3mat_copy_f(&tmat1, &(staple[i]));
  }

  // Gather staple from direction -dir2 to "home" site
  d2[dir2] = -1 * bl;
  tag = start_general_gather_field(staple, sizeof(su3_matrix_f), d2,
                                   EVENANDODD, gen_pt[0]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_su3_nn_f((su3_matrix_f *)F_PT(s, lnk2), &(tempmat1[i]), &tmat1);
    mult_su3_na_f(&tmat1, &(tempmat2[i]), &tmat2);
    add_su3_matrix_f(stp + i, &tmat2, stp + i);
  }

  // Finally add the lower staple
  wait_general_gather(tag);
  FORALLSITES(i, s)
    add_su3_matrix_f(stp + i, (su3_matrix_f *)gen_pt[0][i], stp + i);
  cleanup_general_gather(tag);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void blocked_APE(int Nsmear, double alpha, int block) {
  register int i, n, dir, dir2;
  register site *s;
  int j, bl = 2;
  Real tr, tr2;
  su3_matrix_f tmat;

  tr = alpha / (6.0 * (1.0 - alpha));
  tr2 = 1.0 - alpha;

  // Set number of links to stride, bl = 2^block
  // Allow sanity check of reproducing ploop() with this routine
  for (j = 1; j < block; j++)
    bl *= 2;
  if (block <= 0)
    bl = 1;

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
          blocked_staple(bl, dir, dir2, F_OFFSET(linkf[dir]),
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
