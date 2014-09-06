// -----------------------------------------------------------------
// Construct stout-smeared links, following the standard procedure,
// Morningstar and Peardon, hep-lat/0311018
// The application program must define and allocate staples stp,
// and anti-hermitian matrices Q
#include "generic_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Accumulate staple sum in stp
// dir1 is the direction of the original link
// dir2 is the other direction that defines the staple
// Use gather offsets to handle all five links!
void directional_staple(int dir1, int dir2, field_offset lnk1,
                        field_offset lnk2, su3_matrix_f *stp) {

  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  su3_matrix_f tmat1, tmat2;

  // Get blocked_link[dir2] from direction dir1
  tag0 = start_gather_site(lnk2, sizeof(su3_matrix_f), goffset[dir1],
                            EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir1] from direction dir2
  tag1 = start_gather_site(lnk1, sizeof(su3_matrix_f), goffset[dir2],
                            EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat,
  // then gathered to x
  FORALLSITES(i, s)
    mult_su3_an_f((su3_matrix_f *)F_PT(s, lnk2), (su3_matrix_f *)F_PT(s, lnk1),
                  &(tempmat[i]));

  wait_gather(tag0);
  wait_gather(tag1);

  // Finish lower staple
  FORALLSITES(i, s) {
    mult_su3_nn_f(&(tempmat[i]), (su3_matrix_f *)gen_pt[0][i], &tmat1);
    su3mat_copy_f(&tmat1, &(tempmat[i]));
  }

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmat, sizeof(su3_matrix_f),
                            goffset[dir2] + 1, EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_su3_nn_f((su3_matrix_f *)F_PT(s, lnk2),
                  (su3_matrix_f *)gen_pt[1][i], &tmat1);
    mult_su3_na_f(&tmat1, (su3_matrix_f *)gen_pt[0][i], &tmat2);
    add_su3_matrix_f(stp + i, &tmat2, stp + i);
  }

  // Finally add the lower staple
  wait_gather(tag2);
  FORALLSITES(i, s)
    add_su3_matrix_f(stp + i, (su3_matrix_f *)gen_pt[2][i], stp + i);

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate smearedU = exp(A).U
// Goes to eighth order in the exponential:
//   exp(A) * U = ( 1 + A + A^2/2 + A^3/6 ...) * U
//              = U + A*(U + (A/2)*(U + (A/3)*( ... )))
void exp_mult() {
  register int i, dir;
  register site *s;
  register Real t2, t3, t4, t5, t6, t7, t8;
  su3_matrix_f *link, temp1, temp2, htemp;

  // Take divisions out of site loop (can't be done by compiler)
  t2 = 1.0 / 2.0;
  t3 = 1.0 / 3.0;
  t4 = 1.0 / 4.0;
  t5 = 1.0 / 5.0;
  t6 = 1.0 / 6.0;
  t7 = 1.0 / 7.0;
  t8 = 1.0 / 8.0;

  for (dir = XUP; dir < NUMLINK; dir++) {
    FORALLSITES(i, s) {
      uncompress_anti_hermitian(&(Q[dir][i]), &htemp);
      link = &(s->linkf[dir]);

      mult_su3_nn_f(&htemp, link, &temp1);
      scalar_mult_add_su3_matrix_f(link, &temp1, t8, &temp2);

      mult_su3_nn_f(&htemp, &temp2, &temp1);
      scalar_mult_add_su3_matrix_f(link, &temp1, t7, &temp2);

      mult_su3_nn_f(&htemp, &temp2, &temp1);
      scalar_mult_add_su3_matrix_f(link, &temp1, t6, &temp2);

      mult_su3_nn_f(&htemp, &temp2, &temp1);
      scalar_mult_add_su3_matrix_f(link, &temp1, t5, &temp2);

      mult_su3_nn_f(&htemp, &temp2, &temp1);
      scalar_mult_add_su3_matrix_f(link, &temp1, t4, &temp2);

      mult_su3_nn_f(&htemp, &temp2, &temp1);
      scalar_mult_add_su3_matrix_f(link, &temp1, t3, &temp2);

      mult_su3_nn_f(&htemp, &temp2, &temp1);
      scalar_mult_add_su3_matrix_f(link, &temp1, t2, &temp2);

      mult_su3_nn_f(&htemp, &temp2, &temp1);
      add_su3_matrix_f(link, &temp1, &(smeared_link[dir][i]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Do stout smearing
// Overwrite s->linkf and save original values in thin_link field
void block_stout(int Nstout, double rho) {
  register int i, n, dir, dir2;
  register site *s;
  su3_matrix_f tmat;

#ifdef TIMING
  TIC(0)
#endif

  for (dir = XUP; dir < NUMLINK; dir++) {
    FORALLSITES(i, s)
      su3mat_copy_f(&(s->linkf[dir]), &(thin_link[dir][i]));
  }

  for (n = 0; n < Nstout; n++) {
    for (dir = XUP; dir < NUMLINK; dir++) {
      FORALLSITES(i, s)
        clear_su3mat_f(&(stp[dir][i]));    // Initialize staple sum

      // Compute staple sums
      for (dir2 = XUP; dir2 < NUMLINK; dir2++) {
        if (dir != dir2)                 // Accumulate staples
          directional_staple(dir, dir2, F_OFFSET(linkf[dir]),
                             F_OFFSET(linkf[dir2]), stp[dir]);
      }

      // Multiply by rho Udag, take traceless anti-hermitian part
      FORALLSITES(i, s) {
        mult_su3_na_f(&(stp[dir][i]), &(s->linkf[dir]), &tmat);  // C.Udag
        scalar_mult_su3_matrix_f(&tmat, rho, &(stp[dir][i]));
        make_anti_hermitian(&(stp[dir][i]), &(Q[dir][i]));
      }
    }

    // Do all exponentiations at once to reuse divisions
    exp_mult();

    for (dir = XUP; dir < NUMLINK; dir++) {
      FORALLSITES(i, s)
        su3mat_copy_f(&(smeared_link[dir][i]), &(s->linkf[dir]));
    }
  }

#ifdef TIMING
  TOC(0, time_block_stout)
#endif
}
// -----------------------------------------------------------------
