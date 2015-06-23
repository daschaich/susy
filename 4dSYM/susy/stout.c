// -----------------------------------------------------------------
// Construct stout-smeared links, following the standard procedure,
// Morningstar and Peardon, hep-lat/0311018
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Accumulate staple sum in staple (add instead of overwrite)
// Must copy or project desired links to s->mom before calling
// dir is the direction of the original link
// dir2 is the other direction that defines the staple
// Use gather offsets to handle all five links!
// Use tempmat1 for temporary storage
void directional_staple(int dir, int dir2) {
  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  su3_matrix_f tmat, tmat2, *mat0, *mat1;

  // Get mom[dir2] from direction dir
  tag0 = start_gather_site(F_OFFSET(mom[dir2]), sizeof(su3_matrix_f),
                           goffset[dir], EVENANDODD, gen_pt[0]);

  // Get mom[dir] from direction dir2
  tag1 = start_gather_site(F_OFFSET(mom[dir]), sizeof(su3_matrix_f),
                           goffset[dir2], EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat1,
  // then gathered to x
  FORALLSITES(i, s)
    mult_su3_an_f(&(s->mom[dir2]), &(s->mom[dir]), &(tempmat1[i]));

  // Finish lower staple after gather is done
  wait_gather(tag0);
  FORALLSITES(i, s) {
    mat0 = (su3_matrix_f *)gen_pt[0][i];
    mult_su3_nn_f(&(tempmat1[i]), mat0, &tmat);
    su3mat_copy_f(&tmat, &(tempmat1[i]));
  }

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmat1, sizeof(su3_matrix_f),
                            goffset[dir2] + 1, EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  wait_gather(tag1);
  FORALLSITES(i, s) {
    mat0 = (su3_matrix_f *)gen_pt[0][i];
    mat1 = (su3_matrix_f *)gen_pt[1][i];
    mult_su3_nn_f(&(s->mom[dir2]), mat1, &tmat);
    mult_su3_na_f(&tmat, mat0, &tmat2);
    add_su3_matrix_f(&(staple[i]), &tmat2, &(staple[i]));
  }

  // Finally add the lower staple
  wait_gather(tag2);
  FORALLSITES(i, s)
    add_su3_matrix_f(&(staple[i]), (su3_matrix_f *)gen_pt[2][i], &(staple[i]));

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate newU = exp(A).U, overwriting s->linkf
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
      add_su3_matrix_f(link, &temp1, &(s->linkf[dir]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Do stout smearing
// Overwrite s->linkf
// Use staple, s->mom and Q for temporary storage
void stout_smear(int Nsmear, double alpha) {
  register int i, n, dir, dir2;
  register site *s;
  su3_matrix_f tmat;

  for (n = 0; n < Nsmear; n++) {
    FORALLSITES(i, s) {
      // Unmodified links -- no projection or determinant division
      for (dir = XUP; dir < NUMLINK; dir++)
        su3mat_copy_f(&(s->linkf[dir]), &(s->mom[dir]));
    }

    for (dir = XUP; dir < NUMLINK; dir++) {
      FORALLSITES(i, s)
        clear_su3mat_f(&(staple[i]));     // Initialize staple sum

      // Accumulate staple sum in staple
      for (dir2 = XUP; dir2 < NUMLINK; dir2++) {
        if (dir != dir2)
          directional_staple(dir, dir2);
      }

      // Multiply by alpha Udag, take traceless anti-hermitian part
      FORALLSITES(i, s) {
        mult_su3_na_f(&(staple[i]), &(s->linkf[dir]), &tmat);  // C.Udag
        scalar_mult_su3_matrix_f(&tmat, alpha, &(staple[i]));
        make_anti_hermitian(&(staple[i]), &(Q[dir][i]));
      }
    }

    // Do all exponentiations at once to reuse divisions
    // Overwrites s->linkf
    exp_mult();
  }
}
// -----------------------------------------------------------------
