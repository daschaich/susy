// -----------------------------------------------------------------
// Construct HYP-smeared links, without projection or unitarization
// Basically we just want Omega = (1 - alpha) * U + ftmp2 * Gamma
//   where Gamma is the sum of m staples and ftmp2 = alpha / m
// We will smear the three spatial directions into the diagonal link
// !!! For now it is ignored
// Reverse-engineered from nHYP code
// References: hep-lat/0103029, hep-lat/0702028
#include "generic_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dir1 is the direction of the original link
// dir2 is the other direction that defines the staple
// Use gather offsets to handle all five links!
// stp must be cleared before being this function is called!
void staple_hyp(int dir1, int dir2, matrix_f *lnk1,
                matrix_f *lnk2, matrix_f *stp) {

  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  matrix_f tmat1, tmat2;

  // Get blocked_link[dir2] from direction dir1
  tag0 = start_gather_field(lnk2, sizeof(matrix_f), goffset[dir1],
                            EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir1] from direction dir2
  tag1 = start_gather_field(lnk1, sizeof(matrix_f), goffset[dir2],
                            EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat,
  // then gathered to x
  FORALLSITES(i, s)
    mult_an_f(lnk2 + i, lnk1 + i, &(tempmat[i]));

  wait_gather(tag0);
  wait_gather(tag1);

  // Finish lower staple
  FORALLSITES(i, s) {
    mult_nn_f(&(tempmat[i]), (matrix_f *)gen_pt[0][i], &tmat1);
    mat_copy_f(&tmat1, &(tempmat[i]));
  }

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmat, sizeof(matrix_f),
                            goffset[dir2] + 1, EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_nn_f(lnk2 + i, (matrix_f *)gen_pt[1][i], &tmat1);
    mult_na_f(&tmat1, (matrix_f *)gen_pt[0][i], &tmat2);
    add_matrix_f(stp + i, &tmat2, stp + i);
  }

  // Finally add the lower staple
  wait_gather(tag2);
  FORALLSITES(i, s)
    add_matrix_f(stp + i, (matrix_f *)gen_pt[2][i], stp + i);

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_hyp1() {
  register int dir1, dir2, i;
  register site *s;
  Real ftmp1, ftmp2;
  matrix_f Omega, tmat;

  ftmp1 = alpha_smear[2] / (2 * (1 - alpha_smear[2]));
  ftmp2 = (1 - alpha_smear[2]);

  // dir1 is the direction of the original linkf
  // dir2 is the other direction that defines the staple
  // Only smear spatial staples into diagonal link
  for (dir1 = XUP; dir1 < NUMLINK; dir1++) {
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir1 == DIR_5 && dir2 == TUP)
        continue;   // Actually, we should be done at this point
      if (dir1 != dir2) {
       FORALLSITES(i, s)
         clear_mat_f(&(tempmat[i]));

        // Compute the staple
        staple_hyp(dir1, dir2, thin_link[dir1],
                   thin_link[dir2], tempmat);

        FORALLSITES(i, s) {
          // Make Omega
          scalar_mult_add_matrix_f(thin_link[dir1] + i,
                                   tempmat + i, ftmp1, &tmat);
          scalar_mult_matrix_f(&tmat, ftmp2, &Omega);
          hyplink1[dir2][dir1][i] = Omega;
        }
      }
    } // Loop over dir2
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_hyp2() {
  register int dir, dir2, dir3, dir4, i;
  register site *s;
  Real ftmp1, ftmp2;
  matrix_f Omega, tmat;

  ftmp1 = alpha_smear[1] / (4 * (1 - alpha_smear[1]));
  ftmp2 = (1 - alpha_smear[1]);

  for (dir = XUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir == DIR_5 && dir2 == TUP)
        continue;   // Actually, we should be done at this point
      if (dir2 != dir) {
        FORALLSITES(i, s)
          clear_mat_f(tempmat + i);

        // Compute the staple
        for (dir3 = XUP; dir3 <= TUP; dir3++) {
          if (dir3 != dir && dir3 != dir2) {
            for (dir4 = XUP; dir4 <= TUP; dir4++) {
              if (dir4 != dir && dir4 != dir2 && dir4 != dir3)
                break;
            }
            staple_hyp(dir, dir3, hyplink1[dir4][dir],
                       hyplink1[dir4][dir3], tempmat);
          }
        }

        FORALLSITES(i, s) {
          // Make Omega
          scalar_mult_add_matrix_f(thin_link[dir] + i,
                                   tempmat + i, ftmp1, &tmat);

          scalar_mult_matrix_f(&tmat, ftmp2, &Omega);
          hyplink2[dir2][dir][i] = Omega;
        }
      }
    } // Loop over dir2
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_hyp3() {
  register int dir, dir2, i;
  register site *s;
  Real ftmp1, ftmp2;
  matrix_f Omega, tmat;

  ftmp1 = alpha_smear[0] / (6 * (1 - alpha_smear[0]));
  ftmp2 = 1 - alpha_smear[0];

  for (dir = XUP; dir < NUMLINK; dir++) {
    FORALLSITES(i, s)
      clear_mat_f(&tempmat[i]);

    // Compute the staple
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir == DIR_5 && dir2 == TUP)
        continue;   // Actually, we should be done at this point
      if (dir2 != dir)
        staple_hyp(dir, dir2, hyplink2[dir2][dir],
                   hyplink2[dir][dir2], tempmat);
    }

    FORALLSITES(i, s) {
      // Make Omega
      scalar_mult_add_matrix_f(thin_link[dir] + i,
                               tempmat + i, ftmp1, &tmat);

      scalar_mult_matrix_f(&tmat, ftmp2, &Omega);
      smeared_link[dir][i] = Omega;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Do three smearing levels to construct nHYP links
// Overwrite s->linkf and save original values in thin_link field
void block_hyp() {
  register int i, dir;
  register site *s;

  for (dir = XUP; dir < NUMLINK; dir++) {
    FORALLSITES(i, s)
      mat_copy_f(&(s->linkf[dir]), &(thin_link[dir][i]));
  }

  block_hyp1();
  block_hyp2();
  block_hyp3();

  for (dir = XUP; dir < NUMLINK; dir++) {
    FORALLSITES(i, s)
      mat_copy_f(&(smeared_link[dir][i]), &(s->linkf[dir]));
  }
}
// -----------------------------------------------------------------
