// -----------------------------------------------------------------
// Either stout smearing or APE-like smearing without SU(N) projection
// For stout smearing see Morningstar and Peardon, hep-lat/0311018
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Accumulate staple sum in staple (add instead of overwrite)
// Must copy or project desired links to s->mom before calling
// dir is the direction of the original link
// dir2 is the other direction that defines the staple
// Use gather offsets to handle all five links!
// Use tempmat1 and tempmat2 for temporary storage
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
  // The lower staple is prepared at x-dir2 and stored in tempmat2,
  // then gathered to x
  FORALLSITES(i, s)
    mult_su3_an_f(&(s->mom[dir2]), &(s->mom[dir]), &(tempmat2[i]));

  // Finish lower staple after gather is done
  wait_gather(tag0);
  FORALLSITES(i, s) {
    mat0 = (su3_matrix_f *)gen_pt[0][i];
    mult_su3_nn_f(&(tempmat2[i]), mat0, &(tempmat1[i]));
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
// Nsmear steps of stout smearing, overwriting s->linkf
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



// -----------------------------------------------------------------
// Nsmear steps of APE smearing without unitarization, overwriting s->linkf
// Optionally project unsmeared links and / or staple sums
void APE_smear(int Nsmear, double alpha, int project) {
  register int i, n, dir, dir2;
  register site *s;
  Real tr, tr2;
  su3_matrix_f tmat, tmat2;

  tr = alpha / (8.0 * (1.0 - alpha));
  tr2 = 1.0 - alpha;

  for (n = 0; n < Nsmear; n++) {
    FORALLSITES(i, s) {
      for (dir = XUP; dir < NUMLINK; dir++) {
        // Decide what to do with links before smearing
        // Polar project, divide out determinant, or nothing
        if (project == 1) {
//          det_project(&(s->linkf[dir]), &(s->mom[dir]));
          polar(&(s->linkf[dir]), &(s->mom[dir]), &tmat);
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
          directional_staple(dir, dir2);
      }

      // Combine (1 - alpha).link + (alpha / 8).staple
      // optionally projecting the staple, but not the end result
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
