// -----------------------------------------------------------------
// Either stout smearing or APE-like smearing without SU(N) projection
// For stout smearing see Morningstar and Peardon, hep-lat/0311018
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate newU = exp(eps * Q).U, overwriting s->link
// Here Q is the traceless anti-hermitian lattice field from stout_smear
// Go to eighth order in the exponential:
//   exp(x) * U = (1 + x + x^2/2 + x^3/6 ...) * U
//              = U + x*(U + (x/2)*(U + (x/3)*( ... )))
void exp_mult(Real eps) {
  register int i, dir;
  register site *s;
  register Real t2, t3, t4, t5, t6, t7, t8;
  matrix *link, tmat, tmat2, htemp;

  // Take divisions out of site loop (can't be done by compiler)
  t2 = eps / 2.0;
  t3 = eps / 3.0;
  t4 = eps / 4.0;
  t5 = eps / 5.0;
  t6 = eps / 6.0;
  t7 = eps / 7.0;
  t8 = eps / 8.0;

  FORALLDIR(dir) {
    FORALLSITES(i, s) {
      uncompress_anti_hermitian(&(Q[dir][i]), &htemp);
      link = &(s->link[dir]);

      mult_nn(&htemp, link, &tmat);
      scalar_mult_add_matrix(link, &tmat, t8, &tmat2);

      mult_nn(&htemp, &tmat2, &tmat);
      scalar_mult_add_matrix(link, &tmat, t7, &tmat2);

      mult_nn(&htemp, &tmat2, &tmat);
      scalar_mult_add_matrix(link, &tmat, t6, &tmat2);

      mult_nn(&htemp, &tmat2, &tmat);
      scalar_mult_add_matrix(link, &tmat, t5, &tmat2);

      mult_nn(&htemp, &tmat2, &tmat);
      scalar_mult_add_matrix(link, &tmat, t4, &tmat2);

      mult_nn(&htemp, &tmat2, &tmat);
      scalar_mult_add_matrix(link, &tmat, t3, &tmat2);

      mult_nn(&htemp, &tmat2, &tmat);
      scalar_mult_add_matrix(link, &tmat, t2, &tmat2);

      mult_nn(&htemp, &tmat2, &tmat);
      scalar_mult_sum_matrix(&tmat, eps, link);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Accumulate staple sum in staple (add instead of overwrite)
// Must copy or project desired links to s->mom before calling
// dir is the direction of the original link
// dir2 is the other direction that defines the staple
// Use gather offsets to handle all five links!
// Use tempmat and tempmat2 for temporary storage
void directional_staple(int dir, int dir2) {
  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  matrix tmat;

  // Get mom[dir2] from direction dir
  tag0 = start_gather_site(F_OFFSET(mom[dir2]), sizeof(matrix),
                           goffset[dir], EVENANDODD, gen_pt[0]);

  // Get mom[dir] from direction dir2
  tag1 = start_gather_site(F_OFFSET(mom[dir]), sizeof(matrix),
                           goffset[dir2], EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat2,
  // then gathered to x
  FORALLSITES(i, s)
    mult_an(&(s->mom[dir2]), &(s->mom[dir]), &(tempmat2[i]));

  // Finish lower staple after gather is done
  wait_gather(tag0);
  FORALLSITES(i, s)
    mult_nn(&(tempmat2[i]), (matrix *)gen_pt[0][i], &(tempmat[i]));

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmat, sizeof(matrix),
                            goffset[dir2] + 1, EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  wait_gather(tag1);
  FORALLSITES(i, s) {
    mult_nn(&(s->mom[dir2]), (matrix *)gen_pt[1][i], &tmat);
    mult_na_sum(&tmat, (matrix *)gen_pt[0][i], &(staple[i]));
  }

  // Finally add the lower staple
  wait_gather(tag2);
  FORALLSITES(i, s)
    sum_matrix((matrix *)gen_pt[2][i], &(staple[i]));

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Nsmear steps of stout smearing, overwriting s->link
// Use staple, s->mom and Q for temporary storage
// tempmat and tempmat2 also used through directional_staple
void stout_smear(int Nsmear, double alpha) {
  register int i, n, dir, dir2;
  register site *s;
  matrix tmat;

  for (n = 0; n < Nsmear; n++) {
    FORALLSITES(i, s) {
      // Unmodified links -- no projection or determinant division
      FORALLDIR(dir)
        mat_copy(&(s->link[dir]), &(s->mom[dir]));
    }

    FORALLDIR(dir) {
      FORALLSITES(i, s)
        clear_mat(&(staple[i]));     // Initialize staple sum

      // Accumulate staple sum in staple
      FORALLDIR(dir2) {
        if (dir != dir2)
          directional_staple(dir, dir2);
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
    exp_mult(1.0);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Nsmear steps of APE smearing without unitarization, overwriting s->link
// Optionally project unsmeared links and / or staple sums
void APE_smear(int Nsmear, double alpha, int project) {
  register int i, n, dir, dir2;
  register site *s;
  Real tr, tr2;
  matrix tmat, tmat2;

  tr2 = 1.0 - alpha;
  tr = alpha / (8.0 * tr2);

  for (n = 0; n < Nsmear; n++) {
    FORALLSITES(i, s) {
      FORALLDIR(dir) {
        // Decide what to do with links before smearing
        // Polar project, divide out determinant, or nothing
        if (project == 1) {
//          det_project(&(s->link[dir]), &(s->mom[dir]));
          polar(&(s->link[dir]), &(s->mom[dir]), &tmat);
        }
        else
          mat_copy(&(s->link[dir]), &(s->mom[dir]));
      }
    }

    FORALLDIR(dir) {
      FORALLSITES(i, s)
        clear_mat(&(staple[i]));     // Initialize staple sum

      // Accumulate staple sum in staple
      FORALLDIR(dir2) {
        if (dir != dir2)
          directional_staple(dir, dir2);
      }

      // Combine (1 - alpha).link + (alpha / 8).staple
      // optionally projecting the staple, but not the end result
      FORALLSITES(i, s) {
//        if (project == 1)
//          det_project(&(staple[i]), &tmat2);
//        else
          mat_copy(&(staple[i]), &tmat2);

        scalar_mult_add_matrix(&(s->link[dir]), &tmat2, tr, &tmat);
        scalar_mult_matrix(&tmat, tr2, &(s->link[dir]));
      }
    }
  }
}
// -----------------------------------------------------------------
