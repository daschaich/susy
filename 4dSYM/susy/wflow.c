// -----------------------------------------------------------------
// // Run Wilson flow 
// // Based on the method adopted in 1203.4469
// // Runge--Kutta coefficients computed from Eq. 2.4 of arXiv:1203.4469
// // !!! arXiv 1006.4518 provides useful details 
#include "susy_includes.h"


// -----------------------------------------------------------------
// // Calculate newU = exp(Q).U, overwriting s->link
// // Here Q is traceless anti-hermitian 
// // Go to eighth order in the exponential:
// //   exp(Q) * U = (1 + Q + Q^2/2 + Q^3/6 ...) * U
// //              = U + Q*(U + (Q/2)*(U + (Q/3)*( ... )))
void exp_mult(double eps) {
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
      add_matrix(link, &tmat, &(s->link[dir]));
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
// Sum staples for direction dir over all other directions
void staple(matrix *staple[NDIMS]) {
  register int i;
  register site *s;
  int dir, dir2;

  FORALLUPDIR(dir) {
    FORALLSITES(i, s)
      clear_mat(&(stp[dir][i]));

    FORALLUPDIR(dir2) {
      if (dir == dir2)
        continue;
      directional_staple(dir, dir2);
    }
  }
}
// -----------------------------------------------------------------

// Make anti_hermitian : libraries/make_ahmat.c
// Scalar multiply and sum : libraries/s_m_a_mat.c
// Clear matrix : libraries/clear_mat.c

// -----------------------------------------------------------------
// Calculate A = A + f1 * Project_antihermitian_traceless(U.Sdag)
//           U = exp(f2 * A).U
// S is the Lie derivative of the action being used to flow
void update_flow(double f1, double f2) {
  register int i, dir;
  register site *s;
  matrix tmat;
  anti_hermitmat tmat_ah;

  staple(S);

  FORALLUPDIR(dir) { 
    FORALLSITES(i, s) {
    mult_na(&(s->link[dir]), &(S[dir][i]), &tmat);
    make_anti_hermitian(&tmat, &tmat_ah);
    // A += f1 * U.S
    scalar_mult_sum(&tmat_ah, (Real)f1, &(A[dir][i]));
    }
   exp_mult(f2); // U = exp(f2 * A).U
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void stout_step_rk() {
  register int i, dir;
  register site *s;

  // Clear A, just in case
  FORALLSITES(i, s) {
    FORALLUPDIR(dir)
    clear_mat(&A[dir][i]);
  }

  update_flow(17.0 * epsilon / 36.0, -9.0 / 17.0);
  update_flow(-8.0 * epsilon / 9.0, 1.0);
  update_flow(3.0 * epsilon / 4.0, -1.0);
}
// -----------------------------------------------------------------

void wflow(){






}


