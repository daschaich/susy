// -----------------------------------------------------------------
// // Run Wilson flow 
// // Based on the method adopted in 1203.4469
// // Runge--Kutta coefficients computed from Eq. 2.4 of arXiv:1203.4469
// // !!! arXiv 1006.4518 provides useful details 
#include "susy_includes.h"

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


