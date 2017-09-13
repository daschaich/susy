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

void wflow() {
  register int i, dir;
  register site *s;
  
  int j, istep;
  double t = 0.0, E, tSqE, old_tSqE, der_tSqE;
  double ssplaq, stplaq, plaq, check;
  
  // Go
  for (istep = 0; fabs(t) <  fabs(tmax) - 0.5 * fabs(epsilon); istep++) {
    stout_step_rk();
    t += epsilon;
    
    // Find 8F_munu = sum_{clover} (U - Udag)
    // Subtract the (lattice artifact?) trace at each lattice site
    make_field_strength();
    
    // Compute t^2 E and its slope
    E = 0.0;
    FORALLSITES(i, s) {
      for (dir = 0; dir < 10; dir++)
        E -= (double)realtrace_nn(&(FS[dir]), &(FS[dir])); // Why minus sign? Why no FS dagger?
    }
    g_doublesum(&E);
    E /= (volume * 64.0); // Normalization factor of 1/8 for each F_munu
    tSqE = t * t * E;
    der_tSqE = fabs(t) * (tSqE - old_tSqE) / fabs(epsilon);
    // Any negative signs in t and epsilon should cancel out anyway...
    
    // Check with plaquette
    plaquette(&ssplaq, &stplaq);
    plaq = 0.5 * (ssplaq + stplaq);
    check = 20.0 * t * t * ((double)NCOL - plaq); // This is a guess to numerical factor
    
    node0_printf("WFLOW %g %g %g %g %g %g \n",
                 t, plaq, E, tSqE, der_tSqE, check);
  }
}
