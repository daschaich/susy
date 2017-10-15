// -----------------------------------------------------------------
// Run Wilson flow
// Runge--Kutta coefficients computed from Eq. 2.4 of arXiv:1203.4469
// !!! arXiv 1006.4518 provides useful details
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Clear A
void clear_antiH(anti_hermitmat *a) {
  int i;
  for (i = 0; i < NCOL; i++)
    a->im_diag[i] = 0.0;

  for (i = 0; i < NCOL * (NCOL - 1) / 2; i++)
    a->m[i] = cmplx(0.0, 0.0);
}

// c <-- c + s * b (output is always last)
void scalar_mult_sum_antiH(anti_hermitmat *b, Real s, anti_hermitmat *c) {
  int i;
  for (i = 0; i < NCOL; i++)
    c->im_diag[i] += s * b->im_diag[i];

  for (i = 0; i < NCOL * (NCOL - 1) / 2; i++) {
    c->m[i].real += s * b->m[i].real;
    c->m[i].imag += s * b->m[i].imag;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Sum staples for direction dir over all other directions
void compute_staple() {
  register int i;
  register site *s;
  int dir, dir2;

  FORALLSITES(i, s)
    FORALLDIR(dir) {
      mat_copy(&(s->link[dir]), &(s->mom[dir]));
  }

  FORALLDIR(dir) {
    FORALLSITES(i, s)
      clear_mat(&(staple[i]));

    FORALLDIR(dir2) {
      if (dir == dir2)
        continue;
      directional_staple(dir, dir2);
    }

    FORALLSITES(i, s)
      mat_copy(&(staple[i]), &(S[dir][i]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate Q = Q + f1 * Project_antihermitian_traceless(U.Sdag)
//           U = exp(f2 * Q).U
// S is the Lie derivative of the action being used to flow
void update_flow(Real f1, Real f2) {
  register int i, dir;
  register site *s;
  matrix tmat;
  anti_hermitmat tmat_ah;

  // This is where we can change the 'flow action'
  compute_staple();

  FORALLDIR(dir) {
    FORALLSITES(i, s) {
      mult_na(&(s->link[dir]), &(S[dir][i]), &tmat);
      make_anti_hermitian(&tmat, &tmat_ah);
      // Q += f1 * U.S
      scalar_mult_sum_antiH(&tmat_ah, f1, &(Q[dir][i]));
    }
  }
  exp_mult(f2);             // U = exp(f2 * Q).U
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void stout_step_rk() {
  register int i, dir;
  register site *s;

  // Clear Q, just in case
  FORALLSITES(i, s) {
    FORALLDIR(dir)
      clear_antiH(&(Q[dir][i]));
  }

  update_flow(17.0 * wflow_eps / 36.0, -9.0 / 17.0);
  update_flow(-8.0 * wflow_eps / 9.0, 1.0);
  update_flow(3.0 * wflow_eps / 4.0, -1.0);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void wflow() {
  register int i, dir;
  register site *s;
  int istep;
  double t = 0.0, E, tSqE, old_tSqE = 0.0, der_tSqE;
  double ssplaq, stplaq, plaq, check;

  // Go
  for (istep = 0; fabs(t) <  fabs(tmax) - 0.5 * fabs(wflow_eps); istep++) {
    stout_step_rk();
    t += wflow_eps;

    // Find 8F_munu = sum_{clover} (U - Udag)
    // Subtract the (lattice artifact?) trace at each lattice site
    make_field_strength();

    // Compute t^2 E and its slope
    E = 0.0;
    FORALLSITES(i, s) {
      for (dir = 0; dir < 10; dir++)
        E -= (double)realtrace_nn(&(FS[dir][i]), &(FS[dir][i])); // TODO: Check minus sign, no FS dagger
    }
    g_doublesum(&E);
    E /= (volume * 64.0); // Normalization factor of 1/8 for each F_munu
    tSqE = t * t * E;
    der_tSqE = fabs(t) * (tSqE - old_tSqE) / fabs(wflow_eps);
    // Any negative signs in t and wflow_eps should cancel out anyway...

    // Check with plaquette
    plaquette(&ssplaq, &stplaq);
    plaq = 0.5 * (ssplaq + stplaq);
    // TODO: Guessing numerical factor
    check = 20.0 * t * t * fabs((double)NCOL - plaq);

    old_tSqE = tSqE;   // The reason why 'der_tSqE' was exploding before ! 
    node0_printf("WFLOW %g %g %g %g %g %g\n",
                 t, plaq, E, tSqE, der_tSqE, check);
  }
}
// -----------------------------------------------------------------
