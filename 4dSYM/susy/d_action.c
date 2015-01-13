// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the CG should already have been run,
// so that the vector **sol contains (M_adjoint*M+shift[n])^(-1) * src
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge momentum contribution to the action
double d_hmom_action() {
  register int i, mu;
  register site *s;
  double sum = 0;

  for (mu = 0; mu < NUMLINK; mu++) {
    FORALLSITES(i, s)
      sum += (double)realtrace_su3_f(&(s->mom[mu]), &(s->mom[mu]));
  }
  g_doublesum(&sum);
//  node0_printf("gauge momentum %e\n",sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Include tunable coefficient C2 in the d^2 term
double d_gauge_action(int do_det) {
  register int i;
  register site *s;
  int index;
  double g_action = 0.0, norm = 0.5 * C2;
  complex tc;
  su3_matrix_f tmat, tmat2;

  compute_DmuUmu();   // Includes plaqdet if G is non-zero
  compute_Fmunu();
  FORALLSITES(i, s) {
    // d^2 term normalized by C2 / 2
    mult_su3_nn_f(&(DmuUmu[i]), &(DmuUmu[i]), &tmat);
    scalar_mult_su3_matrix_f(&tmat, norm, &tmat);

    // F^2 term
    for (index = 0; index < NPLAQ; index++) {
      mult_su3_an_f(&(Fmunu[index][i]), &(Fmunu[index][i]), &tmat2);
      scalar_mult_add_su3_matrix_f(&tmat, &tmat2, 2.0, &tmat);
    }

    if (do_det == 1)
      det_project(&tmat, &tmat2);
    else
      su3mat_copy_f(&tmat, &tmat2);

    tc = trace_su3_f(&tmat2);
    g_action += tc.real;
  }
  g_action *= kappa;
  g_doublesum(&g_action);
  return g_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Bosonic mass term to regulate flat directions -- note factor of kappa
double d_bmass_action() {
  register int i;
  register site *s;
  int mu;
  double sum = 0.0, dum;

  for (mu = XUP; mu < NUMLINK; mu++) {
    FORALLSITES(i, s) {
      dum = 1.0 / (double)NCOL;
      dum *= realtrace_su3_f(&(s->linkf[mu]), &(s->linkf[mu]));
      dum -= 1.0;
      sum += kappa * bmass * bmass * dum * dum;
    }
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Separate plaquette determinant contribution to the action
double d_det_action() {
  register int i;
  register site *s;
  int a, b;
  double re, im, det_action = 0.0;

  compute_plaqdet();
  FORALLSITES(i, s) {
    for (a = XUP; a < NUMLINK; a++) {
      for (b = a + 1; b < NUMLINK; b++) {
        re = plaqdet[a][b][i].real;
        im = plaqdet[a][b][i].imag;
        det_action += (re - 1.0) * (re - 1.0);
        det_action += im * im;
      }
    }
  }
  det_action *= kappa_u1;
  g_doublesum(&det_action);
  return det_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion contribution to the action
// Ignore the ampdeg term -- since the pseudofermion src is fixed
// throughout the trajectory, it has no effect on Delta S (checked)
// sol, however, depends on the gauge fields through the CG
double d_fermion_action(Twist_Fermion *src, Twist_Fermion **sol) {
  register int i, j;
  register site *s;
  double sum = 0;
  complex ctmp;

#ifdef DEBUG_CHECK  // Check ampdeg4 term
  FORALLSITES(i, s)
    sum += ampdeg4 * (double)magsq_TF(&(src[i]));
  g_doublesum(&sum);
  node0_printf("ampdeg|chi|^2 = %.4g\n", sum);
#endif

  FORALLSITES(i, s) {
    for (j = 0; j < Norder; j++) {
      ctmp = TF_dot(&(src[i]), &(sol[j][i]));   // src^dag.sol[j]
      sum += (double)(amp4[j] * ctmp.real);
    }
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print out zeros if fermion and determinant actions not included
double d_action(Twist_Fermion *src, Twist_Fermion **sol) {
  double g_act, bmass_act, h_act, f_act = 0.0, det_act = 0.0;
  double total;
  g_act = d_gauge_action(NODET);
  bmass_act = d_bmass_action();
  h_act = d_hmom_action();
  if (kappa_u1 > IMAG_TOL)
    det_act = d_det_action();

#ifndef PUREGAUGE
  f_act = d_fermion_action(src, sol);
#endif
  node0_printf("action: ");
  node0_printf("gauge %.8g bmass %.8g ", g_act, bmass_act);
  node0_printf("det %.8g ", det_act);
  node0_printf("fermion %.8g hmom %.8g ", f_act, h_act);
  total = g_act + bmass_act + det_act + f_act + h_act;
  node0_printf("sum %.8g\n", total);
  return total;
}
// -----------------------------------------------------------------
