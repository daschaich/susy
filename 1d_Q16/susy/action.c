// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the CG should already have been run,
// so that the vector **sol contains (M_adjoint*M+shift[n])^(-1) * src
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Bosonic contribution to the action
// Uses some pointer arithmetic that might be susceptible to bugs
double bosonic_action(double * so3_sq, double * so6_sq, double * Myers) {
  register int i,l;
  register site *s;
  double b_action;
  int j, k;
  matrix tmat, tmat2;
  msg_tag *tag;
  *so3_sq = 0.0;
  *so6_sq = 0.0;
  *Myers = 0.0;
  
  
  // Scalar kinetic term
  tag = start_gather_site(F_OFFSET(X[0]), sizeof(matrix) * NSCALAR,
                          TUP, EVENANDODD, gen_pt[0]);
  FORALLSITES(i, s) {
    for (j = 0; j < 3; j++)
      *so3_sq -= realtrace_nn(&s->X[j], &s->X[j]);
    for (j = 3; j < NSCALAR; j++)
      *so6_sq -= realtrace_nn(&s->X[j], &s->X[j]);
  }
  b_action = *so3_sq + *so6_sq;

  // Scalar hopping term [D_t^+ X_i(t)]^2
  wait_gather(tag);
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      mult_nn(&(s->link), (matrix *)(gen_pt[0][i] + j), &tmat);
      mult_nn(&(s->X[j]), &tmat, &tmat2);
      b_action += realtrace(&(s->link), &tmat2);
    }
  }
  b_action *= 2.0;
  cleanup_gather(tag);

  // Commutator term
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      for (k = j + 1; k < NSCALAR; k++) {
        mult_nn(&(s->X[j]), &(s->X[k]), &tmat);
        mult_nn_dif(&(s->X[k]), &(s->X[j]), &tmat);
        b_action -= realtrace_nn(&tmat, &tmat);
      }
    }
  }

  // Scalar potential terms
  // Couplings are set differently depending on BMN/BFSS in setup.c
  *so3_sq *= mass_so3;
  *so6_sq *= mass_so6;
  b_action += *so3_sq + *so6_sq;

  
#ifdef BMN
  FORALLSITES(i, s) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        if (j == k)
          continue;
        for (l = 0; l < 3; l++) {
          if((j == l) || (k == l))
            continue;

          mult_nn(&(s->X[k]), &(s->X[l]), &tmat);
          *Myers -= epsilon[j][k][l] * realtrace_nn(&s->X[j], &tmat);
        }
      }
    }
  }
  *Myers *= mass_Myers;
  
  b_action += *Myers;

#endif

  b_action *= kappa;
  
  *so3_sq *= kappa;
  *so6_sq *= kappa;
  *Myers *= kappa;
  
  g_doublesum(&b_action);
  
  g_doublesum(so3_sq);
  g_doublesum(so6_sq);
  g_doublesum(Myers);
  
  return b_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion contribution to the action
// Include the ampdeg term to allow sanity check that the fermion action
// is NFERMION DIMF volume on average
// Since the pseudofermion src is fixed throughout the trajectory,
// ampdeg actually has no effect on Delta S (checked)
// sol, however, depends on the gauge fields through the CG
double fermion_action(matrix *src[NFERMION], matrix **sol[NFERMION]) {
  register int i, j, k;
  register site *s;
  double sum = 0.0;
  complex ctmp;
#ifdef DEBUG_CHECK
  double im = 0.0;
#endif

  FORALLSITES(i, s) {
    for (k = 0; k < NFERMION; k++) {
      sum += ampdeg4 * (double)realtrace(&(src[k][i]), &(src[k][i]));
      for (j = 0; j < Norder; j++) {
        ctmp = complextrace_an(&(src[k][i]), &(sol[k][j][i]));   // src^dag.sol[j]
        sum += (double)(amp4[j] * ctmp.real);
#ifdef DEBUG_CHECK  // Make sure imaginary part vanishes
        im += (double)(amp4[j] * ctmp.imag);
#endif
      }
    }
  }
  g_doublesum(&sum);
#ifdef DEBUG_CHECK  // Make sure imaginary part vanishes
  g_doublesum(&im);
  node0_printf("S_f = (%.4g, %.4g)\n", sum, im);
#endif
  return sum;
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Gauge and scalar momenta contribution to the action
double mom_action() {
  register int i, j;
  register site *s;
  double sum = 0.0;

  FORALLSITES(i, s) {
    sum += (double)realtrace(&(s->mom), &(s->mom));
    for (j = 0; j < NSCALAR; j++)
      sum += (double)realtrace(&(s->mom_X[j]), &(s->mom_X[j]));
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Print out zeros for pieces of the action that aren't included
double action(matrix **src[NFERMION], matrix ***sol[NFERMION]) {
  double b_act, p_act, so3_act, so6_act, Myers_act;
  double total;

  b_act = bosonic_action(&so3_act, &so6_act, &Myers_act);
  node0_printf("action: so3 %.8g so6 %.8g Myers %.8g ",
               so3_act, so6_act, Myers_act);
  node0_printf("boson %.8g ", b_act);
  total = b_act;

#ifndef PUREGAUGE
  int n, m;
  double f_act;
  for (n = 0; n < Nroot; n++) {
    f_act = 0.0;
    for (m = 0; m < NFERMION; m++)
      f_act += fermion_action(src[m][n], sol[m][n]);

    node0_printf("fermion%d %.8g ", n, f_act);
    total += f_act;
  }
#endif

  p_act = mom_action();
  node0_printf("mom %.8g ", p_act);
  total += p_act;
  node0_printf("sum %.8g\n", total);
  return total;
}
// -----------------------------------------------------------------
