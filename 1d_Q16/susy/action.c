// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the CG should already have been run,
// so that the vector **sol contains (M_adjoint*M+shift[n])^(-1) * src
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Bosonic contribution to the action
double bosonic_action(double *so3_sq, double *so6_sq, double *comm, double *Myers, double *Fadeev) {

  register int i,l;
  register site *s;
  int j, k;
  double b_action = 0.0;
  matrix tmat;
#ifndef UNGAUGED
  matrix tmat2;
#endif
  msg_tag *tag[NSCALAR];

  // Initialize
  *so3_sq = 0.0;
  *so6_sq = 0.0;
  *comm = 0.0;
  *Myers = 0.0;
  *Fadeev = 0.0;

  // Scalar kinetic term -Tr[D_t X(t)]^2
  //   -Tr[U(t) X(t+1) Udag(t) - X(t)]^2
  //     = Tr[2 X(t) U(t) X(t+1) Udag(t) - X(t+1) X(t+1) - X(t) X(t)]
  // Sum over t --> 2 Tr[Udag(t) X(t) U(t) X(t+1) - X(t) X(t)]
  for (j = 0; j < NSCALAR; j++) {
    tag[j] = start_gather_site(F_OFFSET(X[j]), sizeof(matrix),
                               TUP, EVENANDODD, gen_pt[j]);
  }

  // On-site piece of scalar kinetic term
  // (Has same form as some scalar potential terms, so will re-use below)
  FORALLSITES(i, s) {
    for (j = 0; j < 3; j++){
      *so3_sq -= (double)realtrace_nn(&(s->X[j]), &(s->X[j]));
    }
    for (j = 3; j < NSCALAR; j++){
      *so6_sq -= (double)realtrace_nn(&(s->X[j]), &(s->X[j]));    
    }
  }
  b_action = *so3_sq + *so6_sq;

  // Nearest-neighbor piece of scalar kinetic term
  for (j = 0; j < NSCALAR; j++) {
    wait_gather(tag[j]);
    FORALLSITES(i, s) {
#ifndef UNGAUGED
      mult_nn(&(s->link), (matrix *)(gen_pt[j][i]), &tmat);
      mult_na(&tmat, &(s->link), &tmat2);
      b_action += (double)realtrace_nn(&(s->X[j]), &tmat2);
#else
      b_action += (double)realtrace_nn(&(s->X[j]), (matrix *)(gen_pt[j][i]));
#endif
    }
    cleanup_gather(tag[j]);
  }
  b_action *= 2.0;

  // Commutator term
  //   sum_{i<j} -Tr[X_i, X_j]^2 = sum_{i<j} -Tr[X_i X_j - X_j X_i]^2
  //     = sum_{i<j} -Tr[  X_i X_j X_i X_j - X_i X_j X_j X_i
  //                     - X_j X_i X_i X_j + X_j X_i X_j X_i]
  //     = sum_{i<j} -2Tr[X_i X_j X_i X_j - X_i X_j X_j X_i]
  //     = sum_{i != j} -Tr[X_i X_j X_i X_j - X_i X_j X_j X_i]
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      for (k = j + 1; k < NSCALAR; k++) {
        mult_nn(&(s->X[j]), &(s->X[k]), &tmat);
        mult_nn_dif(&(s->X[k]), &(s->X[j]), &tmat);
        *comm -= realtrace_nn(&tmat, &tmat);
      }
    }
  }
  *comm *= 0.5;

  // Scalar potential terms
  // Couplings are set differently in setup.c depending on BMN vs. BFSS
  *so3_sq *= mass_so3;
  *so6_sq *= mass_so6;
  b_action += *so3_sq + *so6_sq + *comm;

#ifdef BMN
  // Myers term -Tr[epsilon_{jkl} X_j(t) X_k(t) X_l(t)]]
  FORALLSITES(i, s) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        if (j == k)
          continue;
        for (l = 0; l < 3; l++) {
          if ((j == l) || (k == l))
            continue;

          // Different sign to follow Costa et al
          // Ref to notes dated Jun 02, 2019
          mult_nn(&(s->X[k]), &(s->X[l]), &tmat);
          if (epsilon[j][k][l] > 0)
            *Myers += realtrace_nn(&s->X[j], &tmat);
          else if (epsilon[j][k][l] < 0)
            *Myers -= realtrace_nn(&s->X[j], &tmat);
        }
      }
    }
  }
  *Myers *= mass_Myers;
  b_action += *Myers;
#endif

// Fadeev Popov Term
FORALLSITES(i,s) {
  for (j = 0; j < NCOL; j++){
    for (k = 0; k < NCOL; k++){
      if (j != k){
        *Fadeev += (double)-0.5*log((sin((&(s->link.e[j][j]) - &(s->link.e[k][k]))/2.0))*(sin((&(s->link.e[j][j]) - &(s->link.e[k][k]))/2.0)));
      }
    }
  }
}
  b_action *= kappa;
  *so3_sq *= kappa;
  *so6_sq *= kappa;
  *comm *= kappa;
  *Myers *= kappa;
  b_action += *Fadeev;
  
  g_doublesum(&b_action);
  g_doublesum(so3_sq);
  g_doublesum(so6_sq);
  g_doublesum(comm);
  g_doublesum(Myers);
  g_doublesum(Fadeev);
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
double fermion_action(matrix **src, matrix ***sol) {
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
        // src^dag.sol[j]
        ctmp = complextrace_an(&(src[k][i]), &(sol[j][k][i]));
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
// Helper routine computes agnitude squared of an anti-hermition matrix
// including the factor of 1/2 in the effective hamiltonian
Real ahmat_mag_sq(anti_hermitmat *ah) {
  register int i;
  register Real sum;

  sum = ah->im_diag[0] * ah->im_diag[0];
  for (i = 1; i < NCOL; i++)
    sum += ah->im_diag[i] * ah->im_diag[i];
  sum *= 0.5;

  for (i = 0; i < N_OFFDIAG; i++) {
    sum += ah->m[i].real * ah->m[i].real;
    sum += ah->m[i].imag * ah->m[i].imag;
  }

  return sum;
}

#ifndef UNGAUGED
double gauge_mom_action() {
  register int i;
  register site *s;
  double sum = 0.0;

  FORALLSITES(i, s)
    sum += (double)ahmat_mag_sq(&(s->mom));

  g_doublesum(&sum);
  return sum;
}
#endif

double scalar_mom_action() {
  register int i, j;
  register site *s;
  double sum = 0.0;

  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++)
      sum += (double)realtrace(&(s->mom_X[j]), &(s->mom_X[j]));
  }
  g_doublesum(&sum);
  return 0.5 * sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print out zeros for pieces of the action that aren't included
double action(matrix ***src, matrix ****sol) {
  double p_act = 0.0, so3_act, so6_act, comm_act, Myers_act, Fadeev_act, total;

  // Includes so3, so6, Myers and kinetic
  total = bosonic_action(&so3_act, &so6_act, &comm_act, &Myers_act, &Fadeev_act);
  node0_printf("action: so3 %.8g so6 %.8g comm %.8g Myers %.8g Fadeev %.8g boson %.8g ",
               so3_act, so6_act, comm_act, Myers_act, Fadeev_act, total);

#ifndef PUREGAUGE
  int n;
  double f_act;
  for (n = 0; n < Nroot; n++) {
    f_act = fermion_action(src[n], sol[n]);

    node0_printf("fermion%d %.8g ", n, f_act);
    total += f_act;
  }
#endif

#ifndef UNGAUGED
  p_act = gauge_mom_action();
#endif
  node0_printf("Umom %.8g ", p_act);
  total += p_act;
  p_act = scalar_mom_action();
  node0_printf("Xmom %.8g ", p_act);
  total += p_act;
  node0_printf("sum %.8g\n", total);
  return total;
}
// -----------------------------------------------------------------
