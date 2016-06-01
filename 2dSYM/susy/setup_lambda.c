// -----------------------------------------------------------------
// Set up generator matrices, normalization Tr(TaTb) = -delta_ab
#include "susy_includes.h"

void setup_lambda() {
  int i, j, k, l, count;
  complex inv_sqrt = cmplx(1.0 / sqrt(2.0), 0.0);
  complex i_inv_sqrt = cmplx(0.0, 1.0 / sqrt(2.0));

#ifdef DEBUG_CHECK
  int a;
  complex trace, tt;
  node0_printf("Computing generators for U(N)\n");
#endif

  // Make sure Lambda matrices are initialized
  for (i = 0; i < DIMF; i++)
    clear_mat(&(Lambda[i]));

  // N * (N - 1) off-diagonal SU(N) generators
  // (T^{ij, +})_{kl} = i * (de_{ki} de_{lj} + de_{kj} de_{li}) / sqrt(2)
  // (T^{ij, -})_{kl} = (de_{ki} de_{lj} - de_{kj} de_{ki}) / sqrt(2)
  // Sign in second chosen to match previous values
  count = 0;
  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        for (l = 0; l < NCOL; l++) {
          if (k == i && l == j) {
            CSUM(Lambda[count].e[k][l], i_inv_sqrt);
            CSUM(Lambda[count + 1].e[k][l], inv_sqrt);
          }
          else if (k == j && l == i) {
            CSUM(Lambda[count].e[k][l], i_inv_sqrt);
            CDIF(Lambda[count + 1].e[k][l], inv_sqrt);
          }
        }
      }
      count += 2;
    }
  }
  if (count != NCOL * (NCOL - 1)) {
    node0_printf("ERROR: Wrong number of off-diagonal generators, ");
    node0_printf("%d vs. %d\n", count, NCOL * (NCOL - 1));
    terminate(1);
  }

  // N - 1 diagonal SU(N) generators
  // T^k = i * diag(1, 1, ..., -k, 0, ..., 0) / sqrt(k * (k + 1))
  for (i = 0; i < NCOL - 1; i++) {
    j = NCOL * (NCOL - 1) + i;    // Index after +/- above
    k = i + 1;
    i_inv_sqrt = cmplx(0.0, 1.0 / sqrt(k * (k + 1.0)));
    for (l = 0; l <= k; l++)
      Lambda[j].e[l][l] = i_inv_sqrt;
    CMULREAL(Lambda[j].e[k][k], -1.0 * k, Lambda[j].e[k][k]);
  }

  // U(1) generator i * I_N / sqrt(N)
  if (DIMF == NCOL * NCOL) {    // Allow SU(N) compilation for now
    i_inv_sqrt = cmplx(0.0, sqrt(one_ov_N));
    clear_mat(&(Lambda[DIMF - 1]));
    for (i = 0; i < NCOL; i++)
      Lambda[DIMF - 1].e[i][i] = i_inv_sqrt;
  }

#ifdef DEBUG_CHECK
  // Print Lambdas
  for (i = 0; i < DIMF; i++){
    node0_printf("Lambda[%d]\n",i);
    if (this_node == 0)
      dumpmat(&(Lambda[i]));
  }

  // Test group theory
  node0_printf("Check group theory ");
  node0_printf("Sum_a Lambda^a_{kl} Lambda^a_{ij} = -delta_kj delta_il\n");
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        for (l = 0; l < NCOL; l++) {
          trace = cmplx(0, 0);
          for (a = 0; a < DIMF; a++) {
            CMUL(Lambda[a].e[k][l], Lambda[a].e[i][j], tt);
            CSUM(trace, tt);
          }
          if (cabs_sq(&trace) > IMAG_TOL)
            node0_printf("Sum_a La^a_{%d%d} La^a_{%d%d} = (%.4g, %.4g)\n",
                         k, j, i, l, trace.real, trace.imag);
        }
      }
    }
  }
#endif

  // Test orthogonality and compute products of Lambdas for fermion forces
#ifdef DEBUG_CHECK
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      mult_nn(&(Lambda[i]), &(Lambda[j]), &tmat);
      trace = trace(&tmat);
      if (trace.real * trace.real > IMAG_TOL)
        node0_printf("Tr[T_%d T_%d] = (%.4g, %.4g)\n",
                     i, j, trace.real, trace.imag);
    }
  }
#endif
}
// -----------------------------------------------------------------
