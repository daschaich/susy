// -----------------------------------------------------------------
// Set up generator matrices, normalization Tr(TaTb) = -delta_ab
#include "susy_includes.h"

void setup_lambda() {
  int i, j, k, l, count;
  complex inv_sqrt = cmplx(1.0 / sqrt(2.0), 0.0);
  complex i_inv_sqrt = cmplx(0.0, 1.0 / sqrt(2.0));

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

#ifdef DEBUG_CHECK
  complex tc;
  matrix tmat;

  // Print Lambdas
  node0_printf("Computing generators for SU(%d)\n", NCOL);
  for (i = 0; i < DIMF; i++){
    node0_printf("Lambda[%d]\n",i);
    if (this_node == 0)
      dumpmat(&(Lambda[i]));
  }

  // Test group theory (useful reference: arXiv:1310.5353)
#if 0   // This test seems to be specific to U(N) rather than SU(N)
  int a;
  node0_printf("Check group theory ");
  node0_printf("Sum_a Lambda^a_{kl} Lambda^a_{ij} = -delta_kj delta_il\n");
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        for (l = 0; l < NCOL; l++) {
          tc.real = Lambda[0].e[k][l].real * Lambda[0].e[i][j].real
                  - Lambda[0].e[k][l].imag * Lambda[0].e[i][j].imag;
          tc.imag = Lambda[0].e[k][l].imag * Lambda[0].e[i][j].real
                  + Lambda[0].e[k][l].real * Lambda[0].e[i][j].imag;
          for (a = 1; a < DIMF; a++) {
            tc.real += Lambda[a].e[k][l].real * Lambda[a].e[i][j].real
                     - Lambda[a].e[k][l].imag * Lambda[a].e[i][j].imag;
            tc.imag += Lambda[a].e[k][l].imag * Lambda[a].e[i][j].real
                     + Lambda[a].e[k][l].real * Lambda[a].e[i][j].imag;
          }
          if (cabs_sq(&tc) > IMAG_TOL)
            node0_printf("Sum_a La^a_{%d%d} La^a_{%d%d} = (%.4g, %.4g)\n",
                         k, l, i, j, tc.real, tc.imag);
        }
      }
    }
  }
#endif

  // Test orthogonality of products of Lambdas
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      mult_nn(&(Lambda[i]), &(Lambda[j]), &tmat);
      tc = trace(&tmat);
      if (tc.real * tc.real > IMAG_TOL)
        node0_printf("Tr[T_%d T_%d] = (%.4g, %.4g)\n",
                     i, j, tc.real, tc.imag);
    }
  }
#endif
}
// -----------------------------------------------------------------
