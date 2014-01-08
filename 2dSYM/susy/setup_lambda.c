// -----------------------------------------------------------------
// Set up generator matrices
#include "susy_includes.h"
#define POS1 (NCOL * (NCOL - 1) / 2)
#define POS2 (NCOL * (NCOL - 1))
#define SUNGEN (NCOL * NCOL - 1)
#define ROOT2 1.41421356237309504880168872421
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Initialize generator matrices
void initgen(int posmat[NCOL][NCOL]) {
  int i, j, a = 0;

  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      posmat[i][j] = a;
      a++;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Translate adjoint vector to fundamental matrix
void vectomat(complex vec[SUNGEN], complex mat[NCOL][NCOL],
              int posmat[NCOL][NCOL]) {

  int k, j;
  complex mult, a, ci = cmplx(0, 1);

  // Off-diagonal elements: ???
  // mat_kj = (vec[posmat_kj] - i * vec[POS1 + posmat_kj]) / 2
  // mat_jk = (vec[posmat_kj] + i * vec[POS1 + posmat_kj]) / 2
  for (k = 0; k < NCOL; k++) {
    for (j = k + 1; j < NCOL; j++) {
      CMUL(ci, vec[POS1 + posmat[k][j]], a);
      CSUB(vec[posmat[k][j]], a, mat[k][j]);
      CMULREAL(mat[k][j], 0.5, mat[k][j]);
      CADD(vec[posmat[k][j]], a, mat[j][k]);
      CMULREAL(mat[j][k], 0.5, mat[j][k]);
    }
  }

  // Diagonal elements: ???
  // mult = vec[POS2 + i]*(1/sqrt(2+2/(1.+i))/(1.+i)) ;
  // mat_i = sum{vec[POS2 + k] / ((1 + k)* sqrt(2 + 2 / (1 + k)))};
  // mat_N = -(1 + N) * mult;
  mat[0][0] = cmplx(0, 0);
  for (k = 0; k < NCOL - 1; k++) {
    CMULREAL(vec[POS2 + k], (1.0 / sqrt(2. + 2. / (1. + k)) / (1. + k)), mult);
    for (j = 0; j < k + 1; j++)
      CSUM(mat[j][j], mult);

    a = cmplx(-1 - k, 0);
    CMUL(a, mult, mat[k + 1][k + 1]);   // Initialize next element
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute generator matrices
void computegen(complex genmat[SUNGEN][NCOL][NCOL], int posmat[NCOL][NCOL]) {
  int a, i, j;
  complex adj[SUNGEN], fund[NCOL][NCOL];

  for (a = 0; a < SUNGEN; a++) {
    for (i = 0; i < SUNGEN; i++)
      adj[i] = cmplx(0, 0);

    adj[a] = cmplx(1, 0);
    vectomat(adj, fund, posmat);
    for (i = 0; i < NCOL; i++) {
      for (j = 0; j < NCOL; j++)
        genmat[a][i][j] = fund[i][j];
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void my_gen() {
  int a, i, j, posmat[NCOL][NCOL];
  complex isr2 = cmplx(0, sqrt(2));
  complex genmat[SUNGEN][NCOL][NCOL];

  initgen(posmat);
  computegen(genmat, posmat);

  // Map to data structures in code
  for (a = 0; a < SUNGEN; a++) {
    for (i = 0; i < NCOL; i++) {
      for (j = 0; j < NCOL; j++)    // Lambda = i sqrt(2) genmat
        CMUL(isr2, genmat[a][i][j], Lambda[a].e[i][j]);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Normalization: Tr(TaTb) = -delta_ab
void setup_lambda() {
  int i, j;

#ifdef DEBUG_CHECK
  int k, l, a;
  complex trace, tt;
  node0_printf("Computing generators for SU(N)\n");
#endif

  my_gen();
  if (NUMGEN == NCOL * NCOL) {
    clear_su3mat_f(&(Lambda[NUMGEN - 1]));
    for (i = 0; i < NCOL; i++)
      Lambda[NUMGEN - 1].e[i][i] = cmplx(0, 1.0 / sqrt(NCOL));
  }

#ifdef DEBUG_CHECK
  // Print Lambdas
  for (i = 0; i < NUMGEN; i++){
    node0_printf("Lambda[%d]\n",i);
    if (this_node == 0)
      dumpmat_f(&(Lambda[i]));
  }

  // Test group theory
  node0_printf("Check group theory ");
  node0_printf("Sum_a Lambda^a_{kl} Lambda^a_{ij} = -delta_kj delta_il\n");
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        for (l = 0; l < NCOL; l++) {
          trace = cmplx(0, 0);
          for (a = 0; a < NUMGEN; a++) {
            CMUL(Lambda[a].e[k][l], Lambda[a].e[i][j], tt);
            CSUM(trace, tt);
          }
          if (cabs_sq(&trace) > 1e-8)
            node0_printf("Sum_a Lambda^a_{%d%d} Lambda^a_{%d%d} = (%.4g, %.4g)\n",
                         k, j, i, l, trace.real, trace.imag);
        }
      }
    }
  }
#endif

  // Test orthogonality and compute products of Lambdas for fermion forces
  for (i = 0; i < NUMGEN; i++) {
    for (j = 0; j < NUMGEN; j++) {
      mult_su3_nn_f(&(Lambda[i]), &(Lambda[j]), &(Lambda_prod[i][j]));
#ifdef DEBUG_CHECK
      trace = trace_su3_f(&(Lambda_prod[i][j]));
      if (trace.real * trace.real > 1e-6)
        node0_printf("Tr[T_%d T_%d] = (%.4g, %.4g)\n", i, j, trace.real, trace.imag);
#endif
    }
  }
}
// -----------------------------------------------------------------
