// -----------------------------------------------------------------
// Scalar multiplied fundamental matrix multiplication with no adjoints
// c <-- c + s * a * b
// c <-- c - s * a * b
// c <-- s * a * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_nn_sum_f(matrix_f *a, matrix_f *b, Real s, matrix_f *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real += s * (a->e[i][k].real * b->e[k][j].real
                              - a->e[i][k].imag * b->e[k][j].imag);
        c->e[i][j].imag += s * (a->e[i][k].imag * b->e[k][j].real
                                + a->e[i][k].real * b->e[k][j].imag);
      }
    }
  }
}

void scalar_mult_nn_dif_f(matrix_f *a, matrix_f *b, Real s, matrix_f *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real -= s * (a->e[i][k].real * b->e[k][j].real
                              - a->e[i][k].imag * b->e[k][j].imag);
        c->e[i][j].imag -= s * (a->e[i][k].imag * b->e[k][j].real
                                + a->e[i][k].real * b->e[k][j].imag);
      }
    }
  }
}

void scalar_mult_nn_f(matrix_f *a, matrix_f *b, Real s, matrix_f *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = s * (a->e[i][0].real * b->e[0][j].real
                           - a->e[i][0].imag * b->e[0][j].imag);
      c->e[i][j].imag = s * (a->e[i][0].imag * b->e[0][j].real
                           + a->e[i][0].real * b->e[0][j].imag);
      for (k = 1; k < NCOL; k++) {
        c->e[i][j].real += s * (a->e[i][k].real * b->e[k][j].real
                              - a->e[i][k].imag * b->e[k][j].imag);
        c->e[i][j].imag += s * (a->e[i][k].imag * b->e[k][j].real
                                + a->e[i][k].real * b->e[k][j].imag);
      }
    }
  }
}
// -----------------------------------------------------------------
