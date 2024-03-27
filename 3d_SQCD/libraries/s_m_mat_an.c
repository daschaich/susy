// -----------------------------------------------------------------
// Scaled matrix multiplication with adjoint of first matrix
// c <-- s * c + adag * b
// c <-- s * c - adag * b
// c <-- s * adag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_an_sum(matrix *a, matrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real += s * (a->e[k][i].real * b->e[k][j].real
                              + a->e[k][i].imag * b->e[k][j].imag);
        c->e[i][j].imag += s * (a->e[k][i].real * b->e[k][j].imag
                              - a->e[k][i].imag * b->e[k][j].real);
      }
    }
  }
}

void scalar_mult_an_dif(matrix *a, matrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real -= s * (a->e[k][i].real * b->e[k][j].real
                              + a->e[k][i].imag * b->e[k][j].imag);
        c->e[i][j].imag -= s * (a->e[k][i].real * b->e[k][j].imag
                              - a->e[k][i].imag * b->e[k][j].real);
      }
    }
  }
}

void scalar_mult_an(matrix *a, matrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = s * (a->e[0][i].real * b->e[0][j].real
                           + a->e[0][i].imag * b->e[0][j].imag);
      c->e[i][j].imag = s * (a->e[0][i].real * b->e[0][j].imag
                           - a->e[0][i].imag * b->e[0][j].real);
      for (k = 1; k < NCOL; k++) {
        c->e[i][j].real += s * (a->e[k][i].real * b->e[k][j].real
                              + a->e[k][i].imag * b->e[k][j].imag);
        c->e[i][j].imag += s * (a->e[k][i].real * b->e[k][j].imag
                              - a->e[k][i].imag * b->e[k][j].real);
      }
    }
  }
}
// -----------------------------------------------------------------
