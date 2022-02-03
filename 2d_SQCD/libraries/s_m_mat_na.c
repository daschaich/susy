// -----------------------------------------------------------------
// Scaled matrix multiplication with adjoint of second matrix
// c <-- c + a * bdag
// c <-- c - a * bdag
// c <-- a * bdag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_na_sum(matrix *a, matrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real += s * (a->e[i][k].real * b->e[j][k].real
                              + a->e[i][k].imag * b->e[j][k].imag);
        c->e[i][j].imag += s * (a->e[i][k].imag * b->e[j][k].real
                              - a->e[i][k].real * b->e[j][k].imag);
      }
    }
  }
}

void scalar_mult_na_dif(matrix *a, matrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real -= s * (a->e[i][k].real * b->e[j][k].real
                              + a->e[i][k].imag * b->e[j][k].imag);
        c->e[i][j].imag -= s * (a->e[i][k].imag * b->e[j][k].real
                              - a->e[i][k].real * b->e[j][k].imag);
      }
    }
  }
}

void scalar_mult_na(matrix *a, matrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = s * (a->e[i][0].real * b->e[j][0].real
                           + a->e[i][0].imag * b->e[j][0].imag);
      c->e[i][j].imag = s * (a->e[i][0].imag * b->e[j][0].real
                           - a->e[i][0].real * b->e[j][0].imag);
      for (k = 1; k < NCOL; k++) {
        c->e[i][j].real += s * (a->e[i][k].real * b->e[j][k].real
                              + a->e[i][k].imag * b->e[j][k].imag);
        c->e[i][j].imag += s * (a->e[i][k].imag * b->e[j][k].real
                              - a->e[i][k].real * b->e[j][k].imag);
      }
    }
  }
}
// -----------------------------------------------------------------
