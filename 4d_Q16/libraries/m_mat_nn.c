// -----------------------------------------------------------------
// Matrix multiplication with no adjoints
// c <-- c + a * b
// c <-- c - a * b
// c <-- a * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void mult_nn_sum(matrix *a, matrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[k][j].real
                         - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[k][j].real
                         + a->e[i][k].real * b->e[k][j].imag;
      }
    }
  }
}

void mult_nn_dif(matrix *a, matrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real -= a->e[i][k].real * b->e[k][j].real
                         - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag -= a->e[i][k].imag * b->e[k][j].real
                         + a->e[i][k].real * b->e[k][j].imag;
      }
    }
  }
}

void mult_nn(matrix *a, matrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = a->e[i][0].real * b->e[0][j].real
                      - a->e[i][0].imag * b->e[0][j].imag;
      c->e[i][j].imag = a->e[i][0].imag * b->e[0][j].real
                      + a->e[i][0].real * b->e[0][j].imag;
      for (k = 1; k < NCOL; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[k][j].real
                         - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[k][j].real
                         + a->e[i][k].real * b->e[k][j].imag;
      }
    }
  }
}
// -----------------------------------------------------------------
