// -----------------------------------------------------------------
// Subtract two matrices -- unlike CDIF, output is always last
// c <-- c - b
// c <-- a - b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void dif_matrix(matrix *b, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= b->e[i][j].real;
      c->e[i][j].imag -= b->e[i][j].imag;
    }
  }
}

void sub_matrix(matrix *a, matrix *b, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real - b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag - b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
