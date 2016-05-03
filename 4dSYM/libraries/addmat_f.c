// -----------------------------------------------------------------
// Add two fundamental matrices
// c <-- c + b
// c <-- a + b
// NB: Output is always last, unlike CSUM
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void sum_matrix_f(matrix_f *b, matrix_f *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += b->e[i][j].real;
      c->e[i][j].imag += b->e[i][j].imag;
    }
  }
}

void add_matrix_f(matrix_f *a, matrix_f *b, matrix_f *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real + b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
