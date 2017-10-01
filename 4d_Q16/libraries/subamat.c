// -----------------------------------------------------------------
// Subtract two matrices using the adjoint of the second
// c <-- a - bdag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void sub_adj_matrix(matrix *a, matrix *b, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real - b->e[j][i].real;
      c->e[i][j].imag = a->e[i][j].imag + b->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------
