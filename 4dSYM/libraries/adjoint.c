// -----------------------------------------------------------------
// Adjoint of a matrix
// b <-- (+/-)adag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void adjoint(matrix *a, matrix *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      CONJG(a->e[j][i], b->e[i][j]);
  }
}

void neg_adjoint(matrix *a, matrix *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      b->e[i][j].real = -a->e[j][i].real;
      b->e[i][j].imag = a->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------
