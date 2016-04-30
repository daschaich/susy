// -----------------------------------------------------------------
// Adjoint of a fundamental matrix
// b <-- adag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void adjoint_f(matrix_f *a, matrix_f *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      CONJG(a->e[j][i], b->e[i][j]);
  }
}
// -----------------------------------------------------------------
