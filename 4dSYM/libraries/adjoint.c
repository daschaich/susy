// -----------------------------------------------------------------
// Adjoint of an irrep matrix
// b <-- adag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void adjoint(matrix *a, matrix *b) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++)
      CONJG(a->e[j][i], b->e[i][j]);
  }
}
// -----------------------------------------------------------------
