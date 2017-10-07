// -----------------------------------------------------------------
// Scalar multiplication on matrix
// b <-- s * a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_matrix(matrix *a, Real s, matrix *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      b->e[i][j].real = s * a->e[i][j].real;
      b->e[i][j].imag = s * a->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
