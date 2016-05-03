// -----------------------------------------------------------------
// Scalar multiplication on adjoint of fundamental matrix
// b <-- s * adag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_adj_matrix_f(matrix_f *a, Real s, matrix_f *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      b->e[i][j].real = s * a->e[j][i].real;
      b->e[i][j].imag = -1.0 * s * a->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------
