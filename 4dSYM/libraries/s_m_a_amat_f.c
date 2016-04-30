// -----------------------------------------------------------------
// Add result of scalar multiplication on adjoint fundamental matrix
// c <-- a + s*b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_add_adj_matrix_f(matrix_f *a, matrix_f *b,
                                  Real s, matrix_f *c) {

  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real + s * b->e[j][i].real;
      c->e[i][j].imag = a->e[i][j].imag - s * b->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------
