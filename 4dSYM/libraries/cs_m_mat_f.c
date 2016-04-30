// -----------------------------------------------------------------
// Complex scalar multiplication on fundamental matrix
// c <-- s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void c_scalar_mult_mat_f(matrix_f *b, complex *s, matrix_f *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      CMUL(b->e[i][j], *s, c->e[i][j]);
  }
}
// -----------------------------------------------------------------
