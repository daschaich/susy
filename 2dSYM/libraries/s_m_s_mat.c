// -----------------------------------------------------------------
// Subtract result of scalar multiplication on matrix
// c <-- c - s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_dif_matrix(matrix *b, Real s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= s * b->e[i][j].real;
      c->e[i][j].imag -= s * b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
