// -----------------------------------------------------------------
// Add result of scalar multiplication on adjoint matrix
// c <-- c + s * bdag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_sum_adj_matrix(matrix *b, Real s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += s * b->e[j][i].real;
      c->e[i][j].imag -= s * b->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------
