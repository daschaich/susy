// -----------------------------------------------------------------
// Add two fundamental matrices
// c <-- a + b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void add_matrix_f(matrix_f *a, matrix_f *b, matrix_f *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      CADD(a->e[i][j], b->e[i][j], c->e[i][j]);
  }
}
// -----------------------------------------------------------------
