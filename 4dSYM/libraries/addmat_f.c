// -----------------------------------------------------------------
// Add two fundamental matrices
// c <-- a + b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void add_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      CADD(a->e[i][j], b->e[i][j], c->e[i][j]);
  }
}
// -----------------------------------------------------------------
