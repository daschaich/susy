// -----------------------------------------------------------------
// Subtract two irrep matrices
// c <-- a - b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void sub_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++)
      CSUB(a->e[i][j], b->e[i][j], c->e[i][j]);
  }
}
// -----------------------------------------------------------------
