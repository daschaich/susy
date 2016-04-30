// -----------------------------------------------------------------
// Subtract result of complex scalar multiplication on irrep matrix
// c <-- a - s*b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void c_scalar_mult_sub_mat(matrix *a, matrix *b, complex *s, matrix *c) {
  register int i, j;
  complex t;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      CMUL(b->e[i][j], *s, t);
      CSUB(a->e[i][j], t, c->e[i][j]);
    }
  }
}
// -----------------------------------------------------------------
