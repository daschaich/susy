// -----------------------------------------------------------------
// Subtract result of complex scalar multiplication on irrep matrix
// c <-- a - s*b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_sub_su3mat(su3_matrix *a, su3_matrix *b,
                              complex *s, su3_matrix *c) {

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
