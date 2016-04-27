// -----------------------------------------------------------------
// Add result of complex scalar multiplication on fundamental matrix
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_add_su3mat_f(su3_matrix_f *a, su3_matrix_f *b,
                                complex *s, su3_matrix_f *c) {

  register int i, j;
  complex t;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      CMUL(b->e[i][j], *s, t);
      CADD(a->e[i][j], t, c->e[i][j]);
    }
  }
}
// -----------------------------------------------------------------
