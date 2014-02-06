// -----------------------------------------------------------------
// Add result of scalar multiplication on fundamental matrix
// c <-- a + s*b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_add_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b,
                                  Real s, su3_matrix_f *c) {

  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real + s * b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + s * b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
