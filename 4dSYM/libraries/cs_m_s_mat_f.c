// -----------------------------------------------------------------
// Subtract result of complex scalar multiplication on fundamental matrix
// c <-- a - s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_sub_su3mat_f(su3_matrix_f *a, su3_matrix_f *b,
                                complex *s, su3_matrix_f *c) {

  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real - b->e[i][j].real * s->real \
                                        + b->e[i][j].imag * s->imag;
      c->e[i][j].imag = a->e[i][j].imag - b->e[i][j].imag * s->real \
                                        - b->e[i][j].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------
