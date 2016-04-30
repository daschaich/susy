// -----------------------------------------------------------------
// Add result of complex scalar multiplication on adjoint fundamental matrix
// c <-- a + s * bdag
// Just flip i<-->j and negate b.imag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void c_scalar_mult_add_mat_adj_f(matrix_f *a, matrix_f *b,
                                 complex *s, matrix_f *c) {

  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real + b->e[j][i].real * s->real \
                                        + b->e[j][i].imag * s->imag;
      c->e[i][j].imag = a->e[i][j].imag - b->e[j][i].imag * s->real \
                                        + b->e[j][i].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------
