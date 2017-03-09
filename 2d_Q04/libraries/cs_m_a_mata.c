// -----------------------------------------------------------------
// Add result of complex scalar multiplication on adjoint matrix
// c <-- c + s * bdag
// Just flip i<-->j and negate b.imag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void c_scalar_mult_sum_mat_adj(matrix *b, complex *s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += b->e[j][i].real * s->real + b->e[j][i].imag * s->imag;
      c->e[i][j].imag -= b->e[j][i].imag * s->real - b->e[j][i].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------
