// -----------------------------------------------------------------
// Subtract result of complex scalar multiplication on fundamental matrix
// c <-- c - s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void c_scalar_mult_dif_mat_f(matrix_f *b, complex *s, matrix_f *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= b->e[i][j].real * s->real - b->e[i][j].imag * s->imag;
      c->e[i][j].imag -= b->e[i][j].imag * s->real + b->e[i][j].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------
