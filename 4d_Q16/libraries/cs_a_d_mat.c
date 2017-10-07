// -----------------------------------------------------------------
// Add complex scalar to diagonal components of matrix
// a <-- a + cI
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void c_scalar_add_diag(matrix *a, complex *c) {
  register int i;

  for (i = 0; i < NCOL; i++) {
    a->e[i][i].real += c->real;
    a->e[i][i].imag += c->imag;
  }
}
// -----------------------------------------------------------------
