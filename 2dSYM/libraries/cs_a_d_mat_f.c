// -----------------------------------------------------------------
// Add complex scalar to diagonal components of fundamental matrix
// a <-- a + cI
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_add_diag_su3_f(su3_matrix_f *a, complex *c) {
  register int i;

  for (i = 0; i < NCOL; i++) {
    a->e[i][i].real += c->real;
    a->e[i][i].imag += c->imag;
  }
}
// -----------------------------------------------------------------
