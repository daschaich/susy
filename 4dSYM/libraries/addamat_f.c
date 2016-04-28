// -----------------------------------------------------------------
// Add two fundamental matrices with one adjoint
// c <-- a + bdag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void add_su3_adj_matrix_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real + b->e[j][i].real;
      c->e[i][j].imag = a->e[i][j].imag - b->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------
