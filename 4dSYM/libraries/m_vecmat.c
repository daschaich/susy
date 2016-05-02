// -----------------------------------------------------------------
// Irrep vector--matrix operation with no adjoints
// c <-- b * a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void mult_vec_mat(vector *b, matrix *a, vector *c) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    // Initialize
    c->c[i].real = a->e[0][i].real * b->c[0].real
                 - a->e[0][i].imag * b->c[0].imag;
    c->c[i].imag = a->e[0][i].imag * b->c[0].real
                 + a->e[0][i].real * b->c[0].imag;
    for (j = 1; j < DIMF; j++) {
      c->c[i].real += a->e[j][i].real * b->c[j].real
                    - a->e[j][i].imag * b->c[j].imag;
      c->c[i].imag += a->e[j][i].imag * b->c[j].real
                    + a->e[j][i].real * b->c[j].imag;
    }
  }
}
// -----------------------------------------------------------------
