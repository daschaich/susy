// -----------------------------------------------------------------
// Subtract result of irrep vector--matrix operation with no adjoints
// c <-- c - b * a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void mult_vec_mat_dif(vector *b, matrix *a, vector *c) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      c->c[i].real -= a->e[j][i].real * b->c[j].real
                    - a->e[j][i].imag * b->c[j].imag;
      c->c[i].imag -= a->e[j][i].imag * b->c[j].real
                    + a->e[j][i].real * b->c[j].imag;
    }
  }
}
// -----------------------------------------------------------------
