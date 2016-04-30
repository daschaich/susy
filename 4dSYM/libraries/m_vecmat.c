// -----------------------------------------------------------------
// Irrep vector--matrix operation with no adjoints
// c <-- b * a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void mult_vec_mat(vector *b, matrix *a, vector *c) {
  register int i, j;
  register complex y;
  for (i = 0; i < DIMF; i++) {
    CMUL(a->e[0][i], b->c[0], c->c[i]);   // Initialize
    for (j = 1; j < DIMF; j++) {
      CMUL(a->e[j][i], b->c[j], y);
      CSUM(c->c[i], y);
    }
  }
}
// -----------------------------------------------------------------
