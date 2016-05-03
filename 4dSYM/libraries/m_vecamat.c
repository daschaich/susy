// -----------------------------------------------------------------
// Irrep vector--matrix operation with adjoint matrix
// c <-- b * adag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void mult_vec_adj_mat(vector *b, matrix *a, vector *c) {
  register int i, j;
  register complex y;
  for (i = 0; i < DIMF; i++) {
    CMULJ_(a->e[i][0], b->c[0], c->c[i]);   // Initialize
    for (j = 1; j < DIMF; j++) {
      CMULJ_(a->e[i][j], b->c[j], y);
      CSUM(c->c[i], y);
    }
  }
}
// -----------------------------------------------------------------
