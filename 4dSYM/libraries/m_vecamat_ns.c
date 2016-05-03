// -----------------------------------------------------------------
// Subtract result of irrep vector--matrix operation with adjoint matrix
// c <-- c - b * adag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void mult_vec_adj_mat_dif(vector *b, matrix *a, vector *c) {
  register int i, j;
  register complex y;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      CMULJ_(a->e[i][j], b->c[j], y);
      CDIF(c->c[i], y);
    }
  }
}
// -----------------------------------------------------------------
