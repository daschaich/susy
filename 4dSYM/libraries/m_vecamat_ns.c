// -----------------------------------------------------------------
// Subtract result of irrep vector--matrix operation with adjoint matrix
// c <-- c - b * adag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_su3_vec_adj_mat_nsum(su3_vector *b, su3_matrix *a, su3_vector *c) {
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
