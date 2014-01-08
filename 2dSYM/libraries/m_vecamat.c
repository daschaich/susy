// -----------------------------------------------------------------
// Irrep vector--matrix operation with adjoint matrix
// c <-- b * adag
// No FAST or REALREP implementation yet
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_su3_vec_adj_mat(su3_vector *b, su3_matrix *a, su3_vector *c) {
  register int i, j;
  register complex x, y, z;
  for (i = 0; i < DIMF; i++) {
    CONJG(a->e[i][0], z);
    CMUL(z, b->c[0], x)
    for (j = 1; j < DIMF; j++) {
      CONJG(a->e[i][j], z);
      CMUL(z, b->c[j], y)
      CSUM(x, y);
    }
    c->c[i] = x;
  }
}
// -----------------------------------------------------------------
