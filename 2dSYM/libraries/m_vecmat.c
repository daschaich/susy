// -----------------------------------------------------------------
// Irrep vector--matrix operation with no adjoints
// c <-- b * a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_su3_vec_mat(su3_vector *b, su3_matrix *a, su3_vector *c) {
  register int i, j;
  register complex x, y;
  for (i = 0; i < DIMF; i++) {
    CMUL(a->e[0][i], b->c[0], x);
    for (j = 1; j < DIMF; j++) {
      CMUL(a->e[j][i], b->c[j], y)
      CSUM(x, y);
    }
    c->c[i] = x;
  }
}
// -----------------------------------------------------------------
