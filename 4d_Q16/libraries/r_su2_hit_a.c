// -----------------------------------------------------------------
// Multiply a matrix by the adjoint of an SU2 matrix from the right
// The p and q columns of link match the 0 and 1 columns of udag, respectively
// link <-- link * udag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void right_su2_hit_a(su2_matrix *u, int p, int q, matrix *link) {
  register int m;
  for (m = 0; m < NCOL; m++)
    mult_su2_mat_vec_elem_a(u, &(link->e[m][p]), &(link->e[m][q]));
}
// -----------------------------------------------------------------
