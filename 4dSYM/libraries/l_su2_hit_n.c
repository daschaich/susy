// -----------------------------------------------------------------
// Multiply an irrep matrix by an SU2 matrix from the left
// The p and q rows of link match the 0 and 1 rows of u, respectively
// link <-- u * link
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void left_su2_hit_n(su2_matrix *u, int p, int q, su3_matrix *link) {
  register int m;

  for (m = 0; m < DIMF; m++)
    mult_su2_mat_vec_elem_n(u, &(link->e[p][m]), &(link->e[q][m]));
}
// -----------------------------------------------------------------
