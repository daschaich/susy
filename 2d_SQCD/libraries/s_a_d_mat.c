// -----------------------------------------------------------------
// Add or subtract real scalar to diagonal components of matrix
// a <-- a + cI
// a <-- a - cI
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_add_diag(matrix *a, Real c) {
  register int i;

  for (i = 0; i < NCOL; i++)
    a->e[i][i].real += c;
}

void scalar_sub_diag(matrix *a, Real c) {
  register int i;

  for (i = 0; i < NCOL; i++)
    a->e[i][i].real -= c;
}
// -----------------------------------------------------------------
