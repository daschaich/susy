// -----------------------------------------------------------------
// Add real scalar to diagonal components of fundamental matrix
// a <-- a + cI
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_add_diag_f(matrix_f *a, Real c) {
  register int i;

  for (i = 0; i < NCOL; i++)
    a->e[i][i].real += c;
}
// -----------------------------------------------------------------
