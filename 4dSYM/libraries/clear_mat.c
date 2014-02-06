// -----------------------------------------------------------------
// Clear an irrep matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clear_su3mat(su3_matrix *dest) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
#ifndef REALREP
      dest->e[i][j].real = 0.0;
      dest->e[i][j].imag = 0.0;
#else
      dest->e[i][j] = 0.0;
#endif
    }
  }
}
// -----------------------------------------------------------------
