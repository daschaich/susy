// -----------------------------------------------------------------
// Clear a fundamental matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clear_su3mat_f(su3_matrix_f *dest) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      dest->e[i][j].real = 0.0;
      dest->e[i][j].imag = 0.0;
    }
  }
}
// -----------------------------------------------------------------
