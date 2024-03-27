// -----------------------------------------------------------------
// Clear the given matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void clear_mat(matrix *m) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      m->e[i][j].real = 0.0;
      m->e[i][j].imag = 0.0;
    }
  }
}
// -----------------------------------------------------------------
