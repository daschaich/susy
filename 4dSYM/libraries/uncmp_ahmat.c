// -----------------------------------------------------------------
// Uncompresses an anti-hermitian matrix into a fundamental matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void uncompress_anti_hermitian(anti_hermitmat *src, matrix_f *dest) {
  int i, j, index = 0;
  Real tr;

  for (i = 0; i < NCOL; i++) {
    dest->e[i][i].imag = src->im_diag[i];
    dest->e[i][i].real = 0.0;
  }
  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      dest->e[i][j].imag = src->m[index].imag;
      dest->e[j][i].imag = src->m[index].imag;

      tr = src->m[index].real;
      dest->e[i][j].real = tr;
      dest->e[j][i].real = -tr;
      index++;
    }
  }
}
// -----------------------------------------------------------------
