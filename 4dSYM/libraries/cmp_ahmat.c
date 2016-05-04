// -----------------------------------------------------------------
// Make a compressed anti-hermitian matrix from a matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void compress_anti_hermitian(matrix *src, anti_hermitmat *dest) {
  int i, j, index = 0;
  for (i = 0; i < NCOL; i++)
    dest->im_diag[i] = src->e[i][i].imag;

  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      dest->m[index].real = src->e[i][j].real;
      dest->m[index].imag = src->e[i][j].imag;
      index++;
    }
  }
}
// -----------------------------------------------------------------
