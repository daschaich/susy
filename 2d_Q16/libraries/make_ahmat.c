// -----------------------------------------------------------------
// Take and compress the traceless anti-hermitian part of a matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void make_anti_hermitian(matrix *src, anti_hermitmat *dest) {
  int i, j, index = 0;
  Real temp;

  temp = src->e[0][0].imag;
  for (i = 1; i < NCOL; i++)
    temp += src->e[i][i].imag;
  temp /= (Real)NCOL;

  for (i = 0; i < NCOL; i++)
    dest->im_diag[i] = src->e[i][i].imag - temp;

  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      dest->m[index].real = (src->e[i][j].real - src->e[j][i].real) * 0.5;
      dest->m[index].imag = (src->e[i][j].imag + src->e[j][i].imag) * 0.5;
      index++;
    }
  }
}
// -----------------------------------------------------------------
