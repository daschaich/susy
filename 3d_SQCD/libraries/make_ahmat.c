// -----------------------------------------------------------------
// Compute and compress the traceless anti-hermitian part of a matrix
// dest = 0.5 * (src - src^dag) - Tr[0.5 * (src - src^dag)]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void make_anti_hermitian(matrix *src, anti_hermitmat *dest) {
  int i, j, index = 0;
  Real tr;

  tr = src->e[0][0].imag;
  for (i = 1; i < NCOL; i++)
    tr += src->e[i][i].imag;
  tr /= (Real)NCOL;

  for (i = 0; i < NCOL; i++)
    dest->im_diag[i] = src->e[i][i].imag - tr;

  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      dest->m[index].real = 0.5 * (src->e[i][j].real - src->e[j][i].real);
      dest->m[index].imag = 0.5 * (src->e[i][j].imag + src->e[j][i].imag);
      index++;
    }
  }
}
// -----------------------------------------------------------------
