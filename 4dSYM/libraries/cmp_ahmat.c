// -----------------------------------------------------------------
// Make a compressed anti-hermitian matrix from a fundamental matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void compress_anti_hermitian(su3_matrix_f *src, anti_hermitmat *dest) {
  dest->m00im = src->e[0][0].imag;
  dest->m11im = src->e[1][1].imag;
  dest->m01.real = src->e[0][1].real;
  dest->m01.imag = src->e[0][1].imag;
#if (NCOL > 2)
  dest->m22im = src->e[2][2].imag;
  dest->m02.real = src->e[0][2].real;
  dest->m12.real = src->e[1][2].real;
  dest->m02.imag = src->e[0][2].imag;
  dest->m12.imag = src->e[1][2].imag;
#endif
#if (NCOL > 3)
  dest->m33im = src->e[3][3].imag;
  dest->m03.real = src->e[0][3].real;
  dest->m13.real = src->e[1][3].real;
  dest->m23.real = src->e[2][3].real;
  dest->m03.imag = src->e[0][3].imag;
  dest->m13.imag = src->e[1][3].imag;
  dest->m23.imag = src->e[2][3].imag;
#endif
#if (NCOL > 4)
  node0_printf("compress_anti_hermitian only works for NCOL<=4!");
#endif
}
// -----------------------------------------------------------------
