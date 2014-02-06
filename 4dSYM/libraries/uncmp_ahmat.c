// -----------------------------------------------------------------
// Uncompresses an anti-hermitian matrix into a fundamental matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void uncompress_anti_hermitian(anti_hermitmat *src, su3_matrix_f *dest) {
  Real tr;

  dest->e[0][0].imag = src->m00im;
  dest->e[0][0].real = 0.0;
  dest->e[1][1].imag = src->m11im;
  dest->e[1][1].real = 0.0;
  dest->e[0][1].imag = src->m01.imag;

  tr = src->m01.real;
  dest->e[0][1].real = tr;
  dest->e[1][0].real = -tr;
  dest->e[1][0].imag = src->m01.imag;
#if (NCOL > 2)
  dest->e[2][2].imag = src->m22im;
  dest->e[2][2].real = 0.0;
  dest->e[0][2].imag = src->m02.imag;

  tr = src->m02.real;
  dest->e[0][2].real = tr;
  dest->e[2][0].real = -tr;
  dest->e[2][0].imag = src->m02.imag;
  dest->e[1][2].imag = src->m12.imag;

  tr = src->m12.real;
  dest->e[1][2].real = tr;
  dest->e[2][1].real = -tr;
  dest->e[2][1].imag = src->m12.imag;
#endif
#if (NCOL > 3)
  dest->e[3][3].imag = src->m33im;
  dest->e[3][3].real = 0.0;
  dest->e[0][3].imag = src->m03.imag;

  tr = src->m03.real;
  dest->e[0][3].real = tr;
  dest->e[3][0].real = -tr;
  dest->e[3][0].imag = src->m03.imag;
  dest->e[1][3].imag = src->m13.imag;

  tr = src->m13.real;
  dest->e[1][3].real = tr;
  dest->e[3][1].real = -tr;
  dest->e[3][1].imag = src->m13.imag;
  dest->e[2][3].imag = src->m23.imag;

  tr = src->m23.real;
  dest->e[2][3].real = tr;
  dest->e[3][2].real = -tr;
  dest->e[3][2].imag = src->m23.imag;
#endif
#if (NCOL > 4)
  node0_printf("uncompress_anti_hermitian only works for NCOL<=4!");
#endif
}
// -----------------------------------------------------------------
