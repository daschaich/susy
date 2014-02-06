// -----------------------------------------------------------------
// Take and compress the traceless anti-hermitian part
// of a fundamental matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Non-FAST version
#ifndef FAST
void make_anti_hermitian(su3_matrix_f *src, anti_hermitmat *dest) {
  Real temp;

#if (NCOL==2)
  temp = (src->e[0][0].imag + src->e[1][1].imag) * 0.5;
#elif (NCOL==3)
  temp = (src->e[0][0].imag + src->e[1][1].imag + src->e[2][2].imag)
         * 0.33333333333333333;
#elif (NCOL==4)
  temp = (src->e[0][0].imag + src->e[1][1].imag
          + src->e[2][2].imag + src->e[3][3].imag) * 0.25;
#endif
  dest->m00im = src->e[0][0].imag - temp;
  dest->m11im = src->e[1][1].imag - temp;
  dest->m01.real = (src->e[0][1].real - src->e[1][0].real) * 0.5;
  dest->m01.imag = (src->e[0][1].imag + src->e[1][0].imag) * 0.5;
#if (NCOL>2)
  dest->m22im = src->e[2][2].imag - temp;
  dest->m02.real = (src->e[0][2].real - src->e[2][0].real) * 0.5;
  dest->m12.real = (src->e[1][2].real - src->e[2][1].real) * 0.5;
  dest->m02.imag = (src->e[0][2].imag + src->e[2][0].imag) * 0.5;
  dest->m12.imag = (src->e[1][2].imag + src->e[2][1].imag) * 0.5;
#endif
#if (NCOL>3)
  dest->src3im = src->e[3][3].imag - temp;
  dest->m03.real = (src->e[0][3].real - src->e[3][0].real) * 0.5;
  dest->m13.real = (src->e[1][3].real - src->e[3][1].real) * 0.5;
  dest->m23.real = (src->e[2][3].real - src->e[3][2].real) * 0.5;
  dest->m03.imag = (src->e[0][3].imag + src->e[3][0].imag) * 0.5;
  dest->m13.imag = (src->e[1][3].imag + src->e[3][1].imag) * 0.5;
  dest->m23.imag = (src->e[2][3].imag + src->e[3][2].imag) * 0.5;
#endif
#if (NCOL>4)
  node0_printf("make_anti_hermitian only works for NCOL<=4!\n");
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// FAST version
#else
void make_anti_hermitian( su3_matrix *src, anti_hermitmat *dest ) {
  Real temp, temp2;
  temp = (src->e[0][0].imag + src->e[1][1].imag);
  temp2 = temp + src->e[2][2].imag;
  temp = temp2 * 0.333333333333333333;
  dest->m00im = src->e[0][0].imag - temp;
  dest->m11im = src->e[1][1].imag - temp;
  dest->m22im = src->e[2][2].imag - temp;

  temp = src->e[0][1].real - src->e[1][0].real;
  dest->m01.real = temp * 0.5;
  temp = src->e[0][2].real - src->e[2][0].real;
  dest->m02.real = temp * 0.5;
  temp = src->e[1][2].real - src->e[2][1].real;
  dest->m12.real = temp * 0.5;
  temp = src->e[0][1].imag + src->e[1][0].imag;
  dest->m01.imag = temp * 0.5;
  temp = src->e[0][2].imag + src->e[2][0].imag;
  dest->m02.imag = temp * 0.5;
  temp = src->e[1][2].imag + src->e[2][1].imag;
  dest->m12.imag = temp * 0.5;
}
#endif
// -----------------------------------------------------------------
