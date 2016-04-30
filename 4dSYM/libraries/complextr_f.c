// -----------------------------------------------------------------
// Return complex trace of fundamental matrix product adag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

complex complextrace_f(matrix_f *a, matrix_f *b) {
  register int i, j;
  register Real sumr = 0.0, sumi = 0.0;
  complex sum;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      sumr += a->e[i][j].real * b->e[i][j].real
            + a->e[i][j].imag * b->e[i][j].imag;
      sumi += a->e[i][j].real * b->e[i][j].imag
            - a->e[i][j].imag * b->e[i][j].real;
    }
  }
  sum.real = sumr;
  sum.imag = sumi;
  return sum;
}
// -----------------------------------------------------------------
