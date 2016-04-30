// -----------------------------------------------------------------
// Return real trace of fundamental matrix product adag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

Real realtrace_f(matrix_f *a, matrix_f *b) {
  register int i, j;
  register Real sum = 0.0;

  for (i = 0; i < NCOL; i++) {
    for(j = 0; j < NCOL; j++) {
      sum += a->e[i][j].real * b->e[i][j].real
           + a->e[i][j].imag * b->e[i][j].imag;
    }
  }
  return sum;
}
// -----------------------------------------------------------------
