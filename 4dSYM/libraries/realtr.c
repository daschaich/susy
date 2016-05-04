// -----------------------------------------------------------------
// Return real trace of matrix products a * b and adag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

Real realtrace_nn_f(matrix *a, matrix *b) {
  register int i, j;
  register Real sum = 0.0;

  for (i = 0; i < NCOL; i++) {
    for(j = 0; j < NCOL; j++) {
      sum += a->e[i][j].real * b->e[j][i].real
           - a->e[i][j].imag * b->e[j][i].imag;
    }
  }
  return sum;
}

Real realtrace(matrix *a, matrix *b) {
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
