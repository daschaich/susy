// -----------------------------------------------------------------
// Return real trace of irrep matrix product adag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

Real realtrace(matrix *a, matrix *b) {
  register int i, j;
  register Real sum = 0.0;

  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      sum += a->e[i][j].real * b->e[i][j].real
           + a->e[i][j].imag * b->e[i][j].imag;
    }
  }
  return sum;
}
// -----------------------------------------------------------------
