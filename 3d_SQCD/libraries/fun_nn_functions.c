// -----------------------------------------------------------------
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"





// Return real trace of matrix products adag * b
Real fun_nn_realtrace(funmatrix *a, funmatrix *b) {
  register int i, j;
  register Real sum = 0.0;

  for (i = 0; i < NCOL; i++) {
    for(j = 0; j < NCOLF; j++) {
      sum += a->e[i][j].real * b->e[i][j].real
           + a->e[i][j].imag * b->e[i][j].imag;
    }
  }
  return sum;
}
// -----------------------------------------------------------------
