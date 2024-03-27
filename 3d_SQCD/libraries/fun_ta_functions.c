//Binary operations where the first argument is a transpose funmatrix
//and the second an adjoint matrix. Potential 3rd argument for the result
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"


// ---------------------------------------
// Matrix multiplication
// c <-- a*b
void fun_ta_mult_nn(funmatrix *a, matrix *b, funmatrix *c) {
  register int i, j, k;
  for(i = 0; i<NCOLF ; i++)
    for(j = 0; j<NCOL ; j++)
      for(k = 0; k<NCOL ; k++) {
        c->e[j][i].real = a->e[k][i].real * b->e[k][j].real -
                          a->e[k][i].imag * b->e[k][j].imag;
        c->e[j][i].imag = a->e[k][i].real * b->e[k][j].imag +
                          a->e[k][i].imag * b->e[k][j].real;
      }
}


