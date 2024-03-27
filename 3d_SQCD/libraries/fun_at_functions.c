//c -= a*b dag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"


void fun_at_mult_an_sub(matrix *a, funmatrix *b, funmatrix *c) {
  register int i, j, k;
  for(i = 0; i<NCOL ; i++)
    for(j = 0; j<NCOLF ; j++)
      for(k = 0; k<NCOL ; k++) {
        c->e[i][j].real -= a->e[i][k].real * b->e[k][j].real +
                           a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[i][k].real * b->e[k][j].imag -
                           a->e[i][k].imag * b->e[k][j].real;
      }
}
