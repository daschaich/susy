//Binary operations between two transpose funmatrices
//Third parameter possible for the result
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"




//Matrix multiplication
//c<-a*b

void fun_tt_mult_an(funmatrix *a, funmatrix *b, matrix *c) {
  register i, j, k;
  for(i = 0; i<NCOL; i++)
    for(j = 0; j<NCOL; j++)
      for(k = 0; k<NCOLF; k++) {
        c->e[i][j].real = a->e[i][k].real * b->e[j][k].real +
                          a->e[i][k].imag * b->e[j][k].imag;
        c->e[i][j].imag = a->e[i][k].real * b->e[j][k].imag -
                          a->e[i][k].imag * b->e[j][k].real;
      }
}

//c+=a*b
void fun_tt_mult_an_sum(funmatrix *a, funmatrix *b, matrix *c) {
  register i, j, k;
  for(i = 0; i<NCOL; i++)
    for(j = 0; j<NCOL; j++)
      for(k = 0; k<NCOLF; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[j][k].real +
                          a->e[i][k].imag * b->e[j][k].imag;
        c->e[i][j].imag += a->e[i][k].real * b->e[j][k].imag -
                          a->e[i][k].imag * b->e[j][k].real;
      }
}
