#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"
//Matrix addition/subtraction

//diff
//c<-a-b

void fun_o_sub(funmatrix *a, funmatrix *b, funmatrix *c) {
  register int i, j;
  for(i = 0; i<NCOL; i++)
    for(j = 0; j<NCOLF; j++) {
      c->e[i][j].real = a->e[i][j].real - b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag - b->e[i][j].imag;
    }
}

//sum
//b=a+b
void fun_o_sum(funmatrix *a, funmatrix *b) {
  register int i, j;
  for(i = 0; i<NCOL; i++)
    for(j = 0; j<NCOLF; j++) {
      b->e[i][j].real += a->e[i][j].real;
      b->e[i][j].imag += a->e[i][j].imag;
    }
}




//b+=R*a
void fun_o_scalar_mult_sum(funmatrix *a, Real R, funmatrix *b) {
  register int i, j;
  for (i = 0; i < NCOL ; i++)
    for (j = 0; j < NCOLF ; j++) {
      b->e[i][j].real+=R*a->e[i][j].real;
      b->e[i][j].imag+=R*a->e[i][j].imag;
    }
}

void funmat_copy(funmatrix *src, funmatrix *dest) {
  *dest = *src;
}
