// -----------------------------------------------------------------
// Adjoint of a funmatrix scaled by a scalar
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"


// b <-- R * adag
void fun_scaled_adjoint(funmatrix *a, Real R, funmatrix *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      b->e[i][j].real = R *a->e[i][j].real;
      b->e[i][j].imag = - R * a->e[i][j].imag;
    }
  }
}

// b += R* adag
void fun_scaled_adjoint_sum(funmatrix *a, Real R, funmatrix *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      b->e[i][j].real += R *a->e[i][j].real;
      b->e[i][j].imag += - R * a->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
