// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force
#include "susy_includes.h"
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
#ifndef PUREGAUGE
void fermion_op(matrix *src[NFERMION], matrix *dest[NFERMION], int sign) {
  register int i, j;
  register site *s;
  matrix tmat;
  
  if (sign == -1) {
    FORALLSITES(i, s) {
      for(j=0;j<NFERMION;j++)
      {
        adjoint(&(src[j][i]), &tmat);
        mat_copy(&tmat, &(src[j][i]));
      }
    }
  }
  
  // TODO: Fermion operator goes here

  
  if (sign == -1) {    // Both negate and conjugate
    FORALLSITES(i, s) {
      for(j=0;j<NFERMION;j++)
      {
        adjoint(&(dest[j][i]), &tmat);
        scalar_mult_matrix(&tmat, -1.0 , &(dest[j][i]));
      }
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Squared matrix--vector operation
//   dest = (D^2).src
// Use temp_ferm for temporary storage
#ifndef PUREGAUGE
void DSq(matrix *src[NFERMION], matrix *dest[NFERMION]) {
  fermion_op(src, temp_ferm, PLUS);
  fermion_op(temp_ferm, dest, MINUS);
}
#endif
// -----------------------------------------------------------------
