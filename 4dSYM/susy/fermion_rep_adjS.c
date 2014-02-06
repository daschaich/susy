// -----------------------------------------------------------------
// Create link in adjoint rep from linkf in fundamental rep
#include "susy_includes.h"

#if (DIMF != NCOL*NCOL)
  #error "Wrong version of fermion_rep!"
#endif

#if (NUMGEN != NCOL*NCOL)
  #error "Wrong version of fermion_rep!"
#endif

#if (FREP != adjoint)
  #error "Wrong version of fermion_rep!"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void make_fermion_rep_matrix(su3_matrix_f *a, su3_matrix *b) {
  int i, j;
  su3_matrix_f tmat1, tmat2;
  for (i = 0; i < NUMGEN; i++) {
    for (j = 0; j < NUMGEN; j++) {
      mult_su3_nn_f(a, &Lambda[j], &tmat1);
      mult_su3_nn_f(&Lambda[i], &tmat1, &tmat2);
      b->e[i][j] = trace_su3_f(&tmat2);
    }
  }
}
// -----------------------------------------------------------------
