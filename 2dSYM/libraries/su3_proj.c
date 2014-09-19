// -----------------------------------------------------------------
// Outer product of irrep vectors
// c_ij <-- A_i Bdag_j
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void su3_projector(su3_vector *a, su3_vector *b, su3_matrix *c) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++)
      CMUL_J(a->c[i], b->c[j], c->e[i][j]);
  }
}
#else   // FAST version for DIMF=3 only
void su3_projector(su3_vector *a, su3_vector *b, su3_matrix *c) {
  register int i, j;
  register Real tmp, tmp2;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      tmp2 = a->c[i].real * b->c[j].real;
      tmp = a->c[i].imag * b->c[j].imag;
      c->e[i][j].real = tmp + tmp2;

      tmp2 = a->c[i].real * b->c[j].imag;
      tmp = a->c[i].imag * b->c[j].real;
      c->e[i][j].imag = tmp - tmp2;
    }
  }
}
#endif
// -----------------------------------------------------------------
