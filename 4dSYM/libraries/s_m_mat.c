// -----------------------------------------------------------------
// Scalar multiplication on irrep matrix
// b <-- s * a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_su3_matrix(su3_matrix *a, Real s, su3_matrix *b) {
#ifndef FAST
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
#ifndef REALREP
      b->e[i][j].real = s * a->e[i][j].real;
      b->e[i][j].imag = s * a->e[i][j].imag;
#else
      b->e[i][j] = s * a->e[i][j];
#endif
    }
  }
#else   // FAST version for DIMF=3 only
  register Real ss = s;
  b->e[0][0].real = ss * a->e[0][0].real;
  b->e[0][0].imag = ss * a->e[0][0].imag;
  b->e[0][1].real = ss * a->e[0][1].real;
  b->e[0][1].imag = ss * a->e[0][1].imag;
  b->e[0][2].real = ss * a->e[0][2].real;
  b->e[0][2].imag = ss * a->e[0][2].imag;

  b->e[1][0].real = ss * a->e[1][0].real;
  b->e[1][0].imag = ss * a->e[1][0].imag;
  b->e[1][1].real = ss * a->e[1][1].real;
  b->e[1][1].imag = ss * a->e[1][1].imag;
  b->e[1][2].real = ss * a->e[1][2].real;
  b->e[1][2].imag = ss * a->e[1][2].imag;

  b->e[2][0].real = ss * a->e[2][0].real;
  b->e[2][0].imag = ss * a->e[2][0].imag;
  b->e[2][1].real = ss * a->e[2][1].real;
  b->e[2][1].imag = ss * a->e[2][1].imag;
  b->e[2][2].real = ss * a->e[2][2].real;
  b->e[2][2].imag = ss * a->e[2][2].imag;
#endif
}
// -----------------------------------------------------------------
