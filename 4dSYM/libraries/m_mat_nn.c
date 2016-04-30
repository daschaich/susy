// -----------------------------------------------------------------
// Irrep matrix multiplication with no adjoints
// c <-- a * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

#ifndef FAST
void mult_nn(matrix *a, matrix *b, matrix *c) {
  register int i, j, k;
  register complex y;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      CMUL(a->e[i][0], b->e[0][j], c->e[i][j]);   // Initialize
      for (k = 1; k < DIMF; k++) {
        CMUL(a->e[i][k], b->e[k][j], y);
        CSUM(c->e[i][j], y);
      }
    }
  }
}
#else   // FAST version for NCOL=3 only
void mult_nn(matrix *a, matrix *b, matrix *c) {
  int i, j;
  register Real t, ar, ai, br, bi, cr, ci;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      ar = a->e[i][0].real;
      ai = a->e[i][0].imag;
      br = b->e[0][j].real;
      bi = b->e[0][j].imag;
      cr = ar * br;
      t = ai * bi;
      cr -= t;
      ci = ar * bi;
      t = ai * br;
      ci += t;

      ar = a->e[i][1].real;
      ai = a->e[i][1].imag;
      br = b->e[1][j].real;
      bi = b->e[1][j].imag;
      t = ar * br;
      cr += t;
      t = ai * bi;
      cr -= t;
      t = ar * bi;
      ci += t;
      t = ai * br;
      ci += t;

      ar = a->e[i][2].real;
      ai = a->e[i][2].imag;
      br = b->e[2][j].real;
      bi = b->e[2][j].imag;
      t = ar * br;
      cr += t;
      t = ai * bi;
      cr -= t;
      t = ar * bi;
      ci += t;
      t = ai * br;
      ci += t;

      c->e[i][j].real = cr;
      c->e[i][j].imag = ci;
    }
  }
}
#endif
// -----------------------------------------------------------------
