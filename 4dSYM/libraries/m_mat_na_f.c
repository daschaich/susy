// -----------------------------------------------------------------
// Fundamental matrix multiplication with adjoint of second matrix
// c <-- a * bdag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void mult_su3_na_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c) {
  register int i, j, k;
  register complex y;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      CMUL_J(a->e[i][0], b->e[j][0], c->e[i][j]);
      for (k = 1; k < NCOL; k++) {
        CMUL_J(a->e[i][k], b->e[j][k], y);
        CSUM(c->e[i][j], y);
      }
    }
  }
}
#else   // FAST version for NCOL=3 only
void mult_su3_na(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c) {
  int i,j;
  register Real t, ar, ai, br, bi, cr, ci;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      ar = a->e[i][0].real;
      ai = a->e[i][0].imag;
      br = b->e[j][0].real;
      bi = b->e[j][0].imag;
      cr = ar * br;
      t = ai * bi;
      cr += t;
      ci = ai * br;
      t = ar * bi;
      ci -= t;

      ar = a->e[i][1].real;
      ai = a->e[i][1].imag;
      br = b->e[j][1].real;
      bi = b->e[j][1].imag;
      t = ar * br;
      cr += t;
      t = ai * bi;
      cr += t;
      t = ar * bi;
      ci -= t;
      t = ai * br;
      ci += t;

      ar = a->e[i][2].real;
      ai = a->e[i][2].imag;
      br = b->e[j][2].real;
      bi = b->e[j][2].imag;
      t = ar * br;
      cr += t;
      t = ai * bi;
      cr += t;
      t = ar * bi;
      ci -= t;
      t = ai * br;
      ci += t;

      c->e[i][j].real = cr;
      c->e[i][j].imag = ci;
    }
  }
}
#endif
// -----------------------------------------------------------------
