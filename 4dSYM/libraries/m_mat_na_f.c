// -----------------------------------------------------------------
// Fundamental matrix multiplication with adjoint of second matrix
// c <-- c + a * bdag
// c <-- c - a * bdag
// c <-- a * bdag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void mult_na_sum_f(matrix_f *a, matrix_f *b, matrix_f *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[j][k].real
                         + a->e[i][k].imag * b->e[j][k].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[j][k].real
                         - a->e[i][k].real * b->e[j][k].imag;
      }
    }
  }
}

void mult_na_dif_f(matrix_f *a, matrix_f *b, matrix_f *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real -= a->e[i][k].real * b->e[j][k].real
                         + a->e[i][k].imag * b->e[j][k].imag;
        c->e[i][j].imag -= a->e[i][k].imag * b->e[j][k].real
                         - a->e[i][k].real * b->e[j][k].imag;
      }
    }
  }
}

#ifndef FAST
void mult_na_f(matrix_f *a, matrix_f *b, matrix_f *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = a->e[i][0].real * b->e[j][0].real
                      + a->e[i][0].imag * b->e[j][0].imag;
      c->e[i][j].imag = a->e[i][0].imag * b->e[j][0].real
                      - a->e[i][0].real * b->e[j][0].imag;
      for (k = 1; k < NCOL; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[j][k].real
                         + a->e[i][k].imag * b->e[j][k].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[j][k].real
                         - a->e[i][k].real * b->e[j][k].imag;
      }
    }
  }
}
#else   // FAST version for NCOL=3 only
void mult_na(matrix_f *a, matrix_f *b, matrix_f *c) {
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
