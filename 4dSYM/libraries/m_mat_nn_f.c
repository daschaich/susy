// -----------------------------------------------------------------
// Fundamental matrix multiplication with no adjoints
// c <-- c + a * b
// c <-- c - a * b
// c <-- a * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void mult_nn_sum_f(matrix_f *a, matrix_f *b, matrix_f *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[k][j].real
                         - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[k][j].real
                         + a->e[i][k].real * b->e[k][j].imag;
      }
    }
  }
}

void mult_nn_dif_f(matrix_f *a, matrix_f *b, matrix_f *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        c->e[i][j].real -= a->e[i][k].real * b->e[k][j].real
                         - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag -= a->e[i][k].imag * b->e[k][j].real
                         + a->e[i][k].real * b->e[k][j].imag;
      }
    }
  }
}

#ifndef FAST
void mult_nn_f(matrix_f *a, matrix_f *b, matrix_f *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = a->e[i][0].real * b->e[0][j].real
                      - a->e[i][0].imag * b->e[0][j].imag;
      c->e[i][j].imag = a->e[i][0].imag * b->e[0][j].real
                      + a->e[i][0].real * b->e[0][j].imag;
      for (k = 1; k < NCOL; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[k][j].real
                         - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[k][j].real
                         + a->e[i][k].real * b->e[k][j].imag;
      }
    }
  }
}
#else   // FAST version for NCOL=3 only
void mult_nn_f(matrix_f *a, matrix_f *b, matrix_f *c) {
  int i, j;
  register Real t, ar, ai, br, bi, cr, ci;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      ar = a->e[i][0].real;
      ai = a->e[i][0].imag;
      br = b->e[0][j].real;
      bi = b->e[0][j].imag;
      cr = ar*br;
      t = ai*bi;
      cr -= t;
      ci = ar*bi;
      t = ai*br;
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
