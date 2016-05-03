// -----------------------------------------------------------------
// Irrep matrix--vector operation with adjoint matrix
// c <-- adag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

#ifndef FAST
void mult_adj_mat_vec(matrix *a, vector *b, vector *c) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    // Initialize
    c->c[i].real = a->e[0][i].real * b->c[0].real
                 + a->e[0][i].imag * b->c[0].imag;
    c->c[i].imag = a->e[0][i].real * b->c[0].imag
                 - a->e[0][i].imag * b->c[0].real;
    for (j = 1; j < DIMF; j++) {
      c->c[i].real += a->e[j][i].real * b->c[j].real
                    + a->e[j][i].imag * b->c[j].imag;
      c->c[i].imag += a->e[j][i].real * b->c[j].imag
                    - a->e[j][i].imag * b->c[j].real;
    }
  }
}
#else   // FAST version for DIMF=3 only
void mult_adj_mat_vec(matrix *a, vector *b, vector *c) {
  int i;
  register Real t, ar, ai, br, bi, cr, ci;
  for (i = 0; i < 3; i++) {
    ar = a->e[0][i].real;
    ai = a->e[0][i].imag;
    br = b->c[0].real;
    bi = b->c[0].imag;
    cr = ar * br;
    t = ai * bi;
    cr += t;
    ci = ar * bi;
    t = ai * br;
    ci -= t;

    ar = a->e[1][i].real;
    ai = a->e[1][i].imag;
    br = b->c[1].real;
    bi = b->c[1].imag;
    t = ar * br;
    cr += t;
    t = ai * bi;
    cr += t;
    t = ar * bi;
    ci += t;
    t = ai * br;
    ci -= t;

    ar = a->e[2][i].real;
    ai = a->e[2][i].imag;
    br = b->c[2].real;
    bi = b->c[2].imag;
    t = ar * br;
    cr += t;
    t = ai * bi;
    cr += t;
    t = ar * bi;
    ci += t;
    t = ai * br;
    ci -= t;

    c->c[i].real = cr;
    c->c[i].imag = ci;
  }
}
#endif
// -----------------------------------------------------------------
