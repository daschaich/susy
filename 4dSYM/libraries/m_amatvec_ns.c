// -----------------------------------------------------------------
// Subtract result of irrep matrix--vector operation with adjoint matrix
// c <-- c - adag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

#ifndef FAST
void mult_adj_mat_vec_dif(matrix *a, vector *b, vector *c) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      c->c[i].real -= a->e[j][i].real * b->c[j].real
                    + a->e[j][i].imag * b->c[j].imag;
      c->c[i].imag -= a->e[j][i].real * b->c[j].imag
                    - a->e[j][i].imag * b->c[j].real;
    }
  }
}
#else   // FAST version for DIMF=3 only
void mult_adj_mat_vec_dif(matrix *a, vector *b, vector *c) {
  register Real c0r, c0i, c1r, c1i, c2r, c2i;
  register Real br, bi, a0, a1, a2;

  br = b->c[0].real;
  bi = b->c[0].imag;
  a0 = a->e[0][0].real;
  a1 = a->e[0][1].real;
  a2 = a->e[0][2].real;

  c0r = a0 * br;
  c1r = a1 * br;
  c2r = a2 * br;
  c0i = a0 * bi;
  c1i = a1 * bi;
  c2i = a2 * bi;

  a0 = a->e[0][0].imag;
  a1 = a->e[0][1].imag;
  a2 = a->e[0][2].imag;

  c0r += a0 * bi;
  c1r += a1 * bi;
  c2r += a2 * bi;
  c0i -= a0 * br;
  c1i -= a1 * br;
  c2i -= a2 * br;

  br = b->c[1].real;
  bi = b->c[1].imag;
  a0 = a->e[1][0].real;
  a1 = a->e[1][1].real;
  a2 = a->e[1][2].real;

  c0r += a0 * br;
  c1r += a1 * br;
  c2r += a2 * br;
  c0i += a0 * bi;
  c1i += a1 * bi;
  c2i += a2 * bi;

  a0 = a->e[1][0].imag;
  a1 = a->e[1][1].imag;
  a2 = a->e[1][2].imag;

  c0r += a0 * bi;
  c1r += a1 * bi;
  c2r += a2 * bi;
  c0i -= a0 * br;
  c1i -= a1 * br;
  c2i -= a2 * br;

  br = b->c[2].real;
  bi = b->c[2].imag;
  a0 = a->e[2][0].real;
  a1 = a->e[2][1].real;
  a2 = a->e[2][2].real;

  c0r += a0 * br;
  c1r += a1 * br;
  c2r += a2 * br;
  c0i += a0 * bi;
  c1i += a1 * bi;
  c2i += a2 * bi;

  a0 = a->e[2][0].imag;
  a1 = a->e[2][1].imag;
  a2 = a->e[2][2].imag;

  c0r += a0 * bi;
  c1r += a1 * bi;
  c2r += a2 * bi;
  c0i -= a0 * br;
  c1i -= a1 * br;
  c2i -= a2 * br;

  c->c[0].real -= c0r;
  c->c[0].imag -= c0i;
  c->c[1].real -= c1r;
  c->c[1].imag -= c1i;
  c->c[2].real -= c2r;
  c->c[2].imag -= c2i;
}
#endif
// -----------------------------------------------------------------
