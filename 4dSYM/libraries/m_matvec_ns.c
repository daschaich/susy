// -----------------------------------------------------------------
// Subtract result of irrep matrix--vector operation with no adjoints
// c <-- c - a * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void mult_su3_mat_vec_nsum(su3_matrix *a, su3_vector *b, su3_vector *c) {
  register int i, j;
  register complex x, y;
  for (i = 0; i < DIMF; i++) {
    x.real = 0.0;
    x.imag = 0.0;
    for (j = 0; j < DIMF; j++) {
      CMUL(a->e[i][j], b->c[j], y);
      CSUM(x, y);
    }
    c->c[i].real -= x.real;
    c->c[i].imag -= x.imag;
  }
}
#else   // FAST version for DIMF=3 only
void mult_su3_mat_vec_nsum(su3_matrix *a, su3_vector *b, su3_vector *c) {
  int i;
  register Real t, ar, ai, br, bi, cr, ci;
  for(i = 0; i < 3; i++) {
    ar = a->e[i][0].real;
    ai = a->e[i][0].imag;
    br = b->c[0].real;
    bi = b->c[0].imag;
    cr = ar * br;
    t = ai * bi;
    cr -= t;
    ci = ar * bi;
    t = ai * br;
    ci += t;

    ar = a->e[i][1].real;
    ai = a->e[i][1].imag;
    br = b->c[1].real;
    bi = b->c[1].imag;
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
    br = b->c[2].real;
    bi = b->c[2].imag;
    t = ar*br;
    cr += t;
    t = ai*bi;
    cr -= t;
    t = ar * bi;
    ci += t;
    t = ai * br;
    ci += t;

    c->c[i].real -= cr;
    c->c[i].imag -= ci;
  }
}
#endif
// -----------------------------------------------------------------
