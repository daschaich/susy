// -----------------------------------------------------------------
// Sum results of four irrep matrix--vector operations with no adjoints
// c <-- a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void mult_su3_mat_vec_sum_4dir(su3_matrix *a, su3_vector *b0,
                               su3_vector *b1, su3_vector *b2,
                               su3_vector *b3, su3_vector *c) {

    mult_su3_mat_vec(a + 0, b0, c);
    mult_su3_mat_vec_sum(a + 1, b1, c);
    mult_su3_mat_vec_sum(a + 2, b2, c);
    mult_su3_mat_vec_sum(a + 3, b3, c);
}
#else   // FAST version for DIMF=3 only
void mult_su3_mat_vec_sum_4dir(su3_matrix *a, su3_vector *b0,
                               su3_vector *b1, su3_vector *b2,
                               su3_vector *b3, su3_vector *c) {

  register int n;
  register Real c0r, c0i, c1r, c1i, c2r, c2i;
  register Real br, bi, a0, a1, a2;
  register su3_matrix *mat;
  register su3_vector *b;

  c0r = 0.0;
  c0i = 0.0;
  c1r = 0.0;
  c1i = 0.0;
  c2r = 0.0;
  c2i = 0.0;
  mat = a;

  for (n = 0; n < 4; n++, mat++) {
    switch(n) {
      case(0): b = b0; break;
      case(1): b = b1; break;
      case(2): b = b2; break;
      case(3): b = b3; break;
      default: b = 0;
    }

    br = b->c[0].real;
    bi = b->c[0].imag;
    a0 = mat->e[0][0].real;
    a1 = mat->e[1][0].real;
    a2 = mat->e[2][0].real;

    c0r += a0 * br;
    c1r += a1 * br;
    c2r += a2 * br;
    c0i += a0 * bi;
    c1i += a1 * bi;
    c2i += a2 * bi;

    a0 = mat->e[0][0].imag;
    a1 = mat->e[1][0].imag;
    a2 = mat->e[2][0].imag;

    c0r -= a0 * bi;
    c1r -= a1 * bi;
    c2r -= a2 * bi;
    c0i += a0 * br;
    c1i += a1 * br;
    c2i += a2 * br;

    br = b->c[1].real;
    bi = b->c[1].imag;
    a0 = mat->e[0][1].real;
    a1 = mat->e[1][1].real;
    a2 = mat->e[2][1].real;

    c0r += a0 * br;
    c1r += a1 * br;
    c2r += a2 * br;
    c0i += a0 * bi;
    c1i += a1 * bi;
    c2i += a2 * bi;

    a0 = mat->e[0][1].imag;
    a1 = mat->e[1][1].imag;
    a2 = mat->e[2][1].imag;

    c0r -= a0 * bi;
    c1r -= a1 * bi;
    c2r -= a2 * bi;
    c0i += a0 * br;
    c1i += a1 * br;
    c2i += a2 * br;

    br = b->c[2].real;    bi = b->c[2].imag;
    a0 = mat->e[0][2].real;
    a1 = mat->e[1][2].real;
    a2 = mat->e[2][2].real;

    c0r += a0 * br;
    c1r += a1 * br;
    c2r += a2 * br;
    c0i += a0 * bi;
    c1i += a1 * bi;
    c2i += a2 * bi;

    a0 = mat->e[0][2].imag;
    a1 = mat->e[1][2].imag;
    a2 = mat->e[2][2].imag;

    c0r -= a0 * bi;
    c1r -= a1 * bi;
    c2r -= a2 * bi;
    c0i += a0 * br;
    c1i += a1 * br;
    c2i += a2 * br;

  }
  c->c[0].real = c0r;
  c->c[0].imag = c0i;
  c->c[1].real = c1r;
  c->c[1].imag = c1i;
  c->c[2].real = c2r;
  c->c[2].imag = c2i;
}
#endif
// -----------------------------------------------------------------
