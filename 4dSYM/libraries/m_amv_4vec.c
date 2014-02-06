// -----------------------------------------------------------------
// Four irrep matrix--vector operations with adjoint matrices
// ci <-- a[i]dag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void mult_adj_su3_mat_4vec(su3_matrix *a, su3_vector *b,
                           su3_vector *c0, su3_vector *c1,
                           su3_vector *c2, su3_vector *c3) {

  mult_adj_su3_mat_vec(a + 0, b, c0);
  mult_adj_su3_mat_vec(a + 1, b, c1);
  mult_adj_su3_mat_vec(a + 2, b, c2);
  mult_adj_su3_mat_vec(a + 3, b, c3);
}
#else   // FAST version for DIMF=3 only
void mult_adj_su3_mat_4vec(su3_matrix *a, su3_vector *b,
                           su3_vector *c0, su3_vector *c1,
                           su3_vector *c2, su3_vector *c3) {

  register int n;
  register Real c0r, c0i, c1r, c1i, c2r, c2i;
  register Real br, bi, a0, a1, a2;
  register su3_matrix *aa;    // To increment
  register su3_vector *c;
  su3_vector *cc[4];

  cc[0] = c0;
  cc[1] = c1;
  cc[2] = c2;
  cc[3] = c3;

  a = mat;
  c = c0;
  for (n = 0; n < 4; n++, a++, c = cc[n]) {
    br = b->c[0].real;
    bi = b->c[0].imag;
    a0 = aa->e[0][0].real;
    a1 = aa->e[0][1].real;
    a2 = aa->e[0][2].real;

    c0r = a0 * br;
    c1r = a1 * br;
    c2r = a2 * br;
    c0i = a0 * bi;
    c1i = a1 * bi;
    c2i = a2 * bi;

    a0 = aa->e[0][0].imag;
    a1 = aa->e[0][1].imag;
    a2 = aa->e[0][2].imag;

    c0r += a0 * bi;
    c1r += a1 * bi;
    c2r += a2 * bi;
    c0i -= a0 * br;
    c1i -= a1 * br;
    c2i -= a2 * br;

    br = b->c[1].real;
    bi = b->c[1].imag;
    a0 = aa->e[1][0].real;
    a1 = aa->e[1][1].real;
    a2 = aa->e[1][2].real;

    c0r += a0 * br;
    c1r += a1 * br;
    c2r += a2 * br;
    c0i += a0 * bi;
    c1i += a1 * bi;
    c2i += a2 * bi;

    a0 = aa->e[1][0].imag;
    a1 = aa->e[1][1].imag;
    a2 = aa->e[1][2].imag;

    c0r += a0 * bi;
    c1r += a1 * bi;
    c2r += a2 * bi;
    c0i -= a0 * br;
    c1i -= a1 * br;
    c2i -= a2 * br;

    br = b->c[2].real;
    bi = b->c[2].imag;
    a0 = aa->e[2][0].real;
    a1 = aa->e[2][1].real;
    a2 = aa->e[2][2].real;

    c0r += a0 * br;
    c1r += a1 * br;
    c2r += a2 * br;
    c0i += a0 * bi;
    c1i += a1 * bi;
    c2i += a2 * bi;

    a0 = aa->e[2][0].imag;
    a1 = aa->e[2][1].imag;
    a2 = aa->e[2][2].imag;

    c0r += a0 * bi;
    c1r += a1 * bi;
    c2r += a2 * bi;
    c0i -= a0 * br;
    c1i -= a1 * br;
    c2i -= a2 * br;

    c->c[0].real = c0r;
    c->c[0].imag = c0i;
    c->c[1].real = c1r;
    c->c[1].imag = c1i;
    c->c[2].real = c2r;
    c->c[2].imag = c2i;
  }
}
#endif
// -----------------------------------------------------------------
