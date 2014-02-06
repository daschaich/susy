// -----------------------------------------------------------------
// Array of four irrep matrix--vector operations with adjoint matrices
// Specialized case of m_amv_4vec.c
// c[i] <-- a[i]dag * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
void mult_adj_su3_mat_vec_4dir(su3_matrix *a, su3_vector *b,
                               su3_vector *c) {

  mult_adj_su3_mat_vec(a + 0, b, c + 0);
  mult_adj_su3_mat_vec(a + 1, b, c + 1);
  mult_adj_su3_mat_vec(a + 2, b, c + 2);
  mult_adj_su3_mat_vec(a + 3, b, c + 3);
}
#else   // FAST version for DIMF=3 only
void mult_adj_su3_mat_vec_4dir(su3_matrix *a, su3_vector *b,
                               su3_vector *c) {

  register int n;
  register Real c0r, c0i, c1r, c1i, c2r, c2i;
  register Real br, bi, a0, a1, a2;
  register su3_matrix *aa = a;    // To increment
  register su3_vector *cc = c;

  for (n = 0; n < 4; n++, aa++, cc++) {
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

    cc->c[0].real = c0r;
    cc->c[0].imag = c0i;
    cc->c[1].real = c1r;
    cc->c[1].imag = c1i;
    cc->c[2].real = c2r;
    cc->c[2].imag = c2i;
  }
}
#endif
// -----------------------------------------------------------------
