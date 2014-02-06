// -----------------------------------------------------------------
// Complex scalar multiplication on irrep matrix
// c <-- s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_su3mat(su3_matrix *b, complex *s, su3_matrix *c) {
#ifndef REALREP // Never used for REALREP
  register int i, j;
  register double sr, si, br, bi, cr, ci;

  sr = (*s).real;
  si = (*s).imag;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
    br = b->e[i][j].real;
    bi = b->e[i][j].imag;

    cr = sr * br - si * bi;
    ci = sr * bi + si * br;

    c->e[i][j].real = cr;
    c->e[i][j].imag = ci;
    }
  }
#endif
}
// -----------------------------------------------------------------
