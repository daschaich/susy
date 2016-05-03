// -----------------------------------------------------------------
// Add result of scalar multiplication on irrep vector
// c <-- c + s * b
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_sum_vector(vector *b, Real s, vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real += s * b->c[i].real;
    c->c[i].imag += s * b->c[i].imag;
  }
}

void scalar_mult_add_vector(vector *a, vector *b, Real s, vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real = a->c[i].real + s * b->c[i].real;
    c->c[i].imag = a->c[i].imag + s * b->c[i].imag;
  }
}
// -----------------------------------------------------------------
