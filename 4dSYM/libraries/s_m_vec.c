// -----------------------------------------------------------------
// Scalar multiplication on irrep vector
// c <-- s*a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_vector(vector *a, Real s, vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real = s * a->c[i].real;
    c->c[i].imag = s * a->c[i].imag;
  }
}
// -----------------------------------------------------------------
