// -----------------------------------------------------------------
// Sum result of scalar multiplication on irrep vector
// a <-- a + s*b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_sum_vector(vector *a, vector *b, Real s) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    a->c[i].real += s * b->c[i].real;
    a->c[i].imag += s * b->c[i].imag;
  }
}
// -----------------------------------------------------------------
