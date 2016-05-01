// -----------------------------------------------------------------
// Subtract result of scalar multiplication on irrep vector
// c <-- c - s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_dif_vector(vector *b, Real s, vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real -= s * b->c[i].real;
    c->c[i].imag -= s * b->c[i].imag;
  }
}
// -----------------------------------------------------------------
