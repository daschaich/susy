// -----------------------------------------------------------------
// Complex scalar multiplication on irrep vector
// c <-- s*a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void c_scalar_mult_vec(vector *a, complex *s, vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++)
    CMUL(a->c[i], *s, c->c[i]);
}
// -----------------------------------------------------------------
