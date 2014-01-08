// -----------------------------------------------------------------
// Complex scalar multiplication on irrep vector
// c <-- s*a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_su3vec(su3_vector *a, complex *s, su3_vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++)
    CMUL(a->c[i], *s, c->c[i]);
}
// -----------------------------------------------------------------
