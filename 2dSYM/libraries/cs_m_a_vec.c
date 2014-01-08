// -----------------------------------------------------------------
// Sum result of complex scalar multiplication on irrep vector
// a <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_add_su3vec(su3_vector *a, complex *s, su3_vector *b) {
  register int i;
  complex t;
  for (i = 0; i < DIMF; i++) {
    CMUL(b->c[i], *s, t);
    CADD(a->c[i], t, a->c[i]);
  }
}
// -----------------------------------------------------------------
