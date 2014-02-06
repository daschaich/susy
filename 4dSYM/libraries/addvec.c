// -----------------------------------------------------------------
// Add two irrep vectors
// c <-- a + b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void add_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++)
    CADD(a->c[i], b->c[i], c->c[i]);
}
// -----------------------------------------------------------------
