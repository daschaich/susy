// -----------------------------------------------------------------
// Add two irrep vectors
// c <-- a + b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void add_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c) {
  CADD(a->c[0], b->c[0], c->c[0]);
  CADD(a->c[1], b->c[1], c->c[1]);
#if (DIMF > 2)
  CADD(a->c[2], b->c[2], c->c[2]);
#if (DIMF > 3)
  CADD(a->c[3], b->c[3], c->c[3]);
#if (DIMF > 4)
  CADD(a->c[4], b->c[4], c->c[4]);
#if (DIMF > 5)
  CADD(a->c[5], b->c[5], c->c[5]);
#if (DIMF > 6)
  CADD(a->c[6], b->c[6], c->c[6]);
#if (DIMF > 7)
  register int i;
  for (i = 7; i < DIMF; i++)
    CADD(a->c[i], b->c[i], c->c[i]);
#endif
#endif
#endif
#endif
#endif
#endif
}
// -----------------------------------------------------------------
