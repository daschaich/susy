// -----------------------------------------------------------------
// Add two irrep vectors
// c <-- a + b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void add_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c) {
  CADD(a->c[0], b->c[0], c->c[0]);
  CADD(a->c[1], b->c[1], c->c[1]);
  CADD(a->c[2], b->c[2], c->c[2]);
  CADD(a->c[3], b->c[3], c->c[3]);
#if (DIMF > 4)
  CADD(a->c[4], b->c[4], c->c[4]);
  CADD(a->c[5], b->c[5], c->c[5]);
  CADD(a->c[6], b->c[6], c->c[6]);
  CADD(a->c[7], b->c[7], c->c[7]);
  CADD(a->c[8], b->c[8], c->c[8]);
#if (DIMF > 9)
  CADD(a->c[9], b->c[9], c->c[9]);
  CADD(a->c[10], b->c[10], c->c[10]);
  CADD(a->c[11], b->c[11], c->c[11]);
  CADD(a->c[12], b->c[12], c->c[12]);
  CADD(a->c[13], b->c[13], c->c[13]);
  CADD(a->c[14], b->c[14], c->c[14]);
  CADD(a->c[15], b->c[15], c->c[15]);
#if (DIMF > 16)
  CADD(a->c[16], b->c[16], c->c[16]);
  CADD(a->c[17], b->c[17], c->c[17]);
  CADD(a->c[18], b->c[18], c->c[18]);
  CADD(a->c[19], b->c[19], c->c[19]);
  CADD(a->c[20], b->c[20], c->c[20]);
  CADD(a->c[21], b->c[21], c->c[21]);
  CADD(a->c[22], b->c[22], c->c[22]);
  CADD(a->c[23], b->c[23], c->c[23]);
  CADD(a->c[24], b->c[24], c->c[24]);
#if (DIMF > 25)
  CADD(a->c[25], b->c[25], c->c[25]);
  CADD(a->c[26], b->c[26], c->c[26]);
  CADD(a->c[27], b->c[27], c->c[27]);
  CADD(a->c[28], b->c[28], c->c[28]);
  CADD(a->c[29], b->c[29], c->c[29]);
  CADD(a->c[30], b->c[30], c->c[30]);
  CADD(a->c[31], b->c[31], c->c[31]);
  CADD(a->c[32], b->c[32], c->c[32]);
  CADD(a->c[33], b->c[33], c->c[33]);
  CADD(a->c[34], b->c[34], c->c[34]);
  CADD(a->c[35], b->c[35], c->c[35]);
#if (DIMF > 36)
  CADD(a->c[36], b->c[36], c->c[36]);
  CADD(a->c[37], b->c[37], c->c[37]);
  CADD(a->c[38], b->c[38], c->c[38]);
  CADD(a->c[39], b->c[39], c->c[39]);
  CADD(a->c[40], b->c[40], c->c[40]);
  CADD(a->c[41], b->c[41], c->c[41]);
  CADD(a->c[42], b->c[42], c->c[42]);
  CADD(a->c[43], b->c[43], c->c[43]);
  CADD(a->c[44], b->c[44], c->c[44]);
  CADD(a->c[45], b->c[45], c->c[45]);
  CADD(a->c[46], b->c[46], c->c[46]);
  CADD(a->c[47], b->c[47], c->c[47]);
  CADD(a->c[48], b->c[48], c->c[48]);
#if (DIMF > 49)
  register int i;
  for (i = 49; i < DIMF; i++)
    CADD(a->c[i], b->c[i], c->c[i]);
#endif
#endif
#endif
#endif
#endif
#endif
}
// -----------------------------------------------------------------
