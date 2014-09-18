// -----------------------------------------------------------------
// Return complex trace of the given irrep matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex trace_su3(su3_matrix *a) {
  register complex tc;
  int i;
  CADD(a->e[0][0], a->e[1][1], tc);
  for (i = 2; i < DIMF; i++)
    CADD(tc, a->e[i][i], tc);
  return tc;
}
// -----------------------------------------------------------------
