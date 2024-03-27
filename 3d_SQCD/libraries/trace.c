// -----------------------------------------------------------------
// Return complex trace of the given matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

complex trace(matrix *a) {
  register complex tc;
  CADD(a->e[0][0], a->e[1][1], tc);
#if (NCOL > 2)
  CSUM(tc, a->e[2][2]);
#if (NCOL > 3)
  CSUM(tc, a->e[3][3]);
#if (NCOL > 4)
  CSUM(tc, a->e[4][4]);
#if (NCOL > 5)
  CSUM(tc, a->e[5][5]);
#if (NCOL > 6)
  CSUM(tc, a->e[6][6]);
#if (NCOL > 7)
  register int i;
  for (i = 7; i < NCOL; i++)
    CSUM(tc, a->e[i][i]);
#endif
#endif
#endif
#endif
#endif
#endif
  return tc;
}
// -----------------------------------------------------------------
