// -----------------------------------------------------------------
// Print the given SU2 matrix
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dumpsu2(su2_matrix *u) {
  int i, j;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++)
      printf("(%.4g, %.4g) ", u->e[i][j].real, u->e[i][j].imag);
    printf("\n");
  }
  printf("\n");
}
// -----------------------------------------------------------------
