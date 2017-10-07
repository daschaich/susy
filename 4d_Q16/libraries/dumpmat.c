// -----------------------------------------------------------------
// Print the given matrix
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/susy.h"

void dumpmat(matrix *m) {
  int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      printf("  (%.4g, %.4g)", m->e[i][j].real, m->e[i][j].imag);
    printf("\n");
  }
}
// -----------------------------------------------------------------
