// -----------------------------------------------------------------
// Print the anti-hermitian part of the given fundamental matrix
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/susy.h"

void dump_ahmat_f(matrix_f *m) {
  int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      printf("(%.4g, %.4g) ",
             (m->e[i][j].real - m->e[j][i].real) * 0.5,
             (m->e[i][j].imag + m->e[j][i].imag) * 0.5);
    }
    printf("\n");
  }
  printf("\n");
}
// -----------------------------------------------------------------
