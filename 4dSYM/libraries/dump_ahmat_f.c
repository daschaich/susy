// -----------------------------------------------------------------
// Print the anti-hermitian part of the given fundamental matrix
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dump_ahmat_f(su3_matrix_f *m) {
  int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      printf("(%.4g, %.4g) ",
             (m->e[i][j].real - m->e[j][i].real) / 2.0,
             (m->e[i][j].imag + m->e[j][i].imag) / 2.0);
    }
    printf("\n");
  }
  printf("\n");
}
// -----------------------------------------------------------------
