// -----------------------------------------------------------------
// Print the given anti-hermitian matrix
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dump_ahmat(anti_hermitmat *ahm) {
  int i, j, index = 0;
  printf("DIAG");
  for (i = 0; i < NCOL; i++)
    printf(" %.4g", ahm->im_diag[i]);
  printf("\nOFFDIAG");
  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      printf(" (%.4g, %.4g)", ahm->m[index].real, ahm->m[index].imag);
      index++;
    }
  }
  printf("\n\n");
}
// -----------------------------------------------------------------
