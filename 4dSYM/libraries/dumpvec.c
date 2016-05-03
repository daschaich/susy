// -----------------------------------------------------------------
// Print the given irrep vector
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/susy.h"

void dumpvec(vector *v) {
  int j;
  printf("(%.4g, %.4g)", v->c[0].real, v->c[0].imag);
  for (j = 1; j < DIMF; j++)
    printf("  (%.4g, %.4g)", v->c[j].real, v->c[j].imag);
  printf("\n");
}
// -----------------------------------------------------------------
