// -----------------------------------------------------------------
// Clear a fundamental vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void clearvec(vector *v) {
  int i;
  for (i = 0; i < NCOL; i++) {
    v->c[i].real = 0.0;
    v->c[i].imag = 0.0;
  }
}
// -----------------------------------------------------------------
