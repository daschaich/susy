// -----------------------------------------------------------------
// Clear a fundamental vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clearvec_f(su3_vector_f *v) {
  int i;
  for (i = 0; i < NCOL; i++) {
    v->c[i].real = 0.0;
    v->c[i].imag = 0.0;
  }
}
// -----------------------------------------------------------------
