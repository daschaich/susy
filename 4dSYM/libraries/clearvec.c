// -----------------------------------------------------------------
// Clear an irrep vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clearvec(su3_vector *v) {
#if (DIMF<=6)
  v->c[0].real = 0.0;
  v->c[0].imag = 0.0;
  v->c[1].real = 0.0;
  v->c[1].imag = 0.0;
#if (DIMF>2)
  v->c[2].real = 0.0;
  v->c[2].imag = 0.0;
#if (DIMF>3)
  v->c[3].real = 0.0;
  v->c[3].imag = 0.0;
#if (DIMF>4)
  v->c[4].real = 0.0;
  v->c[4].imag = 0.0;
#if (DIMF>5)
  v->c[5].real = 0.0;
  v->c[5].imag = 0.0;
#endif
#endif
#endif
#endif
#else  // DIMF > 6
  int i;
  for (i = 0; i < DIMF; i++) {
    v->c[i].real = 0.0;
    v->c[i].imag = 0.0;
  }
#endif
}
// -----------------------------------------------------------------
