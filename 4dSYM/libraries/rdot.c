// -----------------------------------------------------------------
// Real part of dot product of two irrep vectors
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

Real rdot(vector *a, vector *b) {
  register Real temp1,temp2;
  temp2 = a->c[0].real * b->c[0].real;
  temp1 = a->c[0].imag * b->c[0].imag;
  temp2 += temp1;
#if (DIMF<=6)
  temp1 = a->c[1].real * b->c[1].real;
  temp2 += temp1;
  temp1 = a->c[1].imag * b->c[1].imag;
  temp2 += temp1;
#if (DIMF>2)
  temp1 = a->c[2].real * b->c[2].real;
  temp2 += temp1;
  temp1 = a->c[2].imag * b->c[2].imag;
  temp2 += temp1;
#if (DIMF>3)
  temp1 = a->c[3].real * b->c[3].real;
  temp2 += temp1;
  temp1 = a->c[3].imag * b->c[3].imag;
  temp2 += temp1;
#if (DIMF>4)
  temp1 = a->c[4].real * b->c[4].real;
  temp2 += temp1;
  temp1 = a->c[4].imag * b->c[4].imag;
  temp2 += temp1;
#if (DIMF>5)
  temp1 = a->c[5].real * b->c[5].real;
  temp2 += temp1;
  temp1 = a->c[5].imag * b->c[5].imag;
  temp2 += temp1;
#endif
#endif
#endif
#endif
#else // DIMF>6
  int i;
  for (i = 1; i < DIMF; i++) {
    temp1 = a->c[i].real * b->c[i].real;
    temp2 += temp1;
    temp1 = a->c[i].imag * b->c[i].imag;
    temp2 += temp1;
  }
#endif
  return temp2;
}
// -----------------------------------------------------------------
