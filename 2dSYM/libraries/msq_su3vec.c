// -----------------------------------------------------------------
// Squared magnitude of irrep vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
Real magsq_su3vec(su3_vector *a) {
  register Real sum = 0.0;
  register int i;
  for (i = 0; i < DIMF; i++)
    sum += a->c[i].real * a->c[i].real
         + a->c[i].imag * a->c[i].imag;
  return sum;
}
#else   // FAST version for DIMF=3 only
Real magsq_su3vec(su3_vector *a) {
  register Real temp, sum = 0.0;
  temp = a->c[0].real * a->c[0].real;
  sum += temp;
  temp = a->c[0].imag * a->c[0].imag;
  sum += temp;
  temp = a->c[1].real * a->c[1].real;
  sum += temp;
  temp = a->c[1].imag * a->c[1].imag;
  sum += temp;
  temp = a->c[2].real * a->c[2].real;
  sum += temp;
  temp = a->c[2].imag * a->c[2].imag;
  sum += temp;
  return sum;
}
#endif
// -----------------------------------------------------------------
