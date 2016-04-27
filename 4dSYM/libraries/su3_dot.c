// -----------------------------------------------------------------
// Return the dot product of two su3_vectors: adag b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
complex su3_dot(su3_vector *a, su3_vector *b) {
#ifndef FAST
  complex temp1, temp2;
  CMULJ_(a->c[0], b->c[0], temp1);
  CMULJ_(a->c[1], b->c[1], temp2);
  CSUM(temp1, temp2);
#if (DIMF > 2)
  CMULJ_(a->c[2], b->c[2], temp2);
  CSUM(temp1, temp2);
#if (DIMF > 3)
  CMULJ_(a->c[3], b->c[3], temp2);
  CSUM(temp1, temp2);
#if (DIMF > 4)
  CMULJ_(a->c[4], b->c[4], temp2);
  CSUM(temp1, temp2);
#if (DIMF > 5)
  CMULJ_(a->c[5], b->c[5], temp2);
  CSUM(temp1, temp2);
#if (DIMF > 6)
  CMULJ_(a->c[6], b->c[6], temp2);
  CSUM(temp1, temp2);
#if (DIMF > 7)
  register int i;
  for (i = 7; i < DIMF; i++) {
    CMULJ_(a->c[i], b->c[i], temp2);
    CSUM(temp1, temp2);
  }
#endif
#endif
#endif
#endif
#endif
#endif
  return temp1;
#else // FAST version for DIMF=3 only
  register Real ar, ai, br, bi, cr, ci;
  register complex cc;

  ar = a->c[0].real;
  ai = a->c[0].imag;
  br = b->c[0].real;
  bi = b->c[0].imag;
  cr = ar * br + ai * bi;
  ci = ar * bi - ai * br;

  ar = a->c[1].real;
  ai = a->c[1].imag;
  br = b->c[1].real;
  bi = b->c[1].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;

  ar = a->c[2].real;
  ai = a->c[2].imag;
  br = b->c[2].real;
  bi = b->c[2].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;

  cc.real = cr;
  cc.imag = ci;
  return cc;
#endif
}
// -----------------------------------------------------------------
