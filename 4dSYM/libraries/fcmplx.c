// -----------------------------------------------------------------
// Return single-precision complex number
// from two given generic-precision real numbers
#include "../include/config.h"
#include "../include/complex.h"

fcomplex fcmplx(Real x, Real y) {
  fcomplex c;
  c.real = x;
  c.imag = y;
  return c;
}
// -----------------------------------------------------------------
