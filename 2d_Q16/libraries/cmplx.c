// -----------------------------------------------------------------
// Return complex number from two given real numbers
#include "../include/config.h"
#include "../include/complex.h"

complex cmplx(Real x, Real y) {
  complex c;
  c.real = x;
  c.imag = y;
  return c;
}
// -----------------------------------------------------------------
