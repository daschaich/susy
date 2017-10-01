// -----------------------------------------------------------------
// Return double complex number from two given double real numbers
#include "../include/config.h"
#include "../include/complex.h"

double_complex dcmplx(double x, double y) {
  double_complex c;
  c.real = x;
  c.imag = y;
  return c;
}
// -----------------------------------------------------------------
