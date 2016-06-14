// -----------------------------------------------------------------
// Return double complex exponential of given double complex number
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

double_complex dcexp(double_complex *a) {
  double_complex c;
  double mag = (double)exp((double)(*a).real);
  c.real = mag * (double)cos((double)(*a).imag);
  c.imag = mag * (double)sin((double)(*a).imag);
  return c;
}
// -----------------------------------------------------------------
