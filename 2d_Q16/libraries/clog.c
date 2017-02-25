// -----------------------------------------------------------------
// Return complex logarithm of given complex number
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

complex clog(complex *a) {
  complex c;
  c.real = 0.5 * (Real)log((double)((*a).real * (*a).real
                                  + (*a).imag * (*a).imag));
  c.imag = (Real)atan2((double)(*a).imag, (double)(*a).real);
  return c;
}
// -----------------------------------------------------------------
