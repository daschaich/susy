// -----------------------------------------------------------------
// Return double exp[i theta], probably sub-optimal
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

double_complex dce_itheta(double theta) {
  double_complex c;
  c.real = (double)cos((double)theta);
  c.imag = (double)sin((double)theta);
  return c;
}
// -----------------------------------------------------------------
