// -----------------------------------------------------------------
// Return exp[i theta], probably sub-optimal
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

complex ce_itheta(Real theta) {
  complex c;
  c.real = (Real)cos((double)theta);
  c.imag = (Real)sin((double)theta);
  return c;
}
// -----------------------------------------------------------------
