// -----------------------------------------------------------------
// Copy an irrep vector (hardly worth a function)
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void su3vec_copy(su3_vector *src, su3_vector *dest) {
  *dest = *src;
}
// -----------------------------------------------------------------
