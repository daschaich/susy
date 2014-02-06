// -----------------------------------------------------------------
// Copy an irrep matrix (hardly worth a function)
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void su3mat_copy(su3_matrix *src, su3_matrix *dest) {
  *dest = *src;
}
// -----------------------------------------------------------------
