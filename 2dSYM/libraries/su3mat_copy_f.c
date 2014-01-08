// -----------------------------------------------------------------
// Copy a fundamental matrix (hardly worth a function)
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void su3mat_copy_f(su3_matrix_f *src, su3_matrix_f *dest) {
  *dest = *src;
}
// -----------------------------------------------------------------
