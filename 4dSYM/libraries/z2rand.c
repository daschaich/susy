// -----------------------------------------------------------------
// Z2 distributed random number
// Normalized by 1 / sqrt(2) so that x^2 = 1 / 2
// This requires a random number generator, myrand(),
// to return a Real uniformly distributed between zero and one
#include "../include/config.h"
#include <math.h>
#include "../include/susy.h"
#include "../include/random.h"

// prn_pt is a pointer passed to myrand()
Real Z2_rand_no(double_prn *prn_pt) {
  if (myrand(prn_pt) > 0.5)
    return (1.0 / sqrt(2.0));
  return (-1.0 / sqrt(2.0));
}
// -----------------------------------------------------------------
