// -----------------------------------------------------------------
// Gaussian distributed random number
// Probability distribution exp(-x * x), so <x^2> = 1 / 2
// This requires a random number generator, myrand(),
// to return a Real uniformly distributed between zero and one
#include "../include/config.h"
#include <math.h>
#include "../include/su3.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// prn_pt is a pointer passed to myrand()
Real gaussian_rand_no(double_prn *prn_pt) {
  static int iset = 0;
  static Real gset;
  Real fac, r, v1, v2;

  if  (iset == 0) {
    do {
      v1 = 2.0 * myrand(prn_pt) - 1.0;    // [-1, 1]
      v2 = 2.0 * myrand(prn_pt) - 1.0;
      r = v1 * v1 + v2 * v2;
    } while (r >= 1.0);
    fac = sqrt(-log((double)r) / (double)r);
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  }
  else {
    iset = 0;
    return gset;
  }
}
// -----------------------------------------------------------------
