// -----------------------------------------------------------------
#ifndef _RANDOM_H
#define _RANDOM_H
// Random number structure and generic random number generator
// returning a uniformly distributed random value on [0, 1]
// We assume long is at least 32 bits
#include "../include/precision.h"

typedef struct {
  unsigned long r0, r1, r2, r3, r4, r5, r6;
  unsigned long multiplier, addend, ic_state;
  float scale;
} double_prn;

Real myrand(double_prn *prn_pt);

#endif
// -----------------------------------------------------------------
