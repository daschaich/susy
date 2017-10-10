// -----------------------------------------------------------------
// C random number generator for parallel processors
// Exclusive or of feedback shift register and integer congruence generator
// Use a different multiplier on each generator
// Make sure that fsr is initialized differently on each generator

// Usage:
//   int seed;
//   Real x;
//   initialize_prn(prn_pt, seed, index)
//   x = myrand(prn_pt);

// prn_pt is address of a struct double_prn, index selects algorithm
#include "generic_includes.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// index selects which random number generator -- which multiplier
void initialize_prn(double_prn *prn_pt, int seed, int index) {
  seed = (69607 + 8 * index) * seed + 12345;
  prn_pt->r0 = (seed>>8) & 0xffffff;
  seed = (69607 + 8 * index) * seed + 12345;
  prn_pt->r1 = (seed>>8) & 0xffffff;
  seed = (69607 + 8 * index) * seed + 12345;
  prn_pt->r2 = (seed>>8) & 0xffffff;
  seed = (69607 + 8 * index) * seed + 12345;
  prn_pt->r3 = (seed>>8) & 0xffffff;
  seed = (69607 + 8 * index) * seed + 12345;
  prn_pt->r4 = (seed>>8) & 0xffffff;
  seed = (69607 + 8 * index) * seed + 12345;
  prn_pt->r5 = (seed>>8) & 0xffffff;
  seed = (69607 + 8 * index) * seed + 12345;
  prn_pt->r6 = (seed>>8) & 0xffffff;
  seed = (69607 + 8 * index) * seed + 12345;
  prn_pt->ic_state = seed;
  prn_pt->multiplier = 100005 + 8 * index;
  prn_pt->addend = 12345;
  prn_pt->scale = 1.0 / ((Real)0x1000000);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real myrand(double_prn *prn_pt) {
  register int t, s;

  t = (((prn_pt->r5 >> 7) | (prn_pt->r6 << 17)) ^
       ((prn_pt->r4 >> 1) | (prn_pt->r5 << 23))) & 0xffffff;
  prn_pt->r6 = prn_pt->r5;
  prn_pt->r5 = prn_pt->r4;
  prn_pt->r4 = prn_pt->r3;
  prn_pt->r3 = prn_pt->r2;
  prn_pt->r2 = prn_pt->r1;
  prn_pt->r1 = prn_pt->r0;
  prn_pt->r0 = t;
  s = prn_pt->ic_state * prn_pt->multiplier + prn_pt->addend;
  prn_pt->ic_state = s;
  return (prn_pt->scale*(t ^ ((s>>8)&0xffffff)));
}
// -----------------------------------------------------------------
