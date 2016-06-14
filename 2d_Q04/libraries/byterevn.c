// -----------------------------------------------------------------
// Byte reversal routines
// WARNING: Must be compiled with appropriate short32 flag!
#include "../include/config.h"
#include "../include/int32type.h"
#include <assert.h>
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Byte reversal on n contiguous 32-bit words
void byterevn(int32type w[], int n) {
  register int32type old,newv;
  int j;

  assert(sizeof(int32type) == 4);
  for (j = 0; j < n; j++) {
    old = w[j];
    newv = old >> 24 & 0x000000ff;
    newv |= old >> 8 & 0x0000ff00;
    newv |= old << 8 & 0x00ff0000;
    newv |= old << 24 & 0xff000000;
    w[j] = newv;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Byte reversal on n contiguous 32-bit words
void byterevn64(int32type w[], int n) {
  int32type tmp;
  int j;

  assert(sizeof(int32type) == 4);
  // First swap pairs of 32-bit words
  for (j = 0; j < n; j++) {
    tmp = w[2 * j];
    w[2 * j] = w[2 * j + 1];
    w[2 * j + 1] = tmp;
  }

  // Then swap bytes in 32-bit words
  byterevn(w, 2 * n);
}
// -----------------------------------------------------------------
