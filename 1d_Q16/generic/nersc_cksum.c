// -----------------------------------------------------------------
// Computes the mean global sum of the trace of the gauge links
// Used to aid checking lattice file integrity
// chksum itself is set to zero
#include "generic_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
u_int32type nersc_cksum() {
  u_int32type chksum = 0;
  return chksum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void sum_linktr(double_complex *linktrsum) {
  int i;
  site *s;
  matrix *a;

  linktrsum->real = 0.0;
  linktrsum->imag = 0.0;

  FORALLSITES(i, s) {
    a = &(s->link);
    CSUM(*linktrsum, a->e[0][0]);
    CSUM(*linktrsum, a->e[1][1]);
#if (NCOL > 2)
    CSUM(*linktrsum, a->e[2][2]);
#if (NCOL > 3)
    CSUM(*linktrsum, a->e[3][3]);
#if (NCOL > 4)
    int j;
    for (j = 4; j < NCOL; j++)
      CSUM(*linktrsum, a->e[j][j]);
#endif
#endif
#endif
  }

  g_dcomplexsum(linktrsum);
  CDIVREAL(*linktrsum, nt, *linktrsum);
}
// -----------------------------------------------------------------
