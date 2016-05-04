// -----------------------------------------------------------------
// Deconstruct NxN fermion matrices into N^2 vectors
// psi^a = -Tr[psi.Lambda^a]
#include "susy_includes.h"

#if (DIMF != NCOL * NCOL)
  #error "Fermion reconstruction only for adjoint rep!"
#endif

void deconstruct(matrix *in, vector *out) {
  register int a;
  register complex tc = complextrace_nn(in, &(Lambda[0]));
  out->c[0].real = tc.real;
  out->c[0].imag = tc.imag;
  for (a = 1; a < DIMF; a++) {
    tc = complextrace_nn(in, &(Lambda[a]));
    out->c[a].real = tc.real;
    out->c[a].imag = tc.imag;
  }
}
// -----------------------------------------------------------------
