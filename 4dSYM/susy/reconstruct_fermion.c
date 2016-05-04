// -----------------------------------------------------------------
// Reconstruct NxN fermion matrices from N^2 vectors
#include "susy_includes.h"

#if (DIMF != NCOL * NCOL)
  #error "Fermion reconstruction only for adjoint rep!"
#endif

void reconstruct(vector *in, matrix *out) {
  register int a, i, j;

  // Overwrite out
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      out->e[i][j].real = Lambda[0].e[i][j].real * in->c[0].real
                        - Lambda[0].e[i][j].imag * in->c[0].imag;
      out->e[i][j].imag = Lambda[0].e[i][j].imag * in->c[0].real
                        + Lambda[0].e[i][j].real * in->c[0].imag;
    }
  }
  for (a = 1; a < DIMF; a++) {
    for (i = 0; i < NCOL; i++) {
      for (j = 0; j < NCOL; j++) {
        out->e[i][j].real += Lambda[a].e[i][j].real * in->c[a].real
                           - Lambda[a].e[i][j].imag * in->c[a].imag;
        out->e[i][j].imag += Lambda[a].e[i][j].imag * in->c[a].real
                           + Lambda[a].e[i][j].real * in->c[a].imag;
      }
    }
  }
}

// Negate imaginary parts of in-->in^*
void reconstruct_star(vector *in, matrix *out) {
  register int a, i, j;

  // Overwrite out
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      out->e[i][j].real = Lambda[0].e[i][j].real * in->c[0].real
                        + Lambda[0].e[i][j].imag * in->c[0].imag;
      out->e[i][j].imag = Lambda[0].e[i][j].imag * in->c[0].real
                        - Lambda[0].e[i][j].real * in->c[0].imag;
    }
  }
  for (a = 1; a < DIMF; a++) {
    for (i = 0; i < NCOL; i++) {
      for (j = 0; j < NCOL; j++) {
        out->e[i][j].real += Lambda[a].e[i][j].real * in->c[a].real
                           + Lambda[a].e[i][j].imag * in->c[a].imag;
        out->e[i][j].imag += Lambda[a].e[i][j].imag * in->c[a].real
                           - Lambda[a].e[i][j].real * in->c[a].imag;
      }
    }
  }
}
// -----------------------------------------------------------------
