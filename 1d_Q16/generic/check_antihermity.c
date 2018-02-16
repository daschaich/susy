// -----------------------------------------------------------------
// Check anti-hermiticity of the scalar matrices
// Terminate if deviations exceed TOLERANCE
#include "generic_includes.h"

#define TOLERANCE 0.0001
//#define AHDEBUG
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real check_ah(matrix *c) {
  register int i, j;
  register Real ar, ai, ari, max = 0.0;

  // Check real part of diagonal entries
  for (i = 0; i < NCOL; i++) {
    ar = fabs((*c).e[i][i].real);
    if (max < ar)
      max = ar;
  }

  // Check that off-diagonal entries are negative complex conjugates
  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      ar = (*c).e[i][j].real + (*c).e[j][i].real;
      ai = (*c).e[i][j].imag - (*c).e[j][i].imag;
      ari = sqrt(ar * ar + ai * ai);
      if (max < ari)
        max = ari;
    }
  }

  return max;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real check_antihermity() {
  register int i, j;
  register site *s;
  register matrix *mat;
  Real deviation, max_deviation = 0.0;
  double av_deviation = 0.0;

  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      mat = (matrix *)&(s->X[j]);
      deviation = check_ah(mat);
      if (deviation > TOLERANCE) {
        printf("Anti-hermiticity problem on node %d, site %d, ", mynode(), i);
        printf("scalar %d, deviation=%f\n", j, deviation);
        printf("Matrix:\n");
        dumpmat(mat);
        fflush(stdout);
        terminate(1);
      }
      if (max_deviation < deviation)
        max_deviation = deviation;
      av_deviation += deviation * deviation;
    }
  }

  av_deviation = sqrt(av_deviation / (double)(NSCALAR * nt));
#ifdef AHDEBUG
  printf("Deviation from anti-hermiticity on node %d: max %.4g, ave %.4g\n",
         mynode(), max_deviation, av_deviation);
#endif
  if (max_deviation > TOLERANCE) {
    printf("Anti-hermiticity problem on node %d, maximum deviation %.4g\n",
           mynode(), max_deviation);
  }
  return max_deviation;
}
// -----------------------------------------------------------------
