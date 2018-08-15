// -----------------------------------------------------------------
// Check unitarity of the link matrices
// Terminate if deviations exceed TOLERANCE
#include "generic_includes.h"

#define TOLERANCE 0.0001
#define STRONG    // Check row orthogonality as well as norms
//#define UNIDEBUG
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real check_unit(matrix *c) {
  register int i, j;
  register Real ar, max = 0.0;
#ifdef STRONG
  register int k;
  complex sum, tc;
#endif

  // Check normalization of each row
  for (i = 0; i < NCOL; i++) {
    ar = 0.0;
    for (j = 0; j < NCOL; j++) {    // Sum of squares of row
      ar += (*c).e[i][j].real * (*c).e[i][j].real
          + (*c).e[i][j].imag * (*c).e[i][j].imag;
    }
    ar = fabs(sqrt((double)ar) - 1.0);
    if (max < ar)
      max = ar;
  }

#ifdef STRONG
  // Test orthogonality of row i and row j
  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      CMULJ_(c->e[j][0], c->e[i][0], sum);
      for (k = 1; k < NCOL; k++) {
        CMULJ_(c->e[j][k], c->e[i][k], tc);
        CSUM(sum, tc);
      }
      ar = cabs(&sum);
      if (max < ar)
        max = ar;
    }
  }
#endif

  return max;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real check_unitarity() {
  register int i;
  register site *s;
  register matrix *mat;
  Real deviation, max_deviation = 0.0;
  double av_deviation = 0.0;

  FORALLSITES(i, s) {
    mat = &(s->link);
    deviation = check_unit(mat);
    if (deviation > TOLERANCE) {
      printf("Unitarity problem on node %d, site %d, deviation=%f\n",
             mynode(), i, deviation);
      printf("SU(N) matrix:\n");
      dumpmat(mat);
      fflush(stdout);
      terminate(1);
    }
    if (max_deviation < deviation)
      max_deviation = deviation;
    av_deviation += deviation * deviation;
  }

  av_deviation = sqrt(av_deviation / (double)nt);
#ifdef UNIDEBUG
  printf("Deviation from unitarity on node %d: max %.4g, ave %.4g\n",
         mynode(), max_deviation, av_deviation);
#endif
  if (max_deviation > TOLERANCE) {
    printf("Unitarity problem on node %d, maximum deviation %.4g\n",
           mynode(), max_deviation);
  }
  return max_deviation;
}
// -----------------------------------------------------------------
