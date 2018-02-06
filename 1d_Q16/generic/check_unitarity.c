// -----------------------------------------------------------------
// Check unitarity of the link matrices, terminate if not unitary
#include "generic_includes.h"

#define TOLERANCE 0.0001
#define STRONG    // Check row orthogonality as well as norms
//#define UNIDEBUG
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real check_unit(matrix *c) {
  register int i, j, k;
  register Real ar, ai, ari, max = 0.0;

  // Check normalization of each row
  for (i = 0; i < NCOL; i++) {
    ar = 0.0;
    for (j=0;j<NCOL;j++) {
      ar += (*c).e[i][j].real * (*c).e[i][j].real +    /* sum of squares of row */
  (*c).e[i][j].imag * (*c).e[i][j].imag;
    }
    ar =  fabs(sqrt((double)ar) - 1.0);
    if (max < ar)
      max = ar;
  }

#ifdef STRONG
  // Test orthogonality of row i and row j
  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      ar = 0.0;   // Real part of i dot j
      ai = 0.0;   // Imag part of i dot j
      for (k = 0; k < NCOL; k++) {
        ar += (*c).e[i][k].real * (*c).e[j][k].real
            + (*c).e[i][k].imag * (*c).e[j][k].imag;
        ai += (*c).e[i][k].real * (*c).e[j][k].imag
            - (*c).e[i][k].imag * (*c).e[j][k].real;
      }
      ari = sqrt((double)(ar * ar + ai * ai));
      if (max < ari)
        max = ari;
    }
  }
#endif

  return max;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real check_unitarity() {
  register int i,j;
  int ii, jj;
  register site *s;
  register matrix *mat;
  Real deviation, max_deviation = 0.0;
  double av_deviation = 0.0;
  union {
    Real fval;
    int ival;
  } ifval;
  
  FORALLSITES(i, s) {
    mat = (matrix *)&(s->link);
    deviation = check_unit(mat);
    if (deviation > TOLERANCE) {
      printf("Unitarity problem on node %d, site %d, gauge deviation=%f\n",
             mynode(), i, deviation);
      printf("SU(N) matrix:\n");
      for (ii = 0; ii < NCOL; ii++) {
        for (jj = 0; jj < NCOL; jj++) {
          printf("%f ", (*mat).e[ii][jj].real);
          printf("%f ", (*mat).e[ii][jj].imag);
        }
        printf("\n");
      }
      printf("repeat in hex:\n");
      for (ii = 0; ii < NCOL; ii++) {
        for (jj = 0; jj < NCOL; jj++) {
          ifval.fval = (*mat).e[ii][jj].real;
          printf("%08x ", ifval.ival);
          ifval.fval = (*mat).e[ii][jj].imag;
          printf("%08x ", ifval.ival);
        }
        printf("\n");
      }
      printf("  \n\n");
      fflush(stdout);
      terminate(1);
    }
    if (max_deviation < deviation)
      max_deviation = deviation;
    av_deviation += deviation * deviation;
    for(j=0;j<NSCALAR;j++) {
      mat = (matrix *)&(s->X[j]);
      deviation = check_unit(mat);
      if (deviation > TOLERANCE) {
        printf("Unitarity problem on node %d, site %d, scalar %d, deviation=%f\n",
               mynode(), i, j, deviation);
        printf("SU(N) matrix:\n");
        for (ii = 0; ii < NCOL; ii++) {
          for (jj = 0; jj < NCOL; jj++) {
            printf("%f ", (*mat).e[ii][jj].real);
            printf("%f ", (*mat).e[ii][jj].imag);
          }
          printf("\n");
        }
        printf("repeat in hex:\n");
        for (ii = 0; ii < NCOL; ii++) {
          for (jj = 0; jj < NCOL; jj++) {
            ifval.fval = (*mat).e[ii][jj].real;
            printf("%08x ", ifval.ival);
            ifval.fval = (*mat).e[ii][jj].imag;
            printf("%08x ", ifval.ival);
          }
          printf("\n");
        }
        printf("  \n\n");
        fflush(stdout);
        terminate(1);
      }
      if (max_deviation < deviation)
        max_deviation = deviation;
      av_deviation += deviation * deviation;
    }
  }


  av_deviation = sqrt(av_deviation / ((NSCALAR + 1.0) * nt));
#ifdef UNIDEBUG
  printf("Deviation from unitarity on node %d: max %.4g, ave %.4g\n",
         mynode(), max_deviation, av_deviation);
#endif
  if (max_deviation > TOLERANCE)
    printf("Unitarity problem on node %d, maximum deviation %.4g\n",
           mynode(), max_deviation);
  return max_deviation;
}
// -----------------------------------------------------------------
