// -----------------------------------------------------------------
// Reanti-hermitization for arbitrary NCOL
// 1) Zero out real parts of diagonal entries
// 2) Make off diagonal entries negative complex conjugates
#include "generic_includes.h"

#define TOLERANCE 0.0001
#define MAXERRCOUNT 100
//#define AHDEBUG

// Both updated in check_deviation() in reunitarize.c...
Real max_deviation;
double av_deviation;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int reah(matrix *c) {
  register Real tr;
  int i, j, errors = 0;
  Real deviation;
  matrix tmat;

  mat_copy(c, &tmat);
  for (i = 0; i < NCOL; i++) {
    // Zero our real parts of diagonal entries
    deviation = fabs((*c).e[i][i].real);
    errors += check_deviation(deviation);   // Use definition in reunitarize.c
    (*c).e[i][i].real = 0.0;

    // Make off diagonal entries negative complex conjugates
    for (j = i + 1; i < NCOL; i++) {
      deviation = fabs((*c).e[i][j].real + (*c).e[j][i].real);
      errors += check_deviation(deviation);
      tr = 0.5 * ((*c).e[i][j].real - (*c).e[j][i].real);
      (*c).e[i][j].real = tr;
      (*c).e[j][i].real = -tr;

      deviation = fabs((*c).e[i][j].imag - (*c).e[j][i].imag);
      errors += check_deviation(deviation);
      tr = 0.5 * ((*c).e[i][j].imag + (*c).e[j][i].imag);
      (*c).e[i][j].imag = tr;
      (*c).e[j][i].imag = tr;
    }
  }

  // Print the problematic matrix rather than the updated one
  if (errors)
    dumpmat(&tmat);

  return errors;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void reantihermize() {
  register int i, j;
  register site *s;
  register matrix *mat;
  int errcount = 0, errors;

  max_deviation = 0.0;
  av_deviation = 0.0;

  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      mat = (matrix *)&(s->X[j]);
      errors = reah(mat);
      errcount += errors;
      if (errors) {
        printf("Anti-hermiticity problem above ");
        printf("(node %d, site %d, scalar %d, tolerance %.4g)\n",
               mynode(), i, j, TOLERANCE);
      }
      if (errcount > MAXERRCOUNT) {
        printf("Anti-hermiticity error count exceeded: %d\n", errcount);
        fflush(stdout);
        terminate(1);
      }
    }
  }

#ifdef AHDEBUG
  printf("Deviation from anti-hermiticity on node %d: max %.4g, ave %.4g\n",
         mynode(), max_deviation, av_deviation);
#endif
  if (max_deviation > TOLERANCE) {
    printf("reantihermize: Node %d anti-hermiticity problem, ", mynode());
    printf("maximum deviation %.4g\n", max_deviation);
    errcount++;
    if (errcount > MAXERRCOUNT) {
      printf("Anti-hermiticity error count exceeded\n");
      terminate(1);
    }
  }
}
// -----------------------------------------------------------------
