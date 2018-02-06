// -----------------------------------------------------------------
// Average and extremal scalar eigenvalues, and widths of distributions
// Use LAPACK for arbitrary NCOL
// Note two-color trace subtraction --> +/- eigenvalue pairs

// #define SCALAR_EIG_DIST prints out all eigenvalues for plotting distribution
// CAUTION: Do not run SCALAR_EIG_DIST with MPI!

//#define SCALAR_EIG_DIST
#include "susy_includes.h"

void scalar_eig(double *ave_eigs, double *eig_widths,
                double *min_eigs, double *max_eigs) {

  register int i;
  register site *s;
  char N = 'N';     // Ask LAPACK only for eigenvalues
  char U = 'U';     // Have LAPACK store upper triangle of U.Ubar
  int row, col, Npt = NCOL, stat = 0, Nwork = 2 * NCOL, j, k;
  double sq_eigs[NCOL], norm = 1.0 / (double)(NSCALAR * nt);

#ifdef SCALAR_EIG_DIST
  if (this_node != 0) {
    printf("scalar_eig: don't run SCALAR_EIG_DIST in parallel\n");
    fflush(stdout);
    terminate(1);
  }
#endif

  // Initialize averages and extrema
  for (j = 0; j < NCOL; j++) {
    ave_eigs[j] = 0.0;
    min_eigs[j] =  99.0;
    max_eigs[j] = -99.0;
    sq_eigs[j] = 0.0;
  }

  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
      // Convert X[j] to column-major double array used by LAPACK
      for (row = 0; row < NCOL; row++) {
        for (col = 0; col < NCOL; col++) {
          store[2 * (col * NCOL + row)] = s->X[j].e[row][col].real;
          store[2 * (col * NCOL + row) + 1] = s->X[j].e[row][col].imag;
        }
      }

      // Compute eigenvalues and eigenvectors of X[j]
      zheev_(&N, &U, &Npt, store, &Npt, eigs, work, &Nwork, Rwork, &stat);

      // Make sure eigenvalues are always ordered consistently
      if (stat != 0)
        printf("WARNING: Non-zero return for %d\n", s->t);

#ifdef SCALAR_EIG_DIST
      printf("SCALAR EIG DIST ");
      printf("%d %d", s->t, j);
      for (k = 0; k < NCOL; k++)
        printf(" %.4g", eigs[k]);
      printf("\n");
#endif

      // Average eigenvalues, monitor minimum and maximum
      for (k = 0; k < NCOL; k++) {
        ave_eigs[k] += eigs[k];
        sq_eigs[k] += eigs[k] * eigs[k];
        if (eigs[k] > max_eigs[k])
          max_eigs[k] = eigs[k];
        if (eigs[k] < min_eigs[k])
          min_eigs[k] = eigs[k];
      }
    }
  }

  // Finalize averages, extrema and square root of variance
  for (j = 0; j < NCOL; j++) {
    g_doublesum(&(ave_eigs[j]));
    ave_eigs[j] *= norm;

    g_doublesum(&(sq_eigs[j]));
    sq_eigs[j] *= norm;

    eig_widths[j] = sqrt(sq_eigs[j] - ave_eigs[j] * ave_eigs[j]);

    g_doublemax(&(max_eigs[j]));

    min_eigs[j] = -min_eigs[j];
    g_doublemax(&(min_eigs[j]));
    min_eigs[j] = -min_eigs[j];
  }
}
// -----------------------------------------------------------------
