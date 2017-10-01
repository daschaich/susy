// -----------------------------------------------------------------
// Average and extremal scalar eigenvalues, and widths of distributions
// Use LAPACK for arbitrary NCOL
// Note two-color trace subtraction --> +/- eigenvalue pairs

// #define SCALAR_EIG_DIST prints out all eigenvalues for plotting distribution
// CAUTION: Do not run SCALAR_EIG_DIST with MPI!

//#define SCALAR_EIG_DIST
#include "susy_includes.h"

// project == 1 tells us to consider scalar fields from polar decomposition
// rather than scalar fields from the U.Ubar product
void scalar_eig(int project, double *ave_eigs, double *eig_widths,
                double *min_eigs, double *max_eigs) {

  register int i, dir;
  register site *s;
  char N = 'N';     // Ask LAPACK only for eigenvalues
  char U = 'U';     // Have LAPACK store upper triangle of U.Ubar
  int row, col, Npt = NCOL, stat = 0, Nwork = 2 * NCOL, j;
  Real tr;
  double norm = NUMLINK * volume, sq_eigs[NCOL];
  complex tc;
  matrix USq, tmat;

#ifdef SCALAR_EIG_DIST
  if (this_node != 0) {
    printf("plaquette: don't run SCALAR_EIG_DIST in parallel\n");
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
    FORALLDIR(dir) {
      if (project == 1) {   // Consider polar-projected scalar fields
        polar(&(s->link[dir]), &USq, &tmat);
        // Take log
        matrix_log(&tmat, &USq);
      }
      else {                // Consider U.Ubar scalar fields
        mult_na(&(s->link[dir]), &(s->link[dir]), &USq);
        // Take traceless part
        tc = trace(&USq);
        tr = one_ov_N * tc.real;
        for (j = 0; j < NCOL; j++)
          USq.e[j][j].real -= tr;
      }

      // Convert USq to column-major double array used by LAPACK
      for (row = 0; row < NCOL; row++) {
        for (col = 0; col < NCOL; col++) {
          store[2 * (col * NCOL + row)] = USq.e[row][col].real;
          store[2 * (col * NCOL + row) + 1] = USq.e[row][col].imag;
        }
      }

      // Compute eigenvalues and eigenvectors of USq
      zheev_(&N, &U, &Npt, store, &Npt, eigs, work, &Nwork, Rwork, &stat);

      // Make sure eigenvalues are always ordered consistently
      if (stat != 0) {
        printf("WARNING: Non-zero return for %d %d %d %d %d\n",
               s->x, s->y, s->z, s->t, dir);
      }

#ifdef SCALAR_EIG_DIST
      if (project == 1)
        printf("POLAR_EIG_DIST ");
      else
        printf("UUBAR_EIG_DIST ");
      printf("%d %d %d %d %d", s->x, s->y, s->z, s->t, dir);
      for (j = 0; j < NCOL; j++)
        printf(" %.4g", eigs[j]);
      printf("\n");
#endif

      // Average eigenvalues, monitor minimum and maximum
      for (j = 0; j < NCOL; j++) {
        ave_eigs[j] += eigs[j];
        sq_eigs[j] += eigs[j] * eigs[j];
        if (eigs[j] > max_eigs[j])
          max_eigs[j] = eigs[j];
        if (eigs[j] < min_eigs[j])
          min_eigs[j] = eigs[j];
      }
    }
  }

  // Finalize averages, extrema and square root of variance
  for (j = 0; j < NCOL; j++) {
    g_doublesum(&(ave_eigs[j]));
    ave_eigs[j] /= norm;

    g_doublesum(&(sq_eigs[j]));
    sq_eigs[j] /= norm;

    eig_widths[j] = sqrt(sq_eigs[j] - ave_eigs[j] * ave_eigs[j]);

    g_doublemax(&(max_eigs[j]));

    min_eigs[j] = -min_eigs[j];
    g_doublemax(&(min_eigs[j]));
    min_eigs[j] = -min_eigs[j];
  }
}
// -----------------------------------------------------------------
