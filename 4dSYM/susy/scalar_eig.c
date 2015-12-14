// -----------------------------------------------------------------
// Average and extremal scalar eigenvalues, and widths of distributions
// Use LAPACK for arbitrary NCOL
// Note two-color trace subtraction --> +/- eigenvalue pairs

// #define SCALAR_EIG_DIST prints out all eigenvalues for plotting distribution
// CAUTION: Do not run PLAQ_DIST with MPI!

//#define SCALAR_EIG_DIST
#include "susy_includes.h"

// project == 1 tells us to consider scalar fields from polar decomposition
// rather than scalar fields from the U.Ubar product
void scalar_eig(int project) {
  register int i, dir;
  register site *s;
  char N = 'N';     // Ask LAPACK only for eigenvalues
  char U = 'U';     // Have LAPACK store upper triangle of U.Ubar
  int row, col, Npt = NCOL, stat = 0, Nwork = 2 * NCOL, j;
  Real tr;
  double *store, *work, *Rwork, *eigs, norm = NUMLINK * volume;
  double *ave_eigs, *sq_eigs, width, *min_eigs, *max_eigs;
  complex tc;
  su3_matrix_f USq, tmat;

#ifdef SCALAR_EIG_DIST
  if (this_node != 0) {
    printf("plaquette: don't run SCALAR_EIG_DIST in parallel\n");
    fflush(stdout);
    terminate(1);
  }
#endif

  // Allocate double arrays to be used by LAPACK
  store = malloc(2 * NCOL * NCOL * sizeof(*store));
  work = malloc(2 * Nwork * sizeof(*work));
  Rwork = malloc((3 * NCOL - 2) * sizeof(*Rwork));
  eigs = malloc(NCOL * sizeof(*eigs));
  ave_eigs = malloc(NCOL * sizeof(*ave_eigs));
  max_eigs = malloc(NCOL * sizeof(*max_eigs));
  min_eigs = malloc(NCOL * sizeof(*min_eigs));
  sq_eigs = malloc(NCOL * sizeof(*ave_eigs));
  for (j = 0; j < NCOL; j++) {
    ave_eigs[j] = 0.0;
    min_eigs[j] =  99.0;
    max_eigs[j] = -99.0;
    sq_eigs[j] = 0.0;
  }

  FORALLSITES(i, s) {
    for (dir = XUP; dir < NUMLINK; dir++) {
      if (project == 1)     // Consider polar-projected scalar fields
        polar(&(s->linkf[dir]), &tmat, &USq);
      else                  // Consider U.Ubar scalar fields
        mult_su3_na_f(&(s->linkf[dir]), &(s->linkf[dir]), &USq);

      // In either case, take traceless part
      tc = trace_su3_f(&USq);
      tr = tc.real / (Real)NCOL;
      for (j = 0; j < NCOL; j++)
        USq.e[j][j].real -= tr;
//      dumpmat_f(&USq);

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

  // Print averages, extrema and square root of variance
  // Format: SCALAR_EIG ave1 ... aveN min max width
  for (j = 0; j < NCOL; j++) {
    g_doublesum(&(ave_eigs[j]));
    g_doublesum(&(sq_eigs[j]));
    ave_eigs[j] /= norm;
    sq_eigs[j] /= norm;
  }

  for (j = 0; j < NCOL; j++) {
    width = sqrt(sq_eigs[j] - ave_eigs[j] * ave_eigs[j]);
    if (project == 1) {
      node0_printf("POLAR_EIG ");
    }
    else
      node0_printf("UUBAR_EIG ");
    node0_printf("%d %.6g %.6g %.6g %.6g\n",
                 j, ave_eigs[j], width, min_eigs[j], max_eigs[j]);
  }

  // Free double arrays used by LAPACK
  free(store);
  free(work);
  free(Rwork);
  free(eigs);
  free(ave_eigs);
  free(min_eigs);
  free(max_eigs);
  free(sq_eigs);
}
// -----------------------------------------------------------------
