// -----------------------------------------------------------------
// Average eigenvalues of traceless part of U.Ubar
// Use LAPACK for arbitrary NCOL
// TODO: Add widths and distributions, maybe minimum and maximum

// #define SCALAR_EIG_DIST prints out all eigenvalues for plotting distribution
// CAUTION: Do not run PLAQ_DIST with MPI!
#define SCALAR_EIG_DIST
#include "susy_includes.h"

void scalar_eig() {
  register int i, dir;
  register site *s;
  char N = 'N';     // Ask LAPACK only for eigenvalues
  char U = 'U';     // Have LAPACK store upper triangle of U.Ubar
  int row, col, Npt = NCOL, stat = 0, Nwork = 2 * NCOL, j;
  Real tr;
  double *store, *work, *Rwork, *eigs;
  complex tc;
  su3_matrix_f USq;

  // Allocate double arrays to be used by LAPACK
  store = malloc(2 * NCOL * NCOL * sizeof(*store));
  work = malloc(2 * Nwork * sizeof(*work));
  Rwork = malloc((3 * NCOL - 2) * sizeof(*Rwork));
  eigs = malloc(NCOL * sizeof(*eigs));

  FORALLSITES(i, s) {
    for (dir = XUP; dir < NUMLINK; dir++) {
      // Construct traceless U.Ubar
      mult_su3_na_f(&(s->linkf[dir]), &(s->linkf[dir]), &USq);
      tc = trace_su3_f(&USq);
      tr = tc.real / (Real)NCOL;
      for (j = 0; j < NCOL; j++)
        USq[j].e[k][k] -= tr;
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
        printf("WARNING: Non-zero return for (%d %d %d %d)[%d]\n",
               s->x, s->y, s->z, s->t, dir);
      }

      // TODO: Average eigenvalues, monitor minimum and maximum
    }
  }

  // Free double arrays used by LAPACK
  free(store);
  free(work);
  free(Rwork);
  free(eigs);
}
// -----------------------------------------------------------------
