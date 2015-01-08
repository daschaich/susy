// -----------------------------------------------------------------
// Divide out determinant to project links from gl(N, C) to sl(N, C)
// Project to polar representation to remove scalar contributions
// Polar projection now uses LAPACK
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Divide the determinant out of the matrix in
void det_project(su3_matrix_f *in, su3_matrix_f *out) {
  Real frac = -1.0 / (Real)NCOL;
  complex tc1, tc2;

  tc1 = find_det(in);
  tc2 = clog(&tc1);
  CMULREAL(tc2, frac, tc1);
  tc2 = cexp(&tc1);
  c_scalar_mult_su3mat_f(in, &tc2, out);

#ifdef DEBUG_CHECK
  // Sanity check
  tc1 = find_det(out);
  if (fabs(tc1.imag) > IMAG_TOL || fabs(1.0 - tc1.real) > IMAG_TOL) {
    printf("node%d WARNING: det = (%.4g, %.4g) after projection...\n",
        this_node, tc1.real, tc1.imag);
  }
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Given matrix in, calculate the unitary polar decomposition element
//   in = out.P --> out = in.[1 / P] and P = sqrt[in^dag.in]
// We diagonalize PSq = in^dag.in using LAPACK,
// then project out its inverse square root
// Not returning the real matrix, but could add that if necessary
void polar(su3_matrix_f *in, su3_matrix_f *out) {
  char V = 'V', U = 'U';
  int row, col, Npt = NCOL, stat = 0, Nwork = 2 * NCOL;
  double *store, *work, *Rwork, *eigs, testtrace;
  su3_matrix_f PSq, Pinv, tmat;

  // Allocate double arrays to be used by LAPACK
  store = malloc(2 * NCOL * NCOL * sizeof(*store));
  work = malloc(2 * Nwork * sizeof(*work));
  Rwork = malloc((3 * NCOL - 2) * sizeof(*Rwork));
  eigs = malloc(NCOL * sizeof(*eigs));

  // Convert PSq to column-major double array used by LAPACK
  mult_su3_an_f(in, in, &PSq);
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      store[2 * (col * NCOL + row)] = PSq.e[row][col].real;
      store[2 * (col * NCOL + row) + 1] = PSq.e[row][col].imag;
    }
  }
//  dumpmat_f(&PSq);

  // Compute eigenvalues and eigenvectors
  zheev_(&V, &U, &Npt, store, &Npt, eigs, work, &Nwork, Rwork, &stat);

  // Move the results back into su3_matrix_f structures
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      PSq.e[row][col].real = store[2 * (col * NCOL + row)];
      PSq.e[row][col].imag = store[2 * (col * NCOL + row) + 1];
      Pinv.e[row][col] = cmplx(0.0, 0.0);
    }
//    node0_printf("%.4g, ", eigs[row]);
    Pinv.e[row][row].real = 1.0 / sqrt(eigs[row]);
  }
//  node0_printf("\n");
//  dumpmat_f(&PSq);
//  node0_printf("\n");

  // Check for degenerate eigenvalues in Pinv
  for (row = 0; row < NCOL; row++) {
    for (col = row + 1; col < NCOL; col++) {
      if (fabs(Pinv.e[row][row].real - Pinv.e[col][col].real) < 1.0e-6) {
        printf("WARNING: w[%d] = w[%d] = %.8g\n",
               row, col, Pinv.e[row][row].real);
      }
    }
  }

  // Now project out 1 / sqrt[in^dag.in] to find out
  mult_su3_na_f(&Pinv, &PSq, &tmat);
  mult_su3_nn_f(&PSq, &tmat, &Pinv);
  mult_su3_nn_f(in, &Pinv, out);

  // Check unitarity of out
  testtrace = realtrace_su3_f(out, out);
  if (fabs(testtrace / (Real)NCOL - 1.0) > 1.0e-6)
    printf("Error getting unitary piece: trace = %.4g\n", testtrace);

  // Free double arrays used by LAPACK
  free(store);
  free(work);
  free(Rwork);
  free(eigs);
}
// -----------------------------------------------------------------
