// -----------------------------------------------------------------
// Polar projection now uses LAPACK to remove scalar contributions
#include "susy_includes.h"

// Given matrix in, calculate the unitary polar decomposition element
//   in = out.P --> out = in.[1 / P] and P = sqrt[in^dag.in]
// We diagonalize PSq = in^dag.in using LAPACK,
// then project out its inverse square root
// Also return the hermitian matrix P
void polar(su3_matrix_f *in, su3_matrix_f *out, su3_matrix_f *P) {
  char V = 'V', U = 'U';
  int row, col, Npt = NCOL, stat = 0, Nwork = 2 * NCOL;
  double *store, *work, *Rwork, *eigs;
  complex minus1 = cmplx(-1.0, 0.0);
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

  // Check for degenerate eigenvalues (broke previous Jacobi algorithm)
  for (row = 0; row < NCOL; row++) {
    for (col = row + 1; col < NCOL; col++) {
      if (fabs(eigs[row] - eigs[col]) < IMAG_TOL) {
        printf("WARNING: w[%d] = w[%d] = %.8g\n",
               row, col, Pinv.e[row][row].real);
      }
    }
  }

  // Move the results back into su3_matrix_f structures
  // Overwriting PSq with its eigenvectors for projection
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      PSq.e[row][col].real = store[2 * (col * NCOL + row)];
      PSq.e[row][col].imag = store[2 * (col * NCOL + row) + 1];
      P->e[row][col] = cmplx(0.0, 0.0);
      Pinv.e[row][col] = cmplx(0.0, 0.0);
    }
//    node0_printf("%.4g, ", eigs[row]);
    P->e[row][row].real = sqrt(eigs[row]);
    Pinv.e[row][row].real = 1.0 / sqrt(eigs[row]);
  }
  mult_su3_na_f(P, &PSq, &tmat);
  mult_su3_nn_f(&PSq, &tmat, P);
//  node0_printf("\n");
//  dumpmat_f(&PSq);
//  node0_printf("\n");

  // Now project out 1 / sqrt[in^dag.in] to find out
  mult_su3_na_f(&Pinv, &PSq, &tmat);
  mult_su3_nn_f(&PSq, &tmat, &Pinv);
  mult_su3_nn_f(in, &Pinv, out);

  // Check unitarity of out
  mult_su3_na_f(out, out, &PSq);
  c_scalar_add_diag_su3_f(&PSq, &minus1);
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      if (cabs_sq(&(PSq.e[row][col])) > SQ_TOL) {
        printf("Error getting unitary piece: ");
        printf("%.4g > %.4g for [%d, %d]\n",
               cabs(&(PSq.e[row][col])), IMAG_TOL, row, col);

        dumpmat_f(in);
        dumpmat_f(out);
        dumpmat_f(P);
        return;
      }
    }
  }

  // Check hermiticity of P
  su3_adjoint_f(P, &tmat);
  c_scalar_mult_add_su3mat_f(P, &tmat, &minus1, &PSq);
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      if (cabs_sq(&(PSq.e[row][col])) > SQ_TOL) {
        printf("Error getting hermitian piece: ");
        printf("%.4g > %.4g for [%d, %d]\n",
               cabs(&(PSq.e[row][col])), IMAG_TOL, row, col);

        dumpmat_f(in);
        dumpmat_f(out);
        dumpmat_f(P);
        return;
      }
    }
  }

  // Check in = out.P
  mult_su3_nn_f(out, P, &tmat);
  c_scalar_mult_add_su3mat_f(in, &tmat, &minus1, &PSq);
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      if (cabs_sq(&(PSq.e[row][col])) > SQ_TOL) {
        printf("Error reconstructing initial matrix: ");
        printf("%.4g > %.4g for [%d, %d]\n",
               cabs(&(PSq.e[row][col])), IMAG_TOL, row, col);

        dumpmat_f(in);
        dumpmat_f(out);
        dumpmat_f(P);
        return;
      }
    }
  }

  // Free double arrays used by LAPACK
  free(store);
  free(work);
  free(Rwork);
  free(eigs);
}
// -----------------------------------------------------------------
