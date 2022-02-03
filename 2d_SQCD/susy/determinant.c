// -----------------------------------------------------------------
// Functions for NCOLxNCOL determinant and inverse
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute complex determinant of given link through LAPACK
complex find_det(matrix *Q) {
  int i, row, col, Npt = NCOL, stat = 0;
  complex det, det2;
  matrix tmat;

  // Convert Q to column-major double array used by LAPACK
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      i = 2 * (col * NCOL + row);
      store[i] = Q->e[row][col].real;
      store[i + 1] = Q->e[row][col].imag;
    }
  }

  // Compute LU decomposition of in
  zgetrf_(&Npt, &Npt, store, &Npt, ipiv, &stat);

  // Move the results into the matrix structure tmat
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      i = 2 * (col * NCOL + row);
      tmat.e[row][col].real = store[i];
      tmat.e[row][col].imag = store[i + 1];
    }
  }

  det = cmplx(1.0, 0.0);
  for (i = 0; i < NCOL; i++) {
    CMUL(det, tmat.e[i][i], det2);
    // Negate if row has been pivoted according to pivot array
    // Apparently LAPACK sets up ipiv so this doesn't double-count pivots
    // Note 1<=ipiv[i]<=N rather than 0<=ipiv[i]<N because Fortran
    if (ipiv[i] != i + 1) {     // Braces suppress compiler error
      CNEGATE(det2, det);
    }
    else
      det = det2;
  }
//  printf("DETER %.4g %.4g\n", det.real, det.imag);
  return det;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute link matrix inverse through LAPACK
void invert(matrix *in, matrix *out) {
  // Use LAPACK for more than 4 colors
  // Checked that this produces the correct results for NCOL <= 4
  int i, row, col, Npt = NCOL, stat = 0, Nwork = 2 * NCOL;

  // Convert in to column-major double array used by LAPACK
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      i = 2 * (col * NCOL + row);
      store[i] = in->e[row][col].real;
      store[i + 1] = in->e[row][col].imag;
    }
  }

  // Compute LU decomposition of in
  zgetrf_(&Npt, &Npt, store, &Npt, ipiv, &stat);

  // Invert in given its LU decomposition
  zgetri_(&Npt, store, &Npt, ipiv, work, &Nwork, &stat);

  // Move the results into the matrix structure for out
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      i = 2 * (col * NCOL + row);
      out->e[row][col].real = store[i];
      out->e[row][col].imag = store[i + 1];
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Average plaquette determinant over lattice volume
// Assume compute_plaqdet() has already been run
void measure_det() {
  register int i;
  register site *s;
  Real norm = 1.0 / ((Real)(volume));
  double WSq = 0.0;
  complex tot = cmplx(0.0, 0.0), tot_sq = cmplx(0.0, 0.0), tc;

  FORALLSITES(i, s) {
    CSUM(tot, plaqdet[1][0][i]);
    tot_sq.real += plaqdet[1][0][i].real * plaqdet[1][0][i].real;
    tot_sq.imag += plaqdet[1][0][i].imag * plaqdet[1][0][i].imag;

    CADD(plaqdet[1][0][i], minus1, tc);
    WSq += cabs_sq(&tc);
  }
  g_complexsum(&tot);
  g_complexsum(&tot_sq);
  g_doublesum(&WSq);
  CMULREAL(tot, norm, tot);
  CMULREAL(tot_sq, norm, tot_sq);
  node0_printf("DET %.6g %.6g %.6g %.6g %.6g\n",
               tot.real, tot.imag, tot_sq.real, tot_sq.imag, WSq * norm);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Divide the determinant out of the matrix in
void det_project(matrix *in, matrix *out) {
  complex tc, tc2;

  tc = find_det(in);
  tc2 = clog(&tc);
  CMULREAL(tc2, -one_ov_N, tc);
  tc2 = cexp(&tc);
  c_scalar_mult_mat(in, &tc2, out);

#ifdef DEBUG_CHECK
  // Sanity check
  tc = find_det(out);
  if (fabs(tc.imag) > IMAG_TOL || fabs(1.0 - tc.real) > IMAG_TOL) {
    printf("node%d WARNING: det = (%.4g, %.4g) after projection...\n",
           this_node, tc.real, tc.imag);
  }
#endif
}
// -----------------------------------------------------------------
