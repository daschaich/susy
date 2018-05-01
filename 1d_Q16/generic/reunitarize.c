// -----------------------------------------------------------------
// Reunitarization for arbitrary NCOL
// 1) Re-unitarize row 0
// 2) For rows 1, ..., NCOL-2, make orthogonal w.r.t previous rows,
//    and orthogonalize
// 3) For the last row, U(NCOL-1, i)=-(epsilon(i, j_0, j_1, j_2..)U(0, j_0)U(1, j_1)...U(NCOL-1, j_{NCOL-1})^*
#include "generic_includes.h"

#define TOLERANCE 0.0001
#define MAXERRCOUNT 100
//#define UNIDEBUG

Real max_deviation;
double av_deviation;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int check_deviation(Real deviation) {
  if (max_deviation < deviation)
    max_deviation = deviation;
  av_deviation += deviation * deviation;

  if (deviation < TOLERANCE)
    return 0;
  else
    return 1;
}
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
int reunit(matrix *c) {
  register Real ar;
  int i, j, k, errors = 0;
  Real deviation;
  complex sum, deter, tc;
  matrix tmat;

  mat_copy(c, &tmat);
  for (i = 0; i < NCOL; i++) {
    // Orthogonalize the ith vector w.r.t. all the lower ones
    sum = cmplx(0.0, 0.0);
    for (j = 0; j < i; j++) {
      for (k = 0; k < NCOL; k++) {
        CMULJ_(c->e[j][k], c->e[i][k], tc);
        CSUM(sum, tc);
      }
      for (k = 0; k < NCOL; k++)
        CMULDIF(sum, c->e[j][k], c->e[i][k]);     // Okay since j < i
    }
    // Normalize ith vector
    ar = (*c).e[i][0].real * (*c).e[i][0].real
       + (*c).e[i][0].imag * (*c).e[i][0].imag;
    for (k = 1; k < NCOL; k++) {
      ar += (*c).e[i][k].real * (*c).e[i][k].real
          + (*c).e[i][k].imag * (*c).e[i][k].imag;
    }
#ifdef UNIDEBUG
    printf("%d %e\n", i, ar);
#endif

    // Usual test
    // For the 0 to NCOL-2 vectors, check that the length didn't shift
    // For the last row, the check is of the determinant
    // Computing the row using the determinant
    deviation = fabs(ar - 1.0);
#ifdef UNIDEBUG
    printf("DEV %e\n", deviation);
#endif
    errors += check_deviation(deviation);

    ar = 1.0 / sqrt((double)ar);             /* used to normalize row */
    for (j = 0; j < NCOL; j++) {
      (*c).e[i][j].real *= ar;
      (*c).e[i][j].imag *= ar;
    }
  }

  // Check the determinant
#ifdef UNIDEBUG
  node0_printf("MATRIX AFTER RE-ORTHONORMAL\n");
  dumpmat(c);
#endif
  deter = find_det(c);
  ar = deter.real * deter.real + deter.imag * deter.imag;
  deviation = fabs(ar - 1.0);
  errors += check_deviation(deviation);

  // Rephase last column
  for (j = 0; j < NCOL; j++) {
    CDIV((*c).e[NCOL - 1][j], deter, tc);
    (*c).e[NCOL - 1][j] = tc;
  }

  // Check new det
  deter = find_det(c);
  ar = deter.real * deter.real + deter.imag * deter.imag;
  deviation = fabs(ar - 1.0);
  errors += check_deviation(deviation);
#ifdef UNIDEBUG
  printf("DET DEV %e\n", deviation);
  node0_printf("NEW\n");
  dumpmat(c);
#endif

  // Print the problematic matrix rather than the updated one
  if (errors)
    dumpmat(&tmat);

  return errors;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void reunitarize() {
  register int i;
  register site *s;
  register matrix *mat;
  int errcount = 0, errors;

  max_deviation = 0.0;
  av_deviation = 0.0;

  FORALLSITES(i, s) {
    mat = &(s->link);
    errors = reunit(mat);
    errcount += errors;
    if (errors) {
      printf("Unitarity problem above (node %d, site %d, tolerance %.4g)\n",
             mynode(), i, TOLERANCE);
    }
    if (errcount > MAXERRCOUNT) {
      printf("Unitarity error count exceeded: %d\n", errcount);
      fflush(stdout);
      terminate(1);
    }
  }

#ifdef UNIDEBUG
  printf("Deviation from unitarity on node %d: max %.4g, ave %.4g\n",
         mynode(), max_deviation, av_deviation);
#endif
  if (max_deviation > TOLERANCE) {
    printf("reunitarize: Node %d unitarity problem, maximum deviation %.4g\n",
           mynode(), max_deviation);
    errcount++;
    if (errcount > MAXERRCOUNT) {
      printf("Unitarity error count exceeded\n");
      terminate(1);
    }
  }
}
// -----------------------------------------------------------------
