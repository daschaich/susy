// -----------------------------------------------------------------
// Reunitarization for arbitrary NCOL
// 1) Re-innitarize row 0
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
void reunit_report_problem_matrix(matrix *mat, int i, int j) {
  int ii, jj;
  union {
    Real fval;
    int ival;
  } ifval;

  printf("Unitarity problem on node %d, site %d, field %d tolerance %.4g\n",
          mynode(), i, j, TOLERANCE);
  printf("SU(N) matrix:\n");
  for (ii = 0; ii < NCOL; ii++) {
    for (jj = 0; jj < NCOL; jj++) {
      printf("%f ",(*mat).e[ii][jj].real);
      printf("%f ",(*mat).e[ii][jj].imag);
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
  fflush(stdout);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void ludcmp_cx(matrix *a, int *indx, Real *d) {
  int i, imax, j, k;
  Real big, fdum;
  complex sum, dum, tc;
  Real vv[NCOL];

  *d = 1.0;
  for (i = 0; i < NCOL; i++) {
    big = 0.0;
    for (j = 0; j < NCOL; j++) {
      tc.real = (*a).e[i][j].real*(*a).e[i][j].real
              + (*a).e[i][j].imag*(*a).e[i][j].imag;
      if (tc.real > big)
        big = tc.real;
    }
    if (big == 0.0) {
      node0_printf("ludcmp_cx: Singular matrix\n");
      exit(1);
    }
    vv[i] = 1.0 / sqrt(big);
  }
  // TODO: May now be more efficient complex macros...
  for (j = 0; j < NCOL; j++) {
    for (i = 0; i < j; i++) {
      sum= (*a).e[i][j];
      for (k = 0; k < i; k++) {
        CMUL((*a).e[i][k],(*a).e[k][j], tc);
        CSUB(sum, tc, sum);
      }
      (*a).e[i][j]=sum;
    }
    big = 0.0;
    for (i = j; i < NCOL; i++) {
      sum = (*a).e[i][j];
      for (k = 0; k < j; k++) {
        CMUL((*a).e[i][k],(*a).e[k][j], tc);
        CSUB(sum, tc, sum);
      }
      (*a).e[i][j] = sum;

      if ((fdum = vv[i] * fabs(sum.real)) >= big) {
        big = fdum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < NCOL; k++) {
        dum = (*a).e[imax][k];
        (*a).e[imax][k] = (*a).e[j][k];
        (*a).e[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (j != NCOL-1) {
      dum = (*a).e[j][j];
      for (i = j + 1; i < NCOL; i++) {
        CDIV((*a).e[i][j], dum, tc);
        (*a).e[i][j] = tc;
      }
    }
  }
}

complex find_det(matrix *Q) {
  complex det, det2;
  int i, indx[NCOL];
  Real d;
  matrix QQ;

  mat_copy(Q, &QQ);
  ludcmp_cx(&QQ, indx, &d);
  det = cmplx(d, 0.0);
  for (i = 0; i < NCOL; i++) {
    CMUL(det, QQ.e[i][i], det2);
    det = det2;
  }
#ifdef UNIDEBUG
  printf("DET %e %e\n", det.real, det.imag);
#endif

  return det;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int reunit(matrix *c) {
  Real deviation;
  int i, j, k, errors = 0;
  register Real ar;
  complex sum, deter, tc, tc2;

#ifdef UNIDEBUG
  matrix cold;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      cold.e[i][j] = (*c).e[i][j];
  }
  node0_printf("OLD MATRIX\n");
  dump_mat(&cold);
#endif

  for (i = 0; i < NCOL; i++) {
    // Orthogonalize the ith vector w.r.t. all the lower ones
    // TODO: May now be more efficient complex macros...
    sum = cmplx(0.0, 0.0);
    for (j = 0; j < i; j++) {
      for (k = 0; k < NCOL; k++) {
        CMULJ_(c->e[j][k], c->e[i][k], tc);
        CSUM(sum, tc);
      }
      for (k = 0; k < NCOL; k++) {
        CMUL(sum, c->e[j][k], tc2);
        CSUB(c->e[i][k], tc2, c->e[i][k]);
      }
    }
    // Normalize ith vector
    ar = 0;
    for (k = 0; k < NCOL; k++) {
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
  dump_mat(c);
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
  dump_mat(c);
#endif

  return errors;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void reunitarize() {
  register int i, j;
  register site *s;
  register matrix *mat;
  int errcount = 0, errors;

  max_deviation = 0.0;
  av_deviation = 0.0;

  FORALLSITES(i, s) {
    mat = (matrix *)&(s->link);
    errors = reunit(mat);
    errcount += errors; // Gauge field printed as scalar -1
    if (errors)
      reunit_report_problem_matrix(mat, i, -1);
    if (errcount > MAXERRCOUNT) {
      printf("Unitarity error count exceeded: %d\n", errcount);
      fflush(stdout);
      terminate(1);
    }
    for (j = 0; j < NSCALAR; j++) {
      mat = (matrix *)&(s->X[j]);
      errors = reunit(mat);
      errcount += errors;
      if (errors)
        reunit_report_problem_matrix(mat, i, j);
      if (errcount > MAXERRCOUNT) {
        printf("Unitarity error count exceeded: %d\n", errcount);
        fflush(stdout);
        terminate(1);
      }
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
