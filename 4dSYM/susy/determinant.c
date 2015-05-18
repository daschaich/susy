// -----------------------------------------------------------------
// Functions for NCOLxNCOL determinants, adjugates and inverses
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// LU decomposition based on Numerical Recipes
void ludcmp_cx(su3_matrix_f *a, int *indx, Real *d) {
  int i, imax, j, k;
  Real big, fdum;
  complex sum, dum, ct;
  Real vv[NCOL];

  *d = 1.0;
  for (i = 0; i < NCOL; i++) {
    big = 0.0;
    for (j = 0; j < NCOL; j++) {
      ct.real = a->e[i][j].real * a->e[i][j].real
              + a->e[i][j].imag * a->e[i][j].imag;
      if (ct.real > big)
        big = ct.real;
    }
    if (big == 0.0) {
      node0_printf("Singular matrix in routine LUDCMP:\n");
      dumpmat_f(a);
      terminate(1);
    }
    vv[i] = 1.0 / sqrt(big);
  }
  for (j = 0; j < NCOL; j++) {
    for (i = 0; i < j; i++) {
      sum = a->e[i][j];
      for (k = 0; k < i; k++) {
        CMUL(a->e[i][k], a->e[k][j], ct);
        CSUB(sum, ct, sum);
      }
      a->e[i][j] = sum;
    }
    big = 0.0;
    imax = 0;
    for (i = j; i < NCOL; i++) {
      sum = a->e[i][j];
      for (k = 0; k < j; k++) {
        CMUL(a->e[i][k], a->e[k][j],ct);
        CSUB(sum, ct, sum);
      }
      a->e[i][j] = sum;

      if ((fdum = vv[i] * fabs(sum.real)) >= big) {
        big = fdum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < NCOL; k++) {
        dum = a->e[imax][k];
        a->e[imax][k] = a->e[j][k];
        a->e[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (j != NCOL - 1) {
      dum = a->e[j][j];
      for (i = j + 1; i < NCOL; i++) {
        CDIV(a->e[i][j], dum, ct);
        a->e[i][j] = ct;
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute complex determinant of given link,
// using Numerical Recipes-based LU decomposition above
complex find_det(su3_matrix_f *Q) {
  complex det;

#if (NCOL == 1)
  det = Q->e[0][0];
#endif
#if (NCOL == 2)
  complex det2;

  CMUL(Q->e[0][0], Q->e[1][1], det);
  CMUL(Q->e[0][1], Q->e[1][0], det2);
  CSUB(det, det2, det);
#endif
#if (NCOL > 2)
  int i, indx[NCOL];
  Real d;
  complex det2;
  su3_matrix_f QQ;

  su3mat_copy_f(Q, &QQ);
  ludcmp_cx(&QQ, indx, &d);
  det = cmplx(d, 0.0);
  for (i = 0; i < NCOL; i++) {
    CMUL(det, QQ.e[i][i], det2);
    det = det2;
  }
#endif
//  printf("DETER %.4g %.4g\n", det.real, det.imag);
  return det;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Cofactor of src matrix omitting given row and column
complex cofactor(su3_matrix_f *src, int row, int col) {
  int a = 0, b = 0, i, j;
  complex submat[NCOL - 1][NCOL - 1], cof = cmplx(0.0, 0.0);

  // Set up submatrix omitting row row and column col
  for (i = 0; i < NCOL; i++) {
    if (i != row) {
      for (j = 0; j < NCOL; j++) {
        if (j != col) {
          set_complex_equal(&(src->e[i][j]), &(submat[a][b]));
          b++;
        }
      }
      a++;
    }
    b = 0;
  }

#if (NCOL == 3)
  // Here things are easy
  complex tc1, tc2;
  CMUL(submat[0][0], submat[1][1], tc1);
  CMUL(submat[1][0], submat[0][1], tc2);
  CSUB(tc1, tc2, cof);
#endif
#if (NCOL == 4)
  // Here things are less fun, but still not tough
  // det = 00(11*22 - 21*12) - 01(10*22 - 20*12) + 02(10*21 - 20*11)
  complex tc1, tc2, tc3, sav;
  CMUL(submat[1][1], submat[2][2], tc1);
  CMUL(submat[2][1], submat[1][2], tc2);
  CSUB(tc1, tc2, tc3);
  CMUL(submat[0][0], tc3, cof);

  CMUL(submat[1][0], submat[2][2], tc1);
  CMUL(submat[2][0], submat[1][2], tc2);
  CSUB(tc2, tc1, tc3);    // Absorb negative sign
  CMUL(submat[0][1], tc3, sav);
  CSUM(cof, sav);

  CMUL(submat[1][0], submat[2][1], tc1);
  CMUL(submat[2][0], submat[1][1], tc2);
  CSUB(tc1, tc2, tc3);
  CMUL(submat[0][2], tc3, sav);
  CSUM(cof, sav);
#endif
#if (NCOL > 4)
  node0_printf("Haven't coded cofactor for more than 4 colors\n");
  terminate(1);
#endif
  return cof;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Transpose of the cofactor matrix
void adjugate(su3_matrix_f *src, su3_matrix_f *dest) {
#if (NCOL == 1)
  dest->e[0][0] = cmplx(1.0, 0.0);
#endif
#if (NCOL == 2)
  dest->e[0][0] = src->e[1][1];
  dest->e[1][1] = src->e[0][0];
  CMULREAL(src->e[0][1], -1.0, dest->e[0][1]);
  CMULREAL(src->e[1][0], -1.0, dest->e[1][0]);
#endif
#if (NCOL == 3 || NCOL == 4)
  // This should work for generic NCOL given corresponding cofactor routine
  int i, j;

  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      dest->e[i][j] = cofactor(src, j, i);    // Note transpose!
      // Extra negative sign for 01, 03, 10, 12, 21, 23, 30 and 32
      if ((i + j) % 2 == 1)
        CNEGATE(dest->e[i][j], dest->e[i][j]);
    }
  }
#endif
#if (NCOL > 4)
  node0_printf("Haven't coded adjugate for more than 4 colors\n");
  exit(1);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute link matrix inverse as adjugate matrix
// normalized by the determinant
void invert(su3_matrix_f *in, su3_matrix_f *out) {
#if (NCOL == 1)
  complex one = cmplx(1.0, 0.0);
  CDIV(one, in->e[0][0], out->e[0][0]);
#endif

#if (NCOL == 2 || NCOL == 3 || NCOL == 4)
  int i, j;
  complex tc, det = find_det(in);
  adjugate(in, out);
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      tc = out->e[i][j];    // Can't re-use out in CDIV
      CDIV(tc, det, out->e[i][j]);
    }
  }
#endif
#if (NCOL > 4)
  node0_printf("Haven't coded invert for more than 4 colors\n");
  exit(1);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Average plaquette determinant over lattice volume
void measure_det() {
  register int i, a, b;
  register site *s;
  Real norm = (Real)(volume * NUMLINK * (NUMLINK - 1) / 2);
  double WSq = 0.0;
  complex tot = cmplx(0.0, 0.0), tot_sq = cmplx(0.0, 0.0), tc;

  compute_plaqdet();
  for (a = YUP; a < NUMLINK; a++) {
    for (b = XUP; b < a; b++) {
      FORALLSITES(i, s) {
        CSUM(tot, plaqdet[a][b][i]);
        tot_sq.real += plaqdet[a][b][i].real * plaqdet[a][b][i].real;
        tot_sq.imag += plaqdet[a][b][i].imag * plaqdet[a][b][i].imag;

        CADD(plaqdet[a][b][i], minus1, tc);
        WSq += cabs_sq(&tc);
      }
    }
  }
  g_complexsum(&tot);
  g_complexsum(&tot_sq);
  g_doublesum(&WSq);
  CDIVREAL(tot, norm, tot);
  CDIVREAL(tot_sq, norm, tot_sq);
  node0_printf("DET %.6g %.6g %.6g %.6g %.6g\n",
               tot.real, tot.imag, tot_sq.real, tot_sq.imag, WSq / norm);
}
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
