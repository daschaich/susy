// -----------------------------------------------------------------
// Measure the average value of Tr[U.Ubar] / N
// as well as the width sqrt(<tr^2> - <tr>^2) of its distribution
// At the same time, measure (real) determinants of traceless parts
#include "susy_includes.h"

double link_trace(double *linktr, double *width,
                  double *link_det, double *det_ave, double *det_width) {

  register int i, j, dir;
  register site *s;
  double ave = 0.0, linktrSq = 0.0, td, detSq = 0.0;//, im_check = 0.0;
  complex tc;
  matrix_f tmat;

  *det_ave = 0.0;
  FORALLDIR(dir) {
    linktr[dir] = 0.0;
    link_det[dir] = 0.0;
    FORALLSITES(i, s) {
      mult_na_f(&(s->linkf[dir]), &(s->linkf[dir]), &tmat);
      tc = trace_f(&tmat);
      linktr[dir] += tc.real;
      linktrSq += tc.real * tc.real;

      // Determinant of traceless part
      tc.real *= one_ov_N;
      for (j = 0; j < NCOL; j++)
        tmat.e[j][j].real -= tc.real;
      tc = find_det(&tmat);
      link_det[dir] += tc.real;
      detSq += tc.real * tc.real;
//      im_check += tc.imag;
    }
    linktr[dir] *= one_ov_N / ((double)volume);
    g_doublesum(&(linktr[dir]));
    ave += linktr[dir];

    link_det[dir] /= (double)volume;
    g_doublesum(&(link_det[dir]));
    *det_ave += link_det[dir];
  }
  ave /= (double)NUMLINK;
  linktrSq *= one_ov_N * one_ov_N / ((double)volume * NUMLINK);
  g_doublesum(&linktrSq);
  *width = sqrt(linktrSq - ave * ave);

  *det_ave /= (double)NUMLINK;
  detSq /= ((double)volume * NUMLINK);
  g_doublesum(&detSq);
  td = *det_ave;
  *det_width = sqrt(detSq - td * td);

  // Check that imaginary parts of determinants vanish
//  im_check /= (double)volume;
//  g_doublesum(&im_check);
//  node0_printf("FLINK_IM %.6g\n", im_check);

  return ave;
}
// -----------------------------------------------------------------
