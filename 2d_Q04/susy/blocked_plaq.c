// -----------------------------------------------------------------
// Evaluate the plaquette after block RG blocking steps
// Also print blocked det and widths sqrt(<P^2> - <P>^2) of distributions
// Allow sanity check of reproducing plaquette() with block <= 0
// Use tempmat, tempmat2 and Fmunu for temporary storage
#include "susy_includes.h"

void blocked_plaq(int Nsmear, int block) {
  register int i;
  register site *s;
  int j, stride = 1;
  double plaq = 0.0, plaqSq = 0.0, re = 0.0, reSq = 0.0, im = 0.0, imSq = 0.0;
  double sum = 0.0, norm = 1.0 / volume, tr;
  complex det = cmplx(0.0, 0.0), tc;
  matrix tmat, tmat2, tmat3;

  // Set number of links to stride, 2^block
  for (j = 0; j < block; j++)
    stride *= 2;

  // Compute the strided plaquette, exploiting a symmetry under TUP<-->XUP
  // Copy links to tempmat and tempmat2 to be shifted
  FORALLSITES(i, s) {
    mat_copy(&(s->link[TUP]), &(tempmat[i]));
    mat_copy(&(s->link[XUP]), &(tempmat2[i]));
  }

  // Get mom[XUP] from TUP and mom[TUP] from XUP, both with stride
  // This order may be easier on cache
  for (j = 0; j < stride; j++)
    shiftmat(tempmat2, Fmunu, goffset[TUP]);
  for (j = 0; j < stride; j++)
    shiftmat(tempmat, Fmunu, goffset[XUP]);

  // Compute tmat  = U_1(x) U_2(x + TUP)
  //     and tmat2 = U_2(x) U_1(x + XUP)
  // then plaq = realtrace(tmat2, tmat)[ U_1(x) ]
  //           = tr[Udag_1(x + XUP) Udag_2(x) U_1(x) U_2(x + TUP)]
  FORALLSITES(i, s) {
    mult_nn(&(s->link[TUP]), &(tempmat2[i]), &tmat);
    mult_nn(&(s->link[XUP]), &(tempmat[i]), &tmat2);
    tr = (double)realtrace(&tmat2, &tmat);
    plaq += tr;
    plaqSq += tr * tr;
    sum += tr;

    // Also monitor determinant
    // (na instead of an to match sign conventions)
    mult_na(&tmat2, &tmat, &tmat3);
    tc = find_det(&tmat3);
    CSUM(det, tc);
    re += tc.real;
    reSq += tc.real * tc.real;
    im += tc.imag;
    imSq += tc.imag * tc.imag;
  }
  g_doublesum(&plaq);
  g_doublesum(&plaqSq);
  g_doublesum(&sum);
  g_complexsum(&det);
  g_doublesum(&re);
  g_doublesum(&reSq);
  g_doublesum(&im);
  g_doublesum(&imSq);

  // Average over volume
  node0_printf("BPLAQ %d %d %.8g\n", Nsmear, block, sum * norm);

  CMULREAL(det, norm, det);
  node0_printf("BDET %d %d %.6g %.6g\n", Nsmear, block, det.real, det.imag);

  plaq *= norm;
  plaqSq *= norm;
  re *= norm;
  reSq *= norm;
  im *= norm;
  imSq *= norm;
  node0_printf("BWIDTHS %d %d %.6g %.6g %.6g\n", Nsmear, block,
               sqrt(plaqSq - plaq * plaq),
               sqrt(reSq - re * re), sqrt(imSq - im * im));
}
// -----------------------------------------------------------------
