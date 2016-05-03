// -----------------------------------------------------------------
// Evaluate the plaquette after block RG blocking steps
// Also print blocked det and widths sqrt(<P^2> - <P>^2) of distributions
// Allow sanity check of reproducing plaquette() with block <= 0
// Use tempmat, tempmat2 and Fmunu for temporary storage
#include "susy_includes.h"

void blocked_plaq(int Nsmear, int block) {
  register int i, dir, dir2;
  register site *s;
  int j, stride = 1;
  double plaq = 0.0, plaqSq = 0.0, re = 0.0, reSq = 0.0, im = 0.0, imSq = 0.0;
  double ss_sum = 0.0, st_sum = 0.0, norm = 10.0 * volume, tr;
  complex det = cmplx(0.0, 0.0), tc;
  matrix_f tmat, tmat2, tmat3;

  // Set number of links to stride, 2^block
  for (j = 0; j < block; j++)
    stride *= 2;

  // Compute the strided plaquette, exploiting a symmetry under dir<-->dir2
  for (dir = YUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      // Copy links to tempmat and tempmat2 to be shifted
      FORALLSITES(i, s) {
        mat_copy_f(&(s->linkf[dir]), &(tempmat[i]));
        mat_copy_f(&(s->linkf[dir2]), &(tempmat2[i]));
      }

      // Get mom[dir2] from dir and mom[dir] from dir2, both with stride
      // This order may be easier on cache
      for (j = 0; j < stride; j++)
        shiftmat(tempmat2, Fmunu[2], goffset[dir]);
      for (j = 0; j < stride; j++)
        shiftmat(tempmat, Fmunu[1], goffset[dir2]);

      // Compute tmat  = U_1(x) U_2(x + dir)
      //     and tmat2 = U_2(x) U_1(x + dir2)
      // then plaq = realtrace(tmat2, tmat)[ U_1(x) ]
      //           = tr[Udag_1(x + dir2) Udag_2(x) U_1(x) U_2(x + dir)]
      FORALLSITES(i, s) {
        mult_nn_f(&(s->linkf[dir]), &(tempmat2[i]), &tmat);
        mult_nn_f(&(s->linkf[dir2]), &(tempmat[i]), &tmat2);
        tr = (double)realtrace_f(&tmat2, &tmat);
        plaq += tr;
        plaqSq += tr * tr;
        if (dir == TUP || dir2 == TUP)
          st_sum += tr;
        else
          ss_sum += tr;

        // Also monitor determinant
        // (na instead of an to match sign conventions)
        mult_na_f(&tmat2, &tmat, &tmat3);
        tc = find_det(&tmat3);
        CSUM(det, tc);
        re += tc.real;
        reSq += tc.real * tc.real;
        im += tc.imag;
        imSq += tc.imag * tc.imag;
      }
    }
  }
  g_doublesum(&plaq);
  g_doublesum(&plaqSq);
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);
  g_complexsum(&det);
  g_doublesum(&re);
  g_doublesum(&reSq);
  g_doublesum(&im);
  g_doublesum(&imSq);

  // Average over four plaquettes that involve the temporal link
  // and six that do not
  ss_sum /= ((double)(6.0 * volume));
  st_sum /= ((double)(4.0 * volume));
  tr = (ss_sum + st_sum) / 2.0;
  node0_printf("BPLAQ %d %d %.8g %.8g %.8g\n",
               Nsmear, block, ss_sum, st_sum, tr);

  CDIVREAL(det, norm, det);
  node0_printf("BDET %d %d %.6g %.6g\n", Nsmear, block, det.real, det.imag);

  plaq /= norm;
  plaqSq /= norm;
  re /= norm;
  reSq /= norm;
  im /= norm;
  imSq /= norm;
  node0_printf("BWIDTHS %d %d %.6g %.6g %.6g\n", Nsmear, block,
               sqrt(plaqSq - plaq * plaq),
               sqrt(reSq - re * re), sqrt(imSq - im * im));
}
// -----------------------------------------------------------------
