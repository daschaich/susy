// -----------------------------------------------------------------
// Print widths sqrt(<P^2> - <P>^2) of distributions for the
// plaquette and plaquette determinant
// Separately consider real and imaginary parts of latter
// Use tempmat and staple as temporary storage
#include "susy_includes.h"

void widths() {
  register int i;
  register site *s;
  double plaq = 0.0, plaqSq = 0.0, norm = 1.0 / volume;
  double re = 0.0, reSq = 0.0, im = 0.0, imSq = 0.0;
  complex tc;
  msg_tag *mtag = NULL, *mtag2 = NULL;
  matrix tmat;

  // gen_pt[0] is U_b(x+a), gen_pt[1] is U_a(x+b)
  mtag = start_gather_site(F_OFFSET(link[1]), sizeof(matrix),
                           goffset[0], EVENANDODD, gen_pt[0]);
  mtag2 = start_gather_site(F_OFFSET(link[0]), sizeof(matrix),
                            goffset[1], EVENANDODD, gen_pt[1]);

  // tempmat = Udag_b(x+a) Udag_a(x) = [U_a(x) U_b(x+a)]^dag
  wait_gather(mtag);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[0]), (matrix *)(gen_pt[0][i]), &tmat);
    adjoint(&tmat, &(tempmat[i]));
  }
  cleanup_gather(mtag);

  // staple = U_a(x+b) Udag_b(x+a) Udag_a(x)
  // tmat = U_b(x) U_a(x+b) Udag_b(x+a) Udag_a(x) = P_ab
  wait_gather(mtag2);
  FORALLSITES(i, s) {
    mult_nn((matrix *)(gen_pt[1][i]), &(tempmat[i]), &(staple[i]));
    mult_nn(&(s->link[1]), &(staple[i]), &tmat);

    tc = trace(&tmat);
    plaq += tc.real;
    plaqSq += tc.real * tc.real;

    tc = find_det(&tmat);
    re += tc.real;
    reSq += tc.real * tc.real;
    im += tc.imag;
    imSq += tc.imag * tc.imag;
  }
  cleanup_gather(mtag2);
  g_doublesum(&plaq);
  g_doublesum(&plaqSq);
  g_doublesum(&re);
  g_doublesum(&reSq);
  g_doublesum(&im);
  g_doublesum(&imSq);

  // Now compute and print square root of variances
  // Format: WIDTHS plaq re im
  plaq *= norm;
  plaqSq *= norm;
  re *= norm;
  reSq *= norm;
  im *= norm;
  imSq *= norm;
  node0_printf("WIDTHS %.6g %.6g %.6g\n", sqrt(plaqSq - plaq * plaq),
               sqrt(reSq - re * re), sqrt(imSq - im * im));
}
// -----------------------------------------------------------------
