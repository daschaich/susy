// -----------------------------------------------------------------
// Print widths sqrt(<P^2> - <P>^2) of distributions for the
// plaquette and plaquette determinant
// Separately consider real and imaginary parts of latter
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use tempmat1 and staple as temporary storage
void widths() {
  register int i;
  register site *s;
  int a, b;
  double plaq = 0.0, plaqSq = 0.0, norm = 10.0 * volume;
  double re = 0.0, reSq = 0.0, im = 0.0, imSq = 0.0;
  complex tc;
  msg_tag *mtag0 = NULL, *mtag1 = NULL;
  su3_matrix_f tmat, *mat;

  for (a = XUP; a < NUMLINK; a++) {
    for (b = a + 1; b < NUMLINK; b++) {
      // gen_pt[0] is U_b(x+a), gen_pt[1] is U_a(x+b)
      mtag0 = start_gather_site(F_OFFSET(linkf[b]), sizeof(su3_matrix_f),
                                goffset[a], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(linkf[a]), sizeof(su3_matrix_f),
                                goffset[b], EVENANDODD, gen_pt[1]);

      // tempmat1 = Udag_b(x+a) Udag_a(x) = [U_a(x) U_b(x+a)]^dag
      wait_gather(mtag0);
      FORALLSITES(i, s) {
        mat = (su3_matrix_f *)(gen_pt[0][i]);
        mult_su3_nn_f(&(s->linkf[a]), mat, &tmat);
        su3_adjoint_f(&tmat, &(tempmat1[i]));
      }
      cleanup_gather(mtag0);

      // staple = U_a(x+b) Udag_b(x+a) Udag_a(x)
      // tmat = U_b(x) U_a(x+b) Udag_b(x+a) Udag_a(x) = P_ab
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mat = (su3_matrix_f *)(gen_pt[1][i]);
        mult_su3_nn_f(mat, &(tempmat1[i]), &(staple[i]));
        mult_su3_nn_f(&(s->linkf[b]), &(staple[i]), &tmat);

        tc = trace_su3_f(&tmat);
        plaq += tc.real;
        plaqSq += tc.real * tc.real;

        tc = find_det(&tmat);
        re += tc.real;
        reSq += tc.real * tc.real;
        im += tc.imag;
        imSq += tc.imag * tc.imag;
      }
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&plaq);
  g_doublesum(&plaqSq);
  g_doublesum(&re);
  g_doublesum(&reSq);
  g_doublesum(&im);
  g_doublesum(&imSq);

  // Now compute and print square root of variances
  // Format: WIDTHS plaq re im
  plaq /= norm;
  plaqSq /= norm;
  re /= norm;
  reSq /= norm;
  im /= norm;
  imSq /= norm;
  node0_printf("WIDTHS %.6g %.6g %.6g\n", sqrt(plaqSq - plaq * plaq),
               sqrt(reSq - re * re), sqrt(imSq - im * im));
}
// -----------------------------------------------------------------
