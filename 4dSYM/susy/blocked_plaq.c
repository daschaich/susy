// -----------------------------------------------------------------
// Evaluate the plaquette after block RG blocking steps
// Also print blocked det and widths sqrt(<P^2> - <P>^2) of distributions
// Use general_gathers; lattice must be divisible by 2^block in all dirs
// Use tempmat1 and tempmat2 for temporary storage
#include "susy_includes.h"

void blocked_plaq(int Nsmear, int block) {
  register int i, dir, dir2;
  register site *s;
  register su3_matrix_f *m1, *m4;
  int j, bl = 2, d1[4] = {0, 0, 0, 0}, d2[4] = {0, 0, 0, 0};
  double plaq = 0.0, plaqSq = 0.0, re = 0.0, reSq = 0.0, im = 0.0, imSq = 0.0;
  double ss_sum = 0.0, st_sum = 0.0, norm = 10.0 * volume, tr;
  complex det = cmplx(0.0, 0.0), tc;
  msg_tag *mtag;
  su3_matrix_f *mat, tmat, tmat2;

  // Set number of links to stride, bl = 2^block
  // Allow sanity check of reproducing d_plaquette() with this routine
  for (j = 1; j < block; j++)
    bl *= 2;
  if (block <= 0)
    bl = 1;

  // Compute the bl-strided plaquette, exploiting a symmetry under dir<-->dir2
  for (dir = YUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      for (j = 0; j < NDIMS; j++) {
        d1[j] = bl * offset[dir][j];
        d2[j] = bl * offset[dir2][j];
      }
      // Can only have one general gather at once...
      mtag = start_general_gather_site(F_OFFSET(linkf[dir2]),
                                       sizeof(su3_matrix_f), d1,
                                       EVENANDODD, gen_pt[0]);

      // tempmat2 = Udag_b(x) U_a(x)
      FORALLSITES(i, s) {
        m1 = &(s->linkf[dir]);
        m4 = &(s->linkf[dir2]);
        mult_su3_an_f(m4, m1, &(tempmat2[i]));
      }

      // Copy first gather to tempmat1
      wait_general_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(tempmat1[i]));
      cleanup_general_gather(mtag);

      mtag = start_general_gather_site(F_OFFSET(linkf[dir]),
                                       sizeof(su3_matrix_f), d2,
                                       EVENANDODD, gen_pt[0]);
      wait_general_gather(mtag);

      // Compute tr[Udag_a(x+d2) Udag_b(x) U_a(x) U_b(x+d1)]
      // Also monitor determinant
      FORALLSITES(i, s) {
        mat = (su3_matrix_f *)(gen_pt[0][i]);
        mult_su3_nn_f(&(tempmat2[i]), &(tempmat1[i]), &tmat);
        mult_su3_an_f(mat, &tmat, &tmat2);
        tc = trace_su3_f(&tmat2);
        tr = (double)tc.real;
        plaq += tr;
        plaqSq += tr * tr;
        if (dir == TUP || dir2 == TUP)
          st_sum += tr;
        else
          ss_sum += tr;

        su3_adjoint_f(&tmat2, &tmat);   // Match sign conventions
        tc = find_det(&tmat);
        CSUM(det, tc);
        re += tc.real;
        reSq += tc.real * tc.real;
        im += tc.imag;
        imSq += tc.imag * tc.imag;
      }
      cleanup_general_gather(mtag);
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
