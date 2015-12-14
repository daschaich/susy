// -----------------------------------------------------------------
// Measure plaquettes after block RG blocking steps
// Allow sanity check of reproducing local_plaquette() with block <= 0
// Use tempmat1 and tempmat2 for temporary storage

// #define MIN_PLAQ turns on measurement of minimum plaquette per config
// for tuning smearing as in Hasenfratz & Knechtli, hep-lat/0103029
#define MIN_PLAQ
#include "susy_includes.h"

void blocked_local_plaq(int Nsmear, int block) {
  register int i, dir, dir2;
  register site *s;
  int j, stride = 1;
  double ss_sum = 0.0, st_sum = 0.0, tr;
#ifdef MIN_PLAQ
  double min_plaq = 200.0 * NCOL, max_plaq = -200.0 * NCOL;
#endif
  su3_matrix_f tmat, tmat2;

  // Set number of links to stride, bl = 2^block
  for (j = 0; j < block; j++)
    stride *= 2;

  // Compute the bl-strided plaquette, exploiting a symmetry under dir<-->dir2
  for (dir = YUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      // Copy links to tempmat1 and tempmat2 to be shifted
      FORALLSITES(i, s) {
        su3mat_copy_f(&(s->linkf[dir]), &(tempmat1[i]));
        su3mat_copy_f(&(s->linkf[dir2]), &(tempmat2[i]));
      }

      // Get mom[dir2] from dir and mom[dir] from dir2, both with stride
      // This order may be easier on cache
      for (j = 0; j < stride; j++)
        shiftmat(tempmat2, Fmunu[2], goffset[dir]);
      for (j = 0; j < stride; j++)
        shiftmat(tempmat1, Fmunu[1], goffset[dir2]);

      // Compute tmat  = U_1(x) U_2(x + dir)
      //     and tmat2 = U_2(x) U_1(x + dir2)
      // then plaq = realtrace(tmat2, tmat)[ U_1(x) ]
      //           = tr[Udag_1(x + dir2) Udag_2(x) U_1(x) U_2(x + dir)]
      FORALLSITES(i, s) {
        mult_su3_nn_f(&(s->linkf[dir]), &(tempmat2[i]), &tmat);
        mult_su3_nn_f(&(s->linkf[dir2]), &(tempmat1[i]), &tmat2);
        tr = (double)realtrace_su3_f(&tmat2, &tmat);
#ifdef MIN_PLAQ
        if (tr < min_plaq)
          min_plaq = tr;
        if (tr > max_plaq)
          max_plaq = tr;
#endif
        if (dir == TUP || dir2 == TUP)
          st_sum += tr;
        else
          ss_sum += tr;
      }
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);

#ifdef MIN_PLAQ
  // Somewhat hacky since we don't have g_doublemin...
  g_doublemax(&max_plaq);
  min_plaq = -min_plaq;
  g_doublemax(&min_plaq);
  min_plaq = -min_plaq;
  node0_printf("BMIN_PLAQ %d %d %.8g", Nsmear, block, min_plaq);
#endif

  // Average over four plaquettes that involve the temporal link
  // and six that do not
  ss_sum /= ((double)(6.0 * volume));
  st_sum /= ((double)(4.0 * volume));
  node0_printf(" %.8g %.8g %.8g\n", ss_sum, st_sum, max_plaq);
}
// -----------------------------------------------------------------
