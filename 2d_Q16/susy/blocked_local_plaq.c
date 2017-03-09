// -----------------------------------------------------------------
// Measure plaquettes after block RG blocking steps
// Allow sanity check of reproducing local_plaquette() with block <= 0
// Use tempmat, tempmat2 and Fmunu for temporary storage

// #define MIN_PLAQ turns on measurement of minimum plaquette per config
// for tuning smearing as in Hasenfratz & Knechtli, hep-lat/0103029
#define MIN_PLAQ
#include "susy_includes.h"

void blocked_local_plaq(int Nsmear, int block) {
  register int i;
  register site *s;
  int j, stride = 1;
  double sum = 0.0, tr, max_plaq = -200.0 * NCOL;
#ifdef MIN_PLAQ
  double min_plaq = 200.0 * NCOL;
#endif
  matrix tmat, tmat2;

  // Set number of links to stride, bl = 2^block
  for (j = 0; j < block; j++)
    stride *= 2;

  // Compute the bl-strided plaquette
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
    if (tr > max_plaq)
      max_plaq = tr;
#ifdef MIN_PLAQ
    if (tr < min_plaq)
      min_plaq = tr;
#endif
    sum += tr;
  }
  g_doublesum(&sum);

  g_doublemax(&max_plaq);
#ifdef MIN_PLAQ
  // Somewhat hacky since we don't have g_doublemin...
  min_plaq = -min_plaq;
  g_doublemax(&min_plaq);
  min_plaq = -min_plaq;
  node0_printf("BMIN_PLAQ %d %d %.8g", Nsmear, block, min_plaq);
#endif

  // Average over volume
  sum /= ((double)(volume));
  node0_printf(" %.8g %.8g\n", sum, max_plaq);
}
// -----------------------------------------------------------------
