// -----------------------------------------------------------------
// Compute Polyakov loop at each spatial site, for correlator
// Use tempmat, tempmat2 and staple for temporary storage
#include "susy_includes.h"

void ploop_c() {
#ifdef PL_CORR
  register int i;
  register site *s;
  int t;

  FORALLSITES(i, s)
    mat_copy(&(s->link[TUP]), &(tempmat[i]));

  for (t = 1; t < nt; t++) {
    shiftmat(tempmat, tempmat2, goffset[TUP]);
    FORALLSITES(i, s) {
      mat_copy(&(tempmat2[i]), &(tempmat[i]));
      if (s->t != 0)
        continue;
      if (t == 1)
        mult_nn(&(s->link[TUP]), (matrix *)gen_pt[0][i], &(staple[i]));
      else {
        mult_nn(&(staple[i]), (matrix *)gen_pt[0][i], &(tempmat[i]));
        mat_copy(&(tempmat[i]), &(staple[i]));
      }
    }
  }
  FORALLSITES(i, s) {
    if (s->t != 0)
      continue;

    s->ploop_corr = trace(&(staple[i]));
    CMULREAL((s->ploop_corr), one_ov_N, (s->ploop_corr));
  }
#endif
}
// -----------------------------------------------------------------
