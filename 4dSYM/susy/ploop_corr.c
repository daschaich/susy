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
    mat_copy_f(&(s->linkf[TUP]), &(tempmat[i]));

  for (t = 1; t < nt; t++) {
    shiftmat(tempmat, tempmat2, goffset[TUP]);
    FORALLSITES(i, s) {
      mat_copy_f(&(tempmat2[i]), &(tempmat[i]));
      if (s->t != 0)
        continue;
      if (t == 1)
        mult_nn_f(&(s->linkf[TUP]), (matrix_f *)gen_pt[0][i], &(staple[i]));
      else {
        mult_nn_f(&(staple[i]), (matrix_f *)gen_pt[0][i], &(tempmat[i]));
        mat_copy_f(&(tempmat[i]), &(staple[i]));
      }
    }
  }
  FORALLSITES(i, s) {
    if (s->t != 0)
      continue;

    s->ploop_corr = trace_f(&(staple[i]));
    CDIVREAL((s->ploop_corr), NCOL, (s->ploop_corr));
  }
#endif
}
// -----------------------------------------------------------------
