// -----------------------------------------------------------------
// Compute Polyakov loop at each spatial site, for correlator
// Use tempmat1, tempmat2 and staple for temporary storage
#include "susy_includes.h"

void ploop_c() {
#ifdef PL_CORR
  register int i;
  register site *s;
  int t;
  su3_matrix_f *mat;

  FORALLSITES(i, s)
    su3mat_copy_f(&(s->linkf[TUP]), &(tempmat1[i]));

  for (t = 1; t < nt; t++) {
    shiftmat(tempmat1, tempmat2, goffset[TUP]);
    FORALLSITES(i, s) {
      su3mat_copy_f(&(tempmat2[i]), &(tempmat1[i]));
      if (s->t != 0)
        continue;
      if (t == 1) {
        mat = (su3_matrix_f *)gen_pt[0][i];
        mult_su3_nn_f(&(s->linkf[TUP]), mat, &(staple[i]));
      }
      else {
        mat = (su3_matrix_f *)gen_pt[0][i];
        mult_su3_nn_f(&(staple[i]), mat, &(tempmat1[i]));
        su3mat_copy_f(&(tempmat1[i]), &(staple[i]));
      }
    }
  }
  FORALLSITES(i, s) {
    if (s->t != 0)
      continue;

    s->ploop_corr = trace_su3_f(&(staple[i]));
    CDIVREAL((s->ploop_corr), NCOL, (s->ploop_corr));
  }
#endif
}
// -----------------------------------------------------------------
