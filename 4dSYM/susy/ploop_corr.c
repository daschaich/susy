// -----------------------------------------------------------------
// Based on ploop2.c
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void ploop_c() {
#ifdef PL_CORR
  register int i, t;
  register site *s;
  msg_tag *tag = NULL;

  FORALLSITES(i, s)
    lattice[i].tempmat2 = lattice[i].linkf[TUP];

  for (t = 1; t < nt; t++) {
    tag = start_gather_site(F_OFFSET(tempmat2), sizeof(su3_matrix_f),
                            TUP, EVENANDODD, gen_pt[0] );
    wait_gather(tag);
    FORALLSITES(i, s) {
      if (s->t != 0)
        continue;
      if (t == 1)
          mult_su3_nn_f(&(s->linkf[TUP]), (su3_matrix_f *)gen_pt[0][i],
                        &(s->staple));
      else {
        mult_su3_nn_f(&(s->staple), (su3_matrix_f *)gen_pt[0][i],
                      &(s->tempmat2));
        lattice[i].staple = lattice[i].tempmat2;
      }
    }

    // Need both of these here -- very strange
    FORALLSITES(i, s)
      s->tempmat1 = *(su3_matrix_f *)(gen_pt[0][i]);
    FORALLSITES(i, s)
      s->tempmat2 = s->tempmat1;
    cleanup_gather(tag);
  }
  FORALLSITES(i, s) {
    if (s->t != 0)
      continue;

    s->ploop_corr = trace_su3_f(&(s->staple));
    CDIVREAL((s->ploop_corr), NCOL, (s->ploop_corr));
  }
#endif
}
// -----------------------------------------------------------------
