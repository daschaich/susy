// -----------------------------------------------------------------
// Create link in adjoint rep from linkf in fundamental rep
#include "susy_includes.h"

#if (DIMF != NCOL * NCOL)
  #error "Fermions must be in adjoint rep!"
#endif

void fermion_rep() {
  register int mu, i, j, k;
  register site *s;
  su3_matrix_f tmat1, tmat2;

  FORALLSITES(i, s) {
    for (mu = 0; mu < NUMLINK; mu++) {
      for (j = 0; j < DIMF; j++) {
        for (k = 0; k < DIMF; k++) {
          mult_su3_nn_f(&(s->linkf[mu]), &Lambda[k], &tmat1);
          mult_su3_nn_f(&Lambda[j], &tmat1, &tmat2);
          s->link[mu].e[j][k] = trace_su3_f(&tmat2);
        }
      }
//      node0_printf("Site %d mu %d\n", i, mu);
//      dumpmat(&(s->link[mu]));
    }
  }
}
// -----------------------------------------------------------------
