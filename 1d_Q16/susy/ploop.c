// -----------------------------------------------------------------
// Evaluate the Polyakov loop using repeated single-timeslice gathers
// Use tempmat and tempmat2 for temporary storage
#include "susy_includes.h"

complex ploop() {
  register int i;
  register site *s;
  int t;
  double norm = 1.0/(double) nt;
  complex plp;
  matrix tmat;

  // Special case: nt == 1
  if (nt == 1) {
    FORALLSITES(i, s)
      plp = trace(&(s->link));

    g_complexsum(&plp);
    CMULREAL(plp, norm, plp);

    return plp;
  }

  // Compute line by steadily shifting links to hyperplane 0
  FORALLSITES(i, s)
    mat_copy(&(s->link), &(tempmat[i]));

  for (t = 1; t < nt; t++) {
    shiftmat(tempmat, tempmat2, TUP);
    FORALLSITES(i, s) {
      if (s->t != 0)
        continue;

      if (t == 1)
        mult_nn(&(s->link), &(tempmat[i]), &(tmat));
      else {
        mult_nn(&(tmat), &(tempmat[i]), &(tempmat2[i]));
        mat_copy(&(tempmat2[i]), &(tmat));
      }
    }
  }

  FORALLSITES(i, s) {
    if (s->t != 0)
      continue;

    plp = trace(&(tmat));
  }
  
  g_complexsum(&plp);
  CMULREAL(plp, norm, plp);

  return plp;
}
// -----------------------------------------------------------------
