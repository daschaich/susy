// -----------------------------------------------------------------
// Evaluate the Polyakov loop using repeated single-timeslice gathers
// Use tempmat and tempmat2 for temporary storage
#include "susy_includes.h"

complex ploop() {
  register int i, index = node_index(0);
  register site *s;
  int t;
  complex plp = cmplx(0.0, 0.0);
  matrix tmat;

  // Special case: nt == 1 (i.e., a single-site 'lattice', s=0)
  if (nt == 1) {
    plp = trace(&(lattice[index].link));
    return plp;
  }

  // Compute line by steadily shifting links to hyperplane 0
  FORALLSITES(i, s)
    mat_copy(&(s->link), &(tempmat[i]));

  for (t = 1; t < nt; t++) {
    shiftmat(tempmat, tempmat2, TUP);
    if (t == 1)
      mult_nn(&(lattice[index].link), &(tempmat[index]), &tmat);
    else {
      mult_nn(&tmat, &(tempmat[index]), &(tempmat2[index]));
      mat_copy(&(tempmat2[index]), &tmat);
    }
  }
  if (mynode() == node_number(0))
    plp = trace(&tmat);
  g_sync();
  g_complexsum(&plp);   // Broadcast same value to all nodes
  return plp;
}
// -----------------------------------------------------------------
