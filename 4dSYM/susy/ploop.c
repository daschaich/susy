// -----------------------------------------------------------------
// Evaluate the Polyakov loop using repeated single-timeslice gathers
// Despite the file name we don't consider the distribution of ploop values
// at each spatial site
// Use tempmat1, tempmat2 and staple for temporary storage
#include "susy_includes.h"

complex ploop(double *plpMod) {
  register int i;
  register site *s;
  int t;
  double norm = (double)(nx * ny * nz);
  complex sum  = cmplx(0.0, 0.0), plp;

  FORALLSITES(i, s)
    su3mat_copy_f(&(s->linkf[TUP]), &(tempmat1[i]));

  for (t = 1; t < nt; t++) {
    shiftmat(tempmat1, tempmat2, goffset[TUP]);
    FORALLSITES(i, s) {
      if (s->t != 0)
        continue;

      if (t == 1)
        mult_su3_nn_f(&(s->linkf[TUP]), &(tempmat1[i]), &(staple[i]));
      else {
        mult_su3_nn_f(&(staple[i]), &(tempmat1[i]), &(tempmat2[i]));
        su3mat_copy_f(&(tempmat2[i]), &(staple[i]));
      }
    }
  }

  *plpMod = 0.0;
  FORALLSITES(i, s) {
    if (s->t != 0)
      continue;

    plp = trace_su3_f(&(staple[i]));
    CSUM(sum, plp);
    *plpMod += cabs(&plp);
  }
  g_complexsum(&sum);
  g_doublesum(plpMod);
  CDIVREAL(sum, norm, sum);
  *plpMod /= norm;
  return sum;
}
// -----------------------------------------------------------------
