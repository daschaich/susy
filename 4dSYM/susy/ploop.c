// -----------------------------------------------------------------
// Evaluate the Polyakov loop using repeated single-timeslice gathers
// Use tempmat1, tempmat2, staple and DmuUmu for temporary storage
#include "susy_includes.h"

// dir tells us the direction of the product: XUP, YUP, ZUP, TUP or DIR_5
// project == 1 tells us to consider links unitarized via polar projection
complex ploop(int dir, int project, double *plpMod) {
  register int i;
  register site *s;
  int t, len = nt;
  double norm = (double)(nx * ny * nz);
  complex sum  = cmplx(0.0, 0.0), plp;
  su3_matrix_f tmat, tmat2;

  // Optionally consider polar-projected links,
  // saving original values in DmuUmu to be reset at end
  if (project == 1) {
    FORALLSITES(i, s) {
      su3mat_copy_f(&(s->linkf[dir]), &(DmuUmu[i]));
      polar(&(s->linkf[dir]), &tmat, &tmat2);
      su3mat_copy_f(&tmat, &(s->linkf[dir]));
    }
  }

  FORALLSITES(i, s)
    su3mat_copy_f(&(s->linkf[dir]), &(tempmat1[i]));

  switch(dir) {
    case XUP:   len = nx; break;
    case YUP:   len = ny; break;
    case ZUP:   len = nz; break;
    case TUP:   len = nt; break;
    case DIR_5: len = nt; break;
    default:
      printf("ploop: unrecognized direction %d\n", dir);
      fflush(stdout);
      terminate(1);
  }

  for (t = 1; t < len; t++) {
    shiftmat(tempmat1, tempmat2, goffset[dir]);
    FORALLSITES(i, s) {
      if (s->t != 0)
        continue;

      if (t == 1)
        mult_su3_nn_f(&(s->linkf[dir]), &(tempmat1[i]), &(staple[i]));
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

  if (project == 1) {       // Reset original links
    FORALLSITES(i, s)
      su3mat_copy_f(&(DmuUmu[i]), &(s->linkf[dir]));
  }

  return sum;
}
// -----------------------------------------------------------------
