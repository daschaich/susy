// -----------------------------------------------------------------
// Evaluate the Polyakov loop using repeated single-timeslice gathers
// Use tempmat, tempmat2, staple and Uinv[0] for temporary storage
#include "susy_includes.h"

// dir tells us the direction of the product: XUP, YUP, ZUP, TUP or DIR_5
// project == 1 tells us to consider links unitarized via polar projection
complex ploop(int dir, int project, double *plpMod) {
  register int i;
  register site *s;
  int j, t, len = nt;
  double norm = 0.0;
  complex sum  = cmplx(0.0, 0.0), plp;
  matrix_f tmat, tmat2;

  // Optionally consider polar-projected links,
  // saving original values in Uinv[0] to be reset at end
  if (project == 1) {
    FORALLSITES(i, s) {
      mat_copy_f(&(s->linkf[dir]), &(Uinv[0][i]));
      polar(&(s->linkf[dir]), &tmat, &tmat2);
      mat_copy_f(&tmat, &(s->linkf[dir]));
    }
  }

  switch(dir) {
    case XUP:
      len = nx;
      norm = (double)(nt * ny * nz);
      break;
    case YUP:
      len = ny;
      norm = (double)(nx * nt * nz);
      break;
    case ZUP:
      len = nz;
      norm = (double)(nx * ny * nt);
      break;
    case TUP:
      len = nt;
      norm = (double)(nx * ny * nz);
      break;
    case DIR_5:
      len = nt;
      norm = (double)(nx * ny * nz);
      break;
    default:
      printf("ploop: unrecognized direction %d\n", dir);
      fflush(stdout);
      terminate(1);
  }

  // Special case: len == 1
  if (len == 1) {
    *plpMod = 0.0;
    FORALLSITES(i, s) {
      plp = trace_f(&(s->linkf[dir]));
      CSUM(sum, plp);
      *plpMod += cabs(&plp);
    }
    g_complexsum(&sum);
    g_doublesum(plpMod);
    CDIVREAL(sum, norm, sum);
    *plpMod /= norm;

    if (project == 1) {       // Reset original links
      FORALLSITES(i, s)
        mat_copy_f(&(Uinv[0][i]), &(s->linkf[dir]));
    }

    return sum;
  }

  // Compute line by steadily shifting links to hyperplane 0
  FORALLSITES(i, s)
    mat_copy_f(&(s->linkf[dir]), &(tempmat[i]));

  for (t = 1; t < len; t++) {
    shiftmat(tempmat, tempmat2, goffset[dir]);
    FORALLSITES(i, s) {
      j = s->x;
      switch(dir) {
        case YUP:   j = s->y; break;
        case ZUP:   j = s->z; break;
        case TUP:   j = s->t; break;
        case DIR_5: j = s->t; break;
      }
      if (j != 0)
        continue;

      if (t == 1)
        mult_nn_f(&(s->linkf[dir]), &(tempmat[i]), &(staple[i]));
      else {
        mult_nn_f(&(staple[i]), &(tempmat[i]), &(tempmat2[i]));
        mat_copy_f(&(tempmat2[i]), &(staple[i]));
      }
    }
  }

  *plpMod = 0.0;
  FORALLSITES(i, s) {
    j = s->x;
    switch(dir) {
      case YUP:   j = s->y; break;
      case ZUP:   j = s->z; break;
      case TUP:   j = s->t; break;
      case DIR_5: j = s->t; break;
    }
    if (j != 0)
      continue;

    plp = trace_f(&(staple[i]));
    CSUM(sum, plp);
    *plpMod += cabs(&plp);
  }
  g_complexsum(&sum);
  g_doublesum(plpMod);
  CDIVREAL(sum, norm, sum);
  *plpMod /= norm;

  if (project == 1) {       // Reset original links
    FORALLSITES(i, s)
      mat_copy_f(&(Uinv[0][i]), &(s->linkf[dir]));
  }

  return sum;
}
// -----------------------------------------------------------------
