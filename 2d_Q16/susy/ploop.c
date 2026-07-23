// -----------------------------------------------------------------
// Evaluate the Polyakov loop using repeated single-timeslice gathers
// Use tempmat, tempmat2, staple and UpsiU[0] for temporary storage
#include "susy_includes.h"

// dir tells us the direction of the product: XUP, TUP or DIR_3
// project == 1 tells us to consider links unitarized via polar projection
complex ploop(int dir, int project, double *plpMod) {
  register int i;
  register site *s;
  int j, t, len = nt;
  double norm = 0.0;
  complex sum  = cmplx(0.0, 0.0), plp;
  matrix tmat, tmat2;

  // Optionally consider polar-projected links,
  // saving original values in UpsiU[0] to be reset at end
  if (project == 1) {
    FORALLSITES(i, s) {
      mat_copy(&(s->link[dir]), &(UpsiU[0][i]));
      polar(&(s->link[dir]), &tmat, &tmat2);
      mat_copy(&tmat, &(s->link[dir]));
    }
  }

  switch(dir) {
    case XUP:
      len = nx;
      norm = 1.0 / (double)(nt);
      break;
    case TUP:
      len = nt;
      norm = 1.0 / (double)(nx);
      break;
    case DIR_3:
      len = nt;
      norm = 1.0 / (double)(nx);
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
      plp = trace(&(s->link[dir]));
      CSUM(sum, plp);
      *plpMod += cabs(&plp);
    }
    g_complexsum(&sum);
    g_doublesum(plpMod);
    CMULREAL(sum, norm, sum);
    *plpMod *= norm;

    if (project == 1) {       // Reset original links
      FORALLSITES(i, s)
        mat_copy(&(UpsiU[0][i]), &(s->link[dir]));
    }

    return sum;
  }

  // Compute line by steadily shifting links to hyperplane 0
  FORALLSITES(i, s)
    mat_copy(&(s->link[dir]), &(tempmat[i]));

  for (t = 1; t < len; t++) {
    shiftmat(tempmat, tempmat2, goffset[dir]);
    FORALLSITES(i, s) {
      j = s->x;
      switch(dir) {
        case TUP:   j = s->t; break;
        case DIR_3: j = s->t; break;
      }
      if (j != 0)
        continue;

      if (t == 1)
        mult_nn(&(s->link[dir]), &(tempmat[i]), &(staple[i]));
      else {
        mult_nn(&(staple[i]), &(tempmat[i]), &(tempmat2[i]));
        mat_copy(&(tempmat2[i]), &(staple[i]));
      }
    }
  }

  *plpMod = 0.0;
  FORALLSITES(i, s) {
    j = s->x;
    switch(dir) {
      case TUP:   j = s->t; break;
      case DIR_3: j = s->t; break;
    }
    if (j != 0)
      continue;

    plp = trace(&(staple[i]));
    CSUM(sum, plp);
    *plpMod += cabs(&plp);
  }
  g_complexsum(&sum);
  g_doublesum(plpMod);
  CMULREAL(sum, norm, sum);  //edited -------
  *plpMod *= norm;

  if (project == 1) {       // Reset original links
    FORALLSITES(i, s)
      mat_copy(&(UpsiU[0][i]), &(s->link[dir]));
  }

  return sum;
}
// -----------------------------------------------------------------
