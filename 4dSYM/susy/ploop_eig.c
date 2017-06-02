// -----------------------------------------------------------------
// Print out all eigenvalues of Wilson line on corresponding slice
// Might as well continue to return Wilson line itself and its magnitude
// Use repeated single-timeslice gathers to construct Wilson line
// Use tempmat, tempmat2, staple and UpsiU[0] for temporary storage
// CAUTION: Do not run with MPI!
#include "susy_includes.h"

// Arguments summarized in susy_includes.h (for -DEIG compilation...)
void zgeev_(char *doL, char *doR, int *N1, double *store, int *N2, double *eigs,
            double *dumL, int *NL, double *dumR, int *NR,
            double *work, int *Nwork, double *Rwork, int *stat);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dir tells us the direction of the product: XUP, YUP, ZUP, TUP or DIR_5
// project == 1 tells us to consider links unitarized via polar projection
complex ploop_eig(int dir, int project, double *plpMod) {
  register int i;
  register site *s;
  char N = 'N';
  int j, k, t, len = nt;
  int size = NCOL, stat = 0, unit = 1, doub = 2 * NCOL;
  double norm = 0.0, mag, ave_phase, diff, back;
  double *phase = malloc(NCOL * sizeof(*phase));
  double *eigs = malloc(2 * NCOL * sizeof(*eigs));
  double *dum = malloc(2 * sizeof(*dum));
  complex sum  = cmplx(0.0, 0.0), plp, tc;
  matrix tmat, tmat2;

  if (this_node != 0) {
    printf("ploop_eig: don't run in parallel\n");
    fflush(stdout);
    terminate(1);
  }

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
      plp = trace(&(s->link[dir]));
      CSUM(sum, plp);
      *plpMod += cabs(&plp);
    }
    g_complexsum(&sum);
    g_doublesum(plpMod);
    CDIVREAL(sum, norm, sum);
    *plpMod /= norm;

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
        case YUP:   j = s->y; break;
        case ZUP:   j = s->z; break;
        case TUP:   j = s->t; break;
        case DIR_5: j = s->t; break;
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
      case YUP:   j = s->y; break;
      case ZUP:   j = s->z; break;
      case TUP:   j = s->t; break;
      case DIR_5: j = s->t; break;
    }
    if (j != 0)
      continue;

    // Convert Wilson line to column-major double array expected by LAPACK
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
      store[2 * (k + NCOL * j)] = staple[i].e[j][k].real;
      store[2 * (k + NCOL * j) + 1] = staple[i].e[j][k].imag;
      }
    }

    // Diagonalize Wilson line using LAPACK
    zgeev_(&N, &N, &size, store, &size, eigs,
           dum, &unit, dum, &unit, work, &doub, work, &stat);

    // Convert eigenvalues to phases in [-pi, pi)
    // Accumulate real and imaginary parts to compute average phase:
    // http://en.wikipedia.org/wiki/Mean_of_circular_quantities
    // If the links should be in SU(N), check that eig magnitudes are unity
    sum = cmplx(0.0, 0.0);
    for (j = 0; j < NCOL; j++) {
      tc.real = eigs[2 * j];
      tc.imag = eigs[2 * j + 1];
      if (project == 1) {
        mag = cabs_sq(&tc);
        if (fabs(mag - 1.0) > IMAG_TOL) {
          printf("WARNING: Non-unitary Wilson line eigenvalue: ");
          printf("%d %d %d %d %d %.4g %.4g --> %.4g %.4g\n",
                 s->x, s->y, s->z, s->t, dir, tc.real, tc.imag, mag, phase[j]);
        }
      }
      CSUM(sum, tc);
      phase[j] = carg(&tc);       // Produces result in [-pi, pi)
    }
    ave_phase = carg(&sum);

    // Make sure phases haven't all canceled out
    mag = cabs_sq(&sum);
    if (mag < IMAG_TOL) {
      printf("ERROR: phases cancelled out, can't average\n");
      fflush(stdout);
      terminate(1);
    }

    // Print resulting phases as fluctuations about average
    // Include all four coords, even though the one matching dir will be zero
    if (project == 1)
      printf("LINES_POLAR_EIG ");
    else
      printf("LINES_EIG ");
    printf("%d %d %d %d %d", s->x, s->y, s->z, s->t, dir);
    for (j = 0; j < NCOL; j++) {
      // Make sure the fluctuation remains in [-pi, pi)
      diff = phase[j] - ave_phase;
      back = TWOPI - fabs(diff);
      if (fabs(diff) < fabs(back))
        printf(" %.4g", diff);
      else
        printf(" %.4g", back);
    }
    printf("\n");

    plp = trace(&(staple[i]));
    CSUM(sum, plp);
    *plpMod += cabs(&plp);
  }
  g_complexsum(&sum);
  g_doublesum(plpMod);
  CDIVREAL(sum, norm, sum);
  *plpMod /= norm;

  if (project == 1) {       // Reset original links
    FORALLSITES(i, s)
      mat_copy(&(UpsiU[0][i]), &(s->link[dir]));
  }

  free(eigs);
  free(dum);
  return sum;
}
// -----------------------------------------------------------------
