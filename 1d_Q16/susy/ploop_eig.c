// -----------------------------------------------------------------
// Print out all eigenvalues of Wilson line on corresponding slice
// Might as well continue to return Wilson line itself and its magnitude
// Use repeated single-timeslice gathers to construct Wilson line
// Use tempmat, tempmat2, staple and UpsiU[0] for temporary storage
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
complex ploop_eig() {
  register int i, index = node_index(0);
  register site *s;
  char N = 'N';
  int j, k, t;
  int size = NCOL, stat = 0, unit = 1, doub = 2 * NCOL;
  double mag;
  double *eigs = malloc(sizeof *eigs * 2 * NCOL);
  double *dum = malloc(sizeof *dum * 2);
  complex ave, plp = cmplx(0.0, 0.0), tc;
  complex *ceigs = malloc(sizeof *ceigs * NCOL);
  matrix tmat, tmat2;
  
  // Special case: nt == 1
  if (nt == 1) {
    printf("ploop_eig: not yet set up for nt=1\n");
    fflush(stdout);
    terminate(1);
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
  
  // Kill off roundoff for NCOL=2
  if (fabs(plp.imag) < IMAG_TOL)
    plp.imag = 0.0;
  
  // Convert Wilson line to column-major double array expected by LAPACK
  for (j = 0; j < NCOL; j++) {
    for (k = 0; k < NCOL; k++) {
      store[2 * (k + NCOL * j)] = tmat.e[j][k].real;
      store[2 * (k + NCOL * j) + 1] = tmat.e[j][k].imag;
    }
  }
  
  // Diagonalize Wilson line using LAPACK
  zgeev_(&N, &N, &size, store, &size, eigs,
         dum, &unit, dum, &unit, work, &doub, work, &stat);
  
  // Extract eigenvalues from LAPACK output
  // Accumulate real and imaginary parts to compute average phase:
  // https://en.wikipedia.org/wiki/Mean_of_circular_quantities
  // If the links should be in SU(N), check that eig magnitudes are unity
  ave = cmplx(0.0, 0.0);
  for (j = 0; j < NCOL; j++) {
    ceigs[j].real = eigs[2 * j];
    ceigs[j].imag = eigs[2 * j + 1];
    mag = cabs(&(ceigs[j]));
    if (fabs(mag - 1.0) > IMAG_TOL) {
      printf("WARNING: Non-unitary Wilson line eigenvalue: ");
      printf("%.4g %.4g --> %.4g %.4g\n",
             ceigs[j].real, ceigs[j].imag, mag, carg(&(ceigs[j])));
    }
    CSUM(ave, ceigs[j]);
  }
  
  // Make sure phases haven't all canceled out
  // Might as well move ave back onto unit circle if possible
  // (No need to average ave by NCOL---doesn't affect phase)
  mag = cabs(&ave);
  if (fabs(mag) < IMAG_TOL) {
    printf("ERROR: phases cancelled out, can't average\n");
    fflush(stdout);
    terminate(1);
  }
  CDIVREAL(ave, mag, ave);
  
  // Divide each eigenvalue by ave to extract relative phase, and print
  node0_printf("LINES_EIG ");
  for (j = 0; j < NCOL; j++) {
    CDIV(ceigs[j], ave, tc);
    printf(" %.4g", carg(&tc));     // Produces result in [-pi, pi)
  }
  node0_printf("\n");
  
  // Don't forget about the Wilson line (trace) itself
  // Might as well print this out as well
  
  node0_printf("LINES_TR   %.4g %.4g\n", plp.real, plp.imag);
  
  free(eigs);
  free(dum);
  free(ceigs);
  return plp;
}
// -----------------------------------------------------------------
