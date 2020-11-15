// -----------------------------------------------------------------
// Measure the eigenvalues of DSq using LAPACK
// This provides a simpler serial check of the PRIMME computation
// Note that we work with DIMF=NCOL^2-1 rather than NCOL^2
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// This is now just a convenience routine to set up the fermion_op matrix
// 'in' is all zero except for one unit in[iter]
void matvec(Real *in, complex *out) {
  register site *s;
  int i, j, k, iter;

  // Copy complex vector into matrix* src
  // with NFERMION * DIMF non-trivial complex components
  iter = 0;
  FORALLSITES(i, s) {
    for (k = 0; k < NFERMION; k++) {
      clear_mat(&(src[k][i]));
      for (j = 0; j < DIMF; j++) {
        if (in[iter] > 0.5)
          sum_matrix(&(Lambda[j]), &(src[k][i]));
        iter++;
      }
    }
  }
#ifdef DEBUG_CHECK
  // Check that we didn't miss any components of the input vector
  int Ndat = NFERMION * DIMF;
  if (iter != sites_on_node * Ndat) {
    printf("eig: cycled over %d of %d input components\n",
           iter, sites_on_node * Ndat);
    terminate(1);
  }
#endif

  DSq(src, res);     // DSq

  // Copy the resulting matrix* res back to complex vector y
  iter = 0;
  FORALLSITES(i, s) {
    for (k = 0; k < NFERMION; k++) {
      for (j = 0; j < DIMF; j++) {
        out[iter] = complextrace_nn(&(res[k][i]), &(Lambda[j]));
        iter++;
      }
    }
  }
#ifdef DEBUG_CHECK
  // Check that we didn't miss any components of the output vector
  if (iter != sites_on_node * Ndat) {
    printf("eig: cycled over %d of %d output components\n",
           iter, sites_on_node * Ndat);
    terminate(1);
  }
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef EIG
void eig() {
  register int i;
  char N = 'N';     // Ask LAPACK only for eigenvalues
  char U = 'U';     // Have LAPACK store upper triangle of U.Ubar
  int Ndat = NFERMION * DIMF, tot_dat = nt * Ndat;
  int row, col, stat = 0, Nwork = 2 * tot_dat;
  Real *Dcol = malloc(sizeof *Dcol * tot_dat);
  double *Dstore = malloc(sizeof *store * 2 * tot_dat * tot_dat);
  double *Deigs = malloc(sizeof *Deigs * tot_dat);
  double *Dwork = malloc(sizeof *Dwork * 4 * tot_dat);
  double *DRwork = malloc(sizeof *DRwork * (3 * tot_dat - 2));
  complex **D = malloc(sizeof(complex*) * tot_dat);

  // This is serial code
  if (this_node != 0) {
    printf("ERROR: run this thing in serial!\n");
    terminate(1);
  }

  // Make sure Dcol has only one non-zero component
  for (i = 0; i < tot_dat; i++)
    Dcol[i] = 0.0;

  // Allocate and construct fermion operator D
  // Working with DIMF=NCOL^2-1 rather than NCOL^2
  for (i = 0; i < tot_dat; i++) {
    D[i] = malloc(sizeof(complex) * tot_dat);
    Dcol[i] = 1.0;
    matvec(Dcol, D[i]);
    Dcol[i] = 0.0;
  }
  free(Dcol);    // Done with this

#ifdef DEBUG_CHECK
  // Check anti-symmetry of D
  int count = 0;
  for (i = 0; i < tot_dat; i++) {
    if (fabs(D[i][i].real) > EIG_TOL || fabs(D[i][i].imag) > EIG_TOL)
      node0_printf("  (%.4g, %.4g)\n", D[i][i].real, D[i][i].imag);
    for (j = i + 1; j < tot_dat; j++) {
      if (fabs(D[i][j].real + D[j][i].real) > EIG_TOL
       || fabs(D[i][j].imag + D[j][i].imag) > EIG_TOL) {
        printf("eig: D[%d][%d] = (%.4g, %.4g) but ",
               i, j, D[i][j].real, D[i][j].imag);
        printf("D[%d][%d] = (%.4g, %.4g)\n",
               j, i, D[j][i].real, D[j][i].imag);
      }
    }

    // Make sure every column of D is non-vanishing
    tr = cabs_sq(&(D[i][0]));
    for (j = 1; j < tot_dat; j++)
      tr += cabs_sq(&(D[i][j]));
    if (tr < EIG_TOL)
      printf("eig: Column %d vanishes: %.4g\n", i, tr);

    // Count number of non-zero elements in lower triangle
    for (j = i; j < tot_dat; j++) {
      if (cabs_sq(&(D[i][j])) > EIG_TOL)
        count++;
    }
  }
  printf("%d of %d elements of the full fermion matrix are non-zero\n",
         2 * count, tot_dat * tot_dat);
#endif

  // Now D is the fermion matrix, so let's feed it to LAPACK...
  // Convert X[j] to column-major double array used by LAPACK
  for (row = 0; row < tot_dat; row++) {
    for (col = 0; col < tot_dat; col++) {
      Dstore[2 * (col * tot_dat + row)] = D[row][col].real;
      Dstore[2 * (col * tot_dat + row) + 1] = D[row][col].imag;
    }
  }

  // Compute eigenvalues and eigenvectors of X[j]
  zheev_(&N, &U, &tot_dat, Dstore, &tot_dat, Deigs,
         Dwork, &Nwork, DRwork, &stat);

  if (stat != 0)
    printf("WARNING: Non-zero return from LAPACK\n");

  // Now the eigenvalues of DSq should just be the squares
  // of the non-zero elements of the non-zero 2x2 diagonal blocks
  for (i = 0; i < tot_dat; i++) {
    node0_printf("EIG %.8g\n", Deigs[i]);
    fflush(stdout);
    free(D[i]);
  }
  free(D);
  free(Dstore);
  free(Deigs);
  free(Dwork);
  free(DRwork);
}
#endif
// -----------------------------------------------------------------
