// -----------------------------------------------------------------
// Measure the eigenvalues of DSq using LAPACK
// This provides a simpler serial check of the PRIMME computation
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// This is now just a convenience routine to set up the fermion_op matrix
// 'in' is all zero except for one unit in[iter]
void matvec(Real *in, complex *out) {
  register site *s;
  int i, j, iter;
  msg_tag *mtag0, *mtag1, *mtag2;
  matrix *plaq23, *plaq13, *plaq12;

  // Copy complex vector into Twist_Fermion src
  // Each Twist_Fermion has Ndat = 16DIMF non-trivial complex components
  iter = 0;
  FORALLSITES(i, s) {
    clear_TF(&(src[i]));
    clear_mat(&(plaq_src[7][i]));
    clear_mat(&(plaq_src[5][i]));
    clear_mat(&(plaq_src[4][i]));
    for (j = 0; j < DIMF; j++) {
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Fsite));
      iter++;
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Flink[4]));
      iter++;

      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Flink[0]));
      iter++;
      // 0, 4 --> 3
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Fplaq[3]));
      iter++;
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Flink[1]));
      iter++;
      // 1, 4 --> 6
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Fplaq[6]));
      iter++;
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Flink[2]));
      iter++;
      // 2, 4 --> 8
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Fplaq[8]));
      iter++;
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Flink[3]));
      iter++;
      // 3, 4 --> 9
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Fplaq[9]));
      iter++;

      // 0, 1 --> 0
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Fplaq[0]));
      iter++;
      // 2, 3 --> 7
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(plaq_src[7][i]));
      iter++;

      // 0, 2 --> 1
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Fplaq[1]));
      iter++;
      // 1, 3 --> 5
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(plaq_src[5][i]));
      iter++;

      // 0, 3 --> 2
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Fplaq[2]));
      iter++;
      // 1, 2 --> 4
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(plaq_src[4][i]));
      iter++;
    }
  }

  // Gather plaq_src[7] (2, 3) from x - 0 - 1 (gather path 23),
  // plaq_src[5] (1, 3) from x - 0 - 2 (gather path 31)
  // and plaq_src[4] (1, 2) from x - 0 - 3 (gather path 35)
  mtag0 = start_gather_field(plaq_src[7], sizeof(matrix),
                             23, EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_field(plaq_src[5], sizeof(matrix),
                             31, EVENANDODD, gen_pt[1]);
  mtag2 = start_gather_field(plaq_src[4], sizeof(matrix),
                             35, EVENANDODD, gen_pt[2]);

  wait_gather(mtag0);
  wait_gather(mtag1);
  wait_gather(mtag2);
  FORALLSITES(i, s) {
    mat_copy((matrix *)(gen_pt[0][i]), &(src[i].Fplaq[7]));  // 2, 3
    mat_copy((matrix *)(gen_pt[1][i]), &(src[i].Fplaq[5]));  // 1, 3
    mat_copy((matrix *)(gen_pt[2][i]), &(src[i].Fplaq[4]));  // 1, 2
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);
  cleanup_gather(mtag2);

#ifdef DEBUG_CHECK
  // Check that we didn't miss any components of the input vector
  int Ndat = 16 * DIMF;
  if (iter != sites_on_node * Ndat) {
    printf("eig: cycled over %d of %d input components\n",
           iter, sites_on_node * Ndat);
    terminate(1);
  }
#endif

  DSq(src, res);     // DSq

  // Copy the resulting Twist_Fermion res back to complex vector y
  // Gather plaq_src[7] (2, 3) from x + 0 + 1 (gather path 22),
  // plaq_src[5] (1, 3) from x + 0 + 2 (gather path 30)
  // and plaq_src[4] (1, 2) from x + 0 + 3 (gather path 34)
  FORALLSITES(i, s) {
    mat_copy(&(res[i].Fplaq[7]), &(plaq_src[7][i]));
    mat_copy(&(res[i].Fplaq[5]), &(plaq_src[5][i]));
    mat_copy(&(res[i].Fplaq[4]), &(plaq_src[4][i]));
  }
  mtag0 = start_gather_field(plaq_src[7], sizeof(matrix),
                             22, EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_field(plaq_src[5], sizeof(matrix),
                             30, EVENANDODD, gen_pt[1]);
  mtag2 = start_gather_field(plaq_src[4], sizeof(matrix),
                             34, EVENANDODD, gen_pt[2]);

  wait_gather(mtag0);
  wait_gather(mtag1);
  wait_gather(mtag2);
  iter = 0;
  FORALLSITES(i, s) {
    plaq23 = (matrix *)(gen_pt[0][i]);
    plaq13 = (matrix *)(gen_pt[1][i]);
    plaq12 = (matrix *)(gen_pt[2][i]);
    for (j = 0; j < DIMF; j++) {
      out[iter] = complextrace_nn(&(res[i].Fsite), &(Lambda[j]));
      iter++;
      out[iter] = complextrace_nn(&(res[i].Flink[4]), &(Lambda[j]));
      iter++;

      out[iter] = complextrace_nn(&(res[i].Flink[0]), &(Lambda[j]));
      iter++;
      // 0, 4 --> 3
      out[iter] = complextrace_nn(&(res[i].Fplaq[3]), &(Lambda[j]));
      iter++;
      out[iter] = complextrace_nn(&(res[i].Flink[1]), &(Lambda[j]));
      iter++;
      // 1, 4 --> 6
      out[iter] = complextrace_nn(&(res[i].Fplaq[6]), &(Lambda[j]));
      iter++;
      out[iter] = complextrace_nn(&(res[i].Flink[2]), &(Lambda[j]));
      iter++;
      // 2, 4 --> 8
      out[iter] = complextrace_nn(&(res[i].Fplaq[8]), &(Lambda[j]));
      iter++;
      out[iter] = complextrace_nn(&(res[i].Flink[3]), &(Lambda[j]));
      iter++;
      // 3, 4 --> 9
      out[iter] = complextrace_nn(&(res[i].Fplaq[9]), &(Lambda[j]));
      iter++;

      // 0, 1 --> 0
      out[iter] = complextrace_nn(&(res[i].Fplaq[0]), &(Lambda[j]));
      iter++;
      out[iter] = complextrace_nn(plaq23, &(Lambda[j]));
      iter++;

      // 0, 2 --> 1
      out[iter] = complextrace_nn(&(res[i].Fplaq[1]), &(Lambda[j]));
      iter++;
      out[iter] = complextrace_nn(plaq13, &(Lambda[j]));
      iter++;

      // 0, 3 --> 2
      out[iter] = complextrace_nn(&(res[i].Fplaq[2]), &(Lambda[j]));
      iter++;
      out[iter] = complextrace_nn(plaq12, &(Lambda[j]));
      iter++;
    }
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);
  cleanup_gather(mtag2);

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
  int Ndat = 16 * DIMF, tot_dat = volume * Ndat;
  int row, col, stat = 0, Nwork = 2 * tot_dat;
  Real *Dcol = malloc(sizeof *Dcol * tot_dat), size;
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

  // Lots of mallocing---enough to be worth warning the user about
  // (Recording DRwork as just 3tot_dat; the -2 should be decimal dust
  size = (Real)(sizeof(Real));
  size += (Real)((8.0 + 2.0 * tot_dat) * sizeof(double));
  size += (Real)(tot_dat * sizeof(complex));
  size *= tot_dat;
  node0_printf("\nWARNING: ");
  node0_printf("Mallocing %.1f MBytes for diagonalization\n", size / 1e6);
  normal_exit(0);

  // Make sure Dcol has only one non-zero component
  for (i = 0; i < tot_dat; i++)
    Dcol[i] = 0.0;

  // Allocate and construct fermion operator D
  // Each Twist_Fermion has Ndat = 16DIMF non-trivial complex components
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
