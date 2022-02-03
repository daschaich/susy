// -----------------------------------------------------------------
// Measure the phase of the pfaffian using gaussian elimination
// This provides a simpler serial check of the parallel computation
#include "susy_includes.h"
#define PFA_TOL 1e-12   // !!! The size of this can affect stability
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// This is now just a convenience routine to set up the fermion_op matrix
// 'in' is all zero except for one unit in[iter]
void matvec(Real *in, complex *out) {
  register site *s;
  int i, j, iter;

  // Copy complex vector into Twist_Fermion src
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  iter = 0;
  FORALLSITES(i, s) {
    clear_TF(&(src[i]));
    for (j = 0; j < DIMF; j++) {
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Fsite));
      iter++;
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Flink[0]));
      iter++;
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Flink[1]));
      iter++;
      if (in[iter] > 0.5)
        sum_matrix(&(Lambda[j]), &(src[i].Fplaq));
      iter++;
    }
  }
#ifdef DEBUG_CHECK
  // Check that we didn't miss any components of the input vector
  int Ndat = 4 * DIMF;
  if (iter != sites_on_node * Ndat) {
    printf("phase: cycled over %d of %d input components\n",
           iter, sites_on_node * Ndat);
    terminate(1);
  }
#endif

  fermion_op(src, res, PLUS);     // D

  // Copy the resulting Twist_Fermion res back to complex vector y
  iter = 0;
  FORALLSITES(i, s) {
    for (j = 0; j < DIMF; j++) {
      out[iter] = complextrace_nn(&(res[i].Fsite), &(Lambda[j]));
      iter++;
      out[iter] = complextrace_nn(&(res[i].Flink[0]), &(Lambda[j]));
      iter++;
      out[iter] = complextrace_nn(&(res[i].Flink[1]), &(Lambda[j]));
      iter++;
      out[iter] = complextrace_nn(&(res[i].Fplaq), &(Lambda[j]));
      iter++;
    }
  }
#ifdef DEBUG_CHECK
  // Check that we didn't miss any components of the output vector
  if (iter != sites_on_node * Ndat) {
    printf("phase: cycled over %d of %d output components\n",
           iter, sites_on_node * Ndat);
    terminate(1);
  }
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef PHASE
void phase() {
  register int i, j, k, pivot;
  int Ndat = 4 * DIMF, tot_dat = volume * Ndat, interchange = 1;
  Real *col = malloc(sizeof *col * tot_dat);
  double phase, log_mag, tr, maxSq;
  complex scale;
  complex **D = malloc(sizeof(complex*) * tot_dat);

  // This is serial code
  if (this_node != 0) {
    printf("ERROR: run this thing in serial!\n");
    terminate(1);
  }

#ifdef DEBUG_CHECK
  // In the past the size of PFA_TOL has affected stability
  node0_printf("Running serial pfaffian computation with tolerance %.4g\n",
               PFA_TOL);
#endif

  // Make sure col has only one non-zero component
  for (i = 0; i < tot_dat; i++)
    col[i] = 0.0;

  // Allocate and construct fermion operator D
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  for (i = 0; i < tot_dat; i++) {
    D[i] = malloc(sizeof(complex) * tot_dat);
    col[i] = 1.0;
    matvec(col, D[i]);
    col[i] = 0.0;
  }
  free(col);    // Done with this

#ifdef DEBUG_CHECK
  // Check anti-symmetry of D
  int count = 0;
  for (i = 0; i < tot_dat; i++) {
    if (fabs(D[i][i].real) > PFA_TOL || fabs(D[i][i].imag) > PFA_TOL)
      node0_printf("  (%.4g, %.4g)\n", D[i][i].real, D[i][i].imag);
    for (j = i + 1; j < tot_dat; j++) {
      if (fabs(D[i][j].real + D[j][i].real) > PFA_TOL
       || fabs(D[i][j].imag + D[j][i].imag) > PFA_TOL) {
        printf("phase: D[%d][%d] = (%.4g, %.4g) but ",
               i, j, D[i][j].real, D[i][j].imag);
        printf("D[%d][%d] = (%.4g, %.4g)\n",
               j, i, D[j][i].real, D[j][i].imag);
      }
    }

    // Make sure every column of D is non-vanishing
    tr = cabs_sq(&(D[i][0]));
    for (j = 1; j < tot_dat; j++)
      tr += cabs_sq(&(D[i][j]));
    if (tr < PFA_TOL)
      printf("phase: Column %d vanishes: %.4g\n", i, tr);

    // Count number of non-zero elements in lower triangle
    for (j = i; j < tot_dat; j++) {
      if (cabs_sq(&(D[i][j])) > PFA_TOL)
        count++;
    }
  }
  printf("%d of %d elements of the full fermion matrix are non-zero\n",
         2 * count, tot_dat * tot_dat);
#endif

  // Loop over all pairs of columns of D
  for (i = 0; i < tot_dat - 1; i += 2) {
    // Find row of the largest element in column i, to use as pivot
    pivot = i + 1;
    maxSq = cabs_sq(&(D[i][pivot]));
    for (j = i + 2; j < tot_dat; j++) {
      tr = cabs_sq(&(D[i][j]));
      if (tr > maxSq) {
        maxSq = tr;
        pivot = j;
      }
    }
#ifdef DEBUG_CHECK
    if (maxSq < PFA_TOL) {    // We have a problem
      printf("phase: maxSq = %.4g for i = %d:\n", maxSq, i);
      for (j = 0; j < tot_dat; j++)
        printf("       (%.4g, %.4g)\n", D[i][j].real, D[i][j].imag);
      terminate(1);
    }
#endif

    // If necessary, interchange row (i + 1) with pivot row
    // Then interchange column (i + 1) with pivot column
    node0_printf("Columns %d--%d of %d: ", i + 1, i + 2, tot_dat);
    if (pivot != i + 1) {
      node0_printf("pivoting %d<-->%d to obtain ", i + 1, pivot);
      interchange *= -1;

      for (j = 0; j < tot_dat; j++) {
        set_complex_equal(&(D[j][i + 1]), &scale);
        set_complex_equal(&(D[j][pivot]), &(D[j][i + 1]));
        set_complex_equal(&scale, &(D[j][pivot]));
      }
      for (j = 0; j < tot_dat; j++) {
        set_complex_equal(&(D[i + 1][j]), &scale);
        set_complex_equal(&(D[pivot][j]), &(D[i + 1][j]));
        set_complex_equal(&scale, &(D[pivot][j]));
      }
    }
    else
      node0_printf("no need to pivot to obtain ");
    node0_printf("%.16g %.16g\n", D[i][i + 1].real, D[i][i + 1].imag);

    // Now that maximum elements are in D[i][i + 1] and D[i + 1][i],
    // progressively zero out elements of columns D[i] & D[i + 1]
    // Then zero out elements of rows D[:][i] & D[:][i + 1]
    for (j = i + 2; j < tot_dat; j++) {
      if (cabs_sq(&(D[i][i + 1])) < PFA_TOL) {
        printf("phase: can't divide by D[%d][%d] = (%.4g, %.4g)\n",
               i, i + 1, D[i][i + 1].real, D[i][i + 1].imag);
        terminate(1);
      }
      // scale = D[i][j] / D[i][i + 1]
      CDIV(D[i][j], D[i][i + 1], scale);

      for (k = 0; k < tot_dat; k++) {
        // D[k][j] -= scale * D[k][i + 1]
        CMULDIF(scale, D[k][i + 1], D[k][j]);
        // D[j][k] -= scale * D[i + 1][k]
        CMULDIF(scale, D[i + 1][k], D[j][k]);
      }
    }
    for (j = i + 2; j < tot_dat; j++) {
      // scale = D[i + 1][j] / D[i + 1][i]
      if (cabs_sq(&(D[i + 1][i])) < PFA_TOL) {
        printf("phase: can't divide by D[%d][%d] = (%.4g, %.4g)\n",
               i + 1, i, D[i + 1][i].real, D[i + 1][i].imag);
        terminate(1);
      }
      CDIV(D[i + 1][j], D[i + 1][i], scale);

      for (k = 0; k < tot_dat; k++) {
        // D[k][j] -= scale * D[k][i]
        CMULDIF(scale, D[k][i], D[k][j]);
        // D[j][k] -= scale * D[i][k]
        CMULDIF(scale, D[i][k], D[j][k]);
      }
    }
#ifdef DEBUG_CHECK
    // Make sure every column of D is still non-vanishing
    for (k = 0; k < tot_dat; k++) {
      tr = cabs_sq(&(D[k][0]));
      for (j = 1; j < tot_dat; j++)
        tr += cabs_sq(&(D[k][j]));
      if (tr < PFA_TOL)
        printf("Column %d vanishes after round %d: %.4g\n", k, i, tr);
    }

    // Check that D is still anti-symmetric
    for (k = 0; k < tot_dat; k++) {
      if (fabs(D[k][k].real) > PFA_TOL || fabs(D[k][k].imag) > PFA_TOL)
        node0_printf("  (%.4g, %.4g)\n", D[k][k].real, D[k][k].imag);
      for (j = k + 1; j < tot_dat; j++) {
        if (fabs(D[k][j].real + D[j][k].real) > PFA_TOL
         || fabs(D[k][j].imag + D[j][k].imag) > PFA_TOL) {
          printf("phase: D[%d][%d] = (%.4g, %.4g) but ",
                 k, j, D[k][j].real, D[k][j].imag);
          printf("D[%d][%d] = (%.4g, %.4g)\n",
                 j, k, D[j][k].real, D[j][k].imag);
          terminate(1);
        }
      }
    }
#endif

    // !!! Restore exact anti-symmetry after each round
    // Improves stability, perhaps related to roundoff accumulation?
    for (k = 0; k < tot_dat; k++) {
      for (j = k; j < tot_dat; j++) {
        CDIF(D[k][j], D[j][k]);
        CMULREAL(D[k][j], 0.5, D[k][j]);
        CNEGATE(D[k][j], D[j][k]);
      }
    }
  }

  // Now we can just compute the pfaffian from D[i][i + 1],
  // accumulating the phase and (the log of) the magnitude
  phase = 0.0;
  log_mag = 0.0;
  for (i = 0; i < tot_dat; i += 2) {
    phase += carg(&(D[i][i + 1]));
    log_mag += log(cabs_sq(&(D[i][i + 1]))) / 2.0;
    free(D[i]);
  }
  free(D);

  // Account for potential negative sign from interchanges
  if (interchange > 1)
    tr = phase;
  else
    tr = phase + PI;
  // Following the C++ code, keep phase in [0, 2pi)
  tr = fmod(phase, TWOPI);
  if (tr < 0)
    phase = tr + TWOPI;
  else
    phase = tr;

  node0_printf("PFAFF %.8g %.8g %.8g %.8g\n", log_mag, phase,
               fabs(cos(phase)), fabs(sin(phase)));
  fflush(stdout);
}
#endif
// -----------------------------------------------------------------
