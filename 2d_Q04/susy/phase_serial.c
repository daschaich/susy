// -----------------------------------------------------------------
// Measure the phase of the pfaffian using gaussian elimination
// This provides a simpler serial check of the parallel computation
#include "susy_includes.h"
#define PFA_TOL 1e-12
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// This is now just a convenience routine to set up the fermion_op matrix
// 'in' is all zero except for one unit in[iter]
void matvec(Real *in, complex *out) {
  register site *s;
  int i, j, mu, iter;

  // Copy complex vector into Twist_Fermion src
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  iter = 0;
  FORALLSITES(i, s) {
    clear_TF(&(src[i]));
    for (j = 0; j < DIMF; j++) {
      if (in[iter] > 0.0)
        sum_matrix(&(Lambda[j]), &(src[i].Fsite));
      iter++;
      if (in[iter] > 0.0)
        sum_matrix(&(Lambda[j]), &(src[i].Flink[0]));
      iter++;
      if (in[iter] > 0.0)
        sum_matrix(&(Lambda[j]), &(src[i].Flink[1]));
      iter++;
      if (in[iter] > 0.0)
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

  fermion_op(src, res, PLUS);    // D

  // Copy the resulting Twist_Fermion res back to complex vector y
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  iter = 0;
  FORALLSITES(i, s) {
    for (j = 0; j < DIMF; j++) {
      out[iter] = complextrace_nn(&(res[i].Fsite), &(Lambda[j]));
      iter++;
      FORALLDIR(mu) {
        out[iter] = complextrace_nn(&(res[i].Flink[mu]), &(Lambda[j]));
        iter++;
      }
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
  int Ndat = 4 * DIMF, interchange = 1;
  Real *col = malloc(volume * Ndat * sizeof(*col));
  double phase, log_mag, tr, maxSq;
  complex scale;
  complex **D = malloc(volume * Ndat * sizeof(complex*));

  // This is serial code
  if (this_node != 0) {
    printf("ERROR: run this thing in serial!\n");
    terminate(1);
  }

  // Make sure col has only one non-zero component
  for (i = 0; i < volume * Ndat; i++)
    col[i] = 0.0;

  // Allocate and construct fermion operator D
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  for (i = 0; i < volume * Ndat; i++) {
    D[i] = malloc(volume * Ndat * sizeof(complex));
    col[i] = 1.0;
    matvec(col, D[i]);
    col[i] = 0.0;
  }
  free(col);    // Done with this

#ifdef DEBUG_CHECK
  // Check anti-symmetry of D
  int count = 0;
  for (i = 0; i < volume * Ndat; i++) {
    for (j = i + 1; j < volume * Ndat; j++) {
      if (D[i][j].real + D[j][i].real > PFA_TOL
       || D[i][j].imag + D[j][i].imag > PFA_TOL) {
        printf("phase: D[%d][%d] = (%.4g, %.4g) but ",
               i, j, D[i][j].real, D[i][j].imag);
        printf("D[%d][%d] = (%.4g, %.4g)\n",
               j, i, D[j][i].real, D[j][i].imag);
      }
    }

    // Make sure every column of D is non-vanishing
    tr = cabs_sq(&(D[i][0]));
    for (j = 1; j < volume * Ndat; j++)
      tr += cabs_sq(&(D[i][j]));
    if (tr < PFA_TOL)
      printf("phase: Column %d vanishes: %.4g\n", i, tr);

    // Count number of non-zero elements in upper triangle
    for (j = i; j < volume * Ndat; j++) {
      if (cabs_sq(&(D[i][j])) > PFA_TOL)
        count++;
    }
  }
  printf("%d of %d elements of the full fermion matrix are non-zero\n",
         2 * count, volume * Ndat * volume * Ndat);
#endif

  // Loop over all pairs of columns of D
  for (i = 0; i < volume * Ndat - 1; i += 2) {
    // Find row of the largest element in column i, to use as pivot
    pivot = i + 1;
    maxSq = cabs_sq(&(D[i][pivot]));
    for (j = i + 2; j < volume * Ndat; j++) {
      tr = cabs_sq(&(D[i][j]));
      if (tr > maxSq) {
        maxSq = tr;
        pivot = j;
      }
    }
#ifdef DEBUG_CHECK
    if (maxSq < PFA_TOL) {    // We have a problem
      printf("phase: maxSq = %.4g for i = %d:\n", maxSq, i);
      for (j = 0; j < volume * Ndat; j++)
        printf("       (%.4g, %.4g)\n", D[i][j].real, D[i][j].imag);
      terminate(1);
    }
#endif

    // If necessary, interchange row (i + 1) with pivot row
    // Then interchange column (i + 1) with pivot column
    node0_printf("Columns %d--%d of %d: ", i + 1, i + 2, volume * Ndat);
    if (pivot != i + 1) {
      node0_printf("pivoting %d<-->%d to obtain ", i + 1, pivot);
      interchange *= -1;

      for (j = 0; j < volume * Ndat; j++) {
        set_complex_equal(&(D[j][i + 1]), &scale);
        set_complex_equal(&(D[j][pivot]), &(D[j][i + 1]));
        set_complex_equal(&scale, &(D[j][pivot]));
      }
      for (j = 0; j < volume * Ndat; j++) {
        set_complex_equal(&(D[i + 1][j]), &scale);
        set_complex_equal(&(D[pivot][j]), &(D[i + 1][j]));
        set_complex_equal(&scale, &(D[pivot][j]));
      }
    }
    else
      node0_printf("no need to pivot to obtain ", i + 1, pivot);
    node0_printf("%.16g %.16g\n", D[i][i + 1].real, D[i][i + 1].imag);

    // Now that maximum elements are in D[i][i + 1] and D[i + 1][i],
    // progressively zero out elements of columns D[i] & D[i + 1]
    // Then zero out elements of rows D[:][i] & D[:][i + 1]
    for (j = i + 2; j < volume * Ndat; j++) {
      if (cabs_sq(&(D[i][i + 1])) < PFA_TOL) {
        printf("phase: can't divide by D[%d][%d] = (%.4g, %.4g)\n",
               i, i + 1, D[i][i + 1].real, D[i][i + 1].imag);
        terminate(1);
      }
      // scale = D[i][j] / D[i][i + 1]
      CDIV(D[i][j], D[i][i + 1], scale);

      for (k = 0; k < volume * Ndat; k++) {
        // D[k][j] -= scale * D[k][i + 1]
        CMULDIF(scale, D[k][i + 1], D[k][j]);
        // D[j][k] -= scale * D[i + 1][k]
        CMULDIF(scale, D[i + 1][k], D[j][k]);
      }
    }
    for (j = i + 2; j < volume * Ndat; j++) {
      // scale = D[i + 1][j] / D[i + 1][i]
      if (cabs_sq(&(D[i + 1][i])) < PFA_TOL) {
        printf("phase: can't divide by D[%d][%d] = (%.4g, %.4g)\n",
               i + 1, i, D[i + 1][i].real, D[i + 1][i].imag);
        terminate(1);
      }
      CDIV(D[i + 1][j], D[i + 1][i], scale);

      for (k = 0; k < volume * Ndat; k++) {
        // D[k][j] -= scale * D[k][i]
        CMULDIF(scale, D[k][i], D[k][j]);
        // D[j][k] -= scale * D[i][k]
        CMULDIF(scale, D[i][k], D[j][k]);
      }
    }

#ifdef DEBUG_CHECK
    // Make sure every column of D is still non-vanishing
    for (k = 0; k < volume * Ndat; k++) {
      tr = cabs_sq(&(D[k][0]));
      for (j = 1; j < volume * Ndat; j++)
        tr += cabs_sq(&(D[k][j]));
      if (tr < PFA_TOL)
        printf("Column %d vanishes after round %d: %.4g\n", k, i, tr);
    }
#endif
  }

#ifdef DEBUG_CHECK
  // Check that D is still anti-symmetric, and now tridiagonal
  for (i = 0; i < volume * Ndat; i++) {
    for (j = i + 1; j < volume * Ndat; j++) {
      if (D[i][j].real + D[j][i].real > PFA_TOL
       || D[i][j].imag + D[j][i].imag > PFA_TOL) {
        printf("phase: D[%d][%d] = (%.4g, %.4g) but ",
               i, j, D[i][j].real, D[i][j].imag);
        printf("D[%d][%d] = (%.4g, %.4g)\n",
               j, i, D[j][i].real, D[j][i].imag);
      }
    }
  }
#endif

  // Now we can just compute the pfaffian from D[i][i + 1],
  // accumulating the phase and (the log of) the magnitude
  phase = 0.0;
  log_mag = 0.0;
  for (i = 0; i < volume * Ndat; i += 2) {
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
