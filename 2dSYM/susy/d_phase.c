// -----------------------------------------------------------------
// Measure the phase of the pfaffian
#include "susy_includes.h"
#define PFA_TOL 1e-8
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wrapper for the fermion_op matvec
// Convert between complex vectors and Twist_Fermions
void matvec(complex *in, complex *out) {
  register site *s;
  int i, j, mu, nu, iter, Ndat = 4 * DIMF;

  // Copy complex vector into Twist_Fermion src
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  FORALLSITES(i, s) {
    iter = i * Ndat;
    for (j = 0; j < DIMF; j++) {
      set_complex_equal(&(in[iter]), &(src[i].Fsite.c[j]));
      iter++;
      for (mu = 0; mu < NUMLINK; mu++) {
        set_complex_equal(&(in[iter]), &(src[i].Flink[mu].c[j]));
        iter++;
        src[i].Fplaq[mu][mu].c[j] = cmplx(0.0, 0.0);  // Clear diagonal Fplaq
        for (nu = mu + 1; nu < NUMLINK; nu++) {
          set_complex_equal(&(in[iter]), &(src[i].Fplaq[mu][nu].c[j]));
          CNEGATE(src[i].Fplaq[mu][nu].c[j], src[i].Fplaq[nu][mu].c[j]);
          iter++;
        }
      }
    }
  }

  fermion_op(src, res, 1);    // D

  // Copy the resulting Twist_Fermion res back to complex vector y
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  FORALLSITES(i, s) {
    iter = i * Ndat;
    for (j = 0; j < DIMF; j++) {
      set_complex_equal(&(res[i].Fsite.c[j]), &(out[iter]));
      iter++;
      for (mu = 0; mu < NUMLINK; mu++) {
        set_complex_equal(&(res[i].Flink[mu].c[j]), &(out[iter]));
        iter++;
        for (nu = mu + 1; nu < NUMLINK; nu++) {
          set_complex_equal(&(res[i].Fplaq[mu][nu].c[j]), &(out[iter]));
          iter++;
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void d_phase() {
  register int i, j, k, pivot;
  int Ndat = 4 * DIMF, interchange = 1;
  double phase, log_mag, tr, maxSq;
  complex tc, scale;
  complex *col = malloc(volume * Ndat * sizeof(*col));
  complex **D = malloc(volume * Ndat * sizeof(complex*));

  // To be safe/explicit
  for (i = 0; i < volume * Ndat; i++)
    col[i] = cmplx(0.0, 0.0);

  // Allocate and construct fermion operator D
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  for (i = 0; i < volume * Ndat; i++) {
    D[i] = malloc(volume * Ndat * sizeof(complex));
    col[i] = cmplx(1.0, 0.0);
    matvec(col, D[i]);
    col[i] = cmplx(0.0, 0.0);
  }
  free(col);    // Done with this

//#ifdef DEBUG_CHECK
  // Check anti-symmetry of D
  int count = 0;
  for (i = 0; i < volume * Ndat; i++) {
    for (j = i + 1; j < volume * Ndat; j++) {
      if (D[i][j].real + D[j][i].real > PFA_TOL
       || D[i][j].imag + D[j][i].imag > PFA_TOL) {
        printf("d_phase: D[%d][%d] = (%.4g, %.4g) but ",
               i, j, D[i][j].real, D[i][j].imag);
        printf("D[%d][%d] = (%.4g, %.4g)\n",
               j, i, D[j][i].real, D[j][i].imag);
      }
    }

    // Make sure every column of D is non-vanishing
    tr = 0.0;
    for (j = 0; j < volume * Ndat; j++)
      tr += cabs_sq(&(D[i][j]));
    if (tr < PFA_TOL)
      printf("d_phase: Column %d vanishes: %.4g\n", i, tr);

    // Count number of non-zero elements in upper triangle
    for (j = i; j < volume * Ndat; j++) {
      if (cabs_sq(&(D[i][j])) > PFA_TOL)
        count++;
    }
  }
  printf("%d of %d elements of the full fermion matrix are non-zero\n",
         2 * count, volume * Ndat * volume * Ndat);
//#endif

  // Now that we have an explicit matrix, loop over all pairs of its columns
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
//    if (maxSq < PFA_TOL) {    // We have a problem
//      printf("d_phase: maxSq = %.4g for i = %d:\n", maxSq, i);
//      for (j = 0; j < volume * Ndat; j++)
//        printf("         (%.4g, %.4g)\n", D[i][j].real, D[i][j].imag);
//      terminate(1);
//    }

    // If necessary, interchange row (i + 1) with pivot row
    // Then interchange column (i + 1) with pivot column
    if (pivot != i + 1) {
      interchange *= -1;

      for (j = 0; j < volume * Ndat; j++) {
        set_complex_equal(&(D[j][i + 1]), &tc);
        set_complex_equal(&(D[j][pivot]), &(D[j][i + 1]));
        set_complex_equal(&tc, &(D[j][pivot]));
      }
      for (j = 0; j < volume * Ndat; j++) {
        set_complex_equal(&(D[i + 1][j]), &tc);
        set_complex_equal(&(D[pivot][j]), &(D[i + 1][j]));
        set_complex_equal(&tc, &(D[pivot][j]));
      }
    }

    // Now that maximum elements are in D[i][i + 1] and D[i + 1][i],
    // progressively zero out elements of columns D[i] & D[i + 1]
    // Then zero out elements of rows D[:][i] & D[:][i + 1]
    for (j = i + 2; j < volume * Ndat; j++) {
      // scale = D[i][j] / D[i][i + 1]
      if (cabs_sq(&(D[i][i + 1])) < PFA_TOL) {
        printf("d_phase: can't divide by D[%d][%d] = (%.4g, %.4g)\n",
               i, i + 1, D[i][i + 1].real, D[i][i + 1].imag);
        terminate(1);
      }
      CDIV(D[i][j], D[i][i + 1], scale);
      for (k = 0; k < volume * Ndat; k++) {
        CMUL(scale, D[k][i + 1], tc);
        CDIF(D[k][j], tc);    // D[k][j] -= scale * D[k][i + 1]

        CMUL(scale, D[i + 1][k], tc);
        CDIF(D[j][k], tc);    // D[j][k] -= scale * D[i + 1][k]
      }
    }
    for (j = i + 2; j < volume * Ndat; j++) {
      // scale = D[i + 1][j] / D[i + 1][i]
      if (cabs_sq(&(D[i + 1][i])) < PFA_TOL) {
        printf("d_phase: can't divide by D[%d][%d] = (%.4g, %.4g)\n",
               i + 1, i, D[i + 1][i].real, D[i + 1][i].imag);
        terminate(1);
      }
      CDIV(D[i + 1][j], D[i + 1][i], scale);
      for (k = 0; k < volume * Ndat; k++) {
        CMUL(scale, D[k][i], tc);
        CDIF(D[k][j], tc);    // D[k][j] -= scale * D[k][i]

        CMUL(scale, D[i][k], tc);
        CDIF(D[j][k], tc);    // D[j][k] -= scale * D[i][k]
      }
    }

#ifdef DEBUG_CHECK
    // Make sure every column of D is still non-vanishing
    for (k = 0; k < volume * Ndat; k++) {
      tr = 0.0;
      for (j = k + 1; j < volume * Ndat; j++)
        tr += cabs_sq(&(D[k][j]));
      if (tr < PFA_TOL) {
        printf("Column %d vanishes after round %d: %.4g\n", k, i, tr);
        break;
      }
    }
#endif
  }

//#ifdef DEBUG_CHECK
  // Check that D is still anti-symmetric, and now tridiagonal
  for (i = 0; i < volume * Ndat; i++) {
    for (j = i + 1; j < volume * Ndat; j++) {
      if (D[i][j].real + D[j][i].real > PFA_TOL
       || D[i][j].imag + D[j][i].imag > PFA_TOL) {
        printf("d_phase: D[%d][%d] = (%.4g, %.4g) but ",
               i, j, D[i][j].real, D[i][j].imag);
        printf("D[%d][%d] = (%.4g, %.4g)\n",
               j, i, D[j][i].real, D[j][i].imag);
      }
    }
  }
//#endif

  // Now we can just compute the pfaffian from D[i][i + 1],
  // accumulating the phase and (the log of) the magnitude
  phase = 0.0;
  log_mag = 0.0;
  for (i = 0; i < volume * Ndat; i += 2) {
//    printf("(%.4g, %.4g)\n", D[i][i + 1].real, D[i][i + 1].imag);
    phase += carg(&(D[i][i + 1]));
    log_mag += log(cabs_sq(&(D[i][i + 1]))) / 2.0;
    free(D[i]);
  }
  free(D);

  // Accumulate phase and (log of) magnitude across all nodes
  g_doublesum(&phase);
  g_doublesum(&log_mag);

  // Account for potential negative sign from interchanges, and mod out 2pi
  if (interchange > 1)
    tr = phase;
  else
    tr = phase + PI;
  phase = fmod(tr, 2.0 * PI);

  node0_printf("PFAFF %.8g %.8g %.8g %.8g\n", log_mag, phase,
               fabs(cos(phase)), fabs(sin(phase)));
}
// -----------------------------------------------------------------
