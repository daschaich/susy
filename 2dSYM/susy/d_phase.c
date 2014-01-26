// -----------------------------------------------------------------
// Measure the phase of the pfaffian
// Based on the algorithm in the appendix of hep-lat/0305002
// hep-lat/9903014 provides important definitions and caveats (!!!)
#include "susy_includes.h"
#define PFA_TOL 1e-8
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return a_i * b_i for two complex vectors (no conjugation!)
double_complex dot(complex *a, complex *b) {
  int i, Ndat = 4 * DIMF;
  complex tc;
  double_complex dot = cmplx(0.0, 0.0);

  for (i = 0; i < sites_on_node * Ndat; i++) {
    CMUL(a[i], b[i], tc);
    CSUM(dot, tc);
  }
  // Accumulate inner product across all nodes
  g_dcomplexsum(&dot);
  return dot;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wrapper for the fermion_op matvec
// Convert between complex vectors and Twist_Fermions
void matvec(complex *in, complex *out) {
  register site *s;
  int i, j, iter;

  // Copy complex vector into Twist_Fermion src
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  // !!! Need to cycle over fields to ensure non-zero Q[i + 1] M Q[i]
  iter = 0;
  FORALLSITES(i, s) {
    for (j = 0; j < DIMF; j++) {
      set_complex_equal(&(in[iter]), &(src[i].Fsite.c[j]));
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Flink[0].c[j]));
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Flink[1].c[j]));
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Fplaq[0][1].c[j]));
      CNEGATE(src[i].Fplaq[0][1].c[j], src[i].Fplaq[1][0].c[j]);
      iter++;
      // Clear diagonal Fplaq
      src[i].Fplaq[0][0].c[j] = cmplx(0.0, 0.0);
      src[i].Fplaq[1][1].c[j] = cmplx(0.0, 0.0);
    }
  }

#ifdef DEBUG_CHECK
  // Check that we didn't miss any components of the input vector
  int Ndat = 4 * DIMF;
  if (iter != sites_on_node * Ndat) {
    printf("d_phase: cycled over %d of %d input components\n",
           iter, sites_on_node * Ndat);
    terminate(1);
  }
#endif

  fermion_op(src, res, 1);    // D
//  fermion_op(src, res, -1);    // Ddag
  Nmatvecs++;

  // Copy the resulting Twist_Fermion res back to complex vector y
  iter = 0;
  FORALLSITES(i, s) {
    for (j = 0; j < DIMF; j++) {
      set_complex_equal(&(res[i].Fsite.c[j]), &(out[iter]));
      iter++;
      set_complex_equal(&(res[i].Flink[0].c[j]), &(out[iter]));
      iter++;
      set_complex_equal(&(res[i].Flink[1].c[j]), &(out[iter]));
      iter++;
      set_complex_equal(&(res[i].Fplaq[0][1].c[j]), &(out[iter]));
      iter++;
    }
  }

#ifdef DEBUG_CHECK
  // Check that we didn't miss any components of the output vector
  if (iter != sites_on_node * Ndat) {
    printf("d_phase: cycled over %d of %d output components\n",
           iter, sites_on_node * Ndat);
    terminate(1);
  }
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void d_phase() {
  register int i, j, k;
  int Ndat = 4 * DIMF, shift = this_node * sites_on_node * Ndat;
  double phase, log_mag, tr, dtime;
  complex temp, temp2;
  complex *diag = malloc(sites_on_node * Ndat * sizeof(*diag));
  complex *MonC = malloc(sites_on_node * Ndat * sizeof(*MonC));
  complex **Q = malloc(volume * Ndat * sizeof(complex*));

  if (Q == NULL) {
    printf("d_phase: can't malloc Q\n");
    fflush(stdout);
    terminate(1);
  }

  // Allocate and initialize Q to unit matrix
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  // We keep part of every column on this node
  for (i = 0; i < volume * Ndat; i++) {
    Q[i] = malloc(sites_on_node * Ndat * sizeof(complex));
    for (j = 0; j < sites_on_node * Ndat; j++)
      Q[i][j] = cmplx(0.0, 0.0);

    // Somewhat hakcy: below we will only set diag[i] on the appropriate node
    // and then scatter it by summing over nodes
    diag[i] = cmplx(0.0, 0.0);
  }
  if (Q[volume * Ndat - 1] == NULL) {
    printf("d_phase: can't malloc Q[i]\n");
    fflush(stdout);
    terminate(1);
  }
  // The diagonal elements are distributed across different nodes
  // according to shift = this_node * sites_on_node * Ndat defined above
  for (i = 0; i < sites_on_node * Ndat; i++)
    Q[shift + i][i] = cmplx(1.0, 0.0);

  // Cycle over ALL pairs of columns
  for (i = 0; i < volume * Ndat - 1; i += 2) {
    node0_printf("Column %d of %d: ", i, volume * Ndat - 2);
    Nmatvecs = 0;
    dtime = -dclock();

    // q_{i + 1} --> q_{i + 1} / <q_{i + 1} | D | q_i>
    // But only if <q_{i + 1} | D | q_i> is non-zero!
    matvec(Q[i], MonC);
    temp = dot(Q[i + 1], MonC);
    // !!! Must have non-vanishing matrix element!
    if (cabs_sq(&temp) < PFA_TOL) {
      printf("d_phase: <i+1 | D | i> = (%.4g, %.4g) too small for i = %d\n",
             temp.real, temp.imag, i);
      terminate(1);
    }

    // All loops over j (the elements of the given column) run in parallel
    // We can shorten many of them since Q is upper-triangular
    if (i + 1 < sites_on_node * Ndat) {
      for (j = 0; j < i + 2; j++) {
        CDIV(Q[i + 1][j], temp, temp2);
        set_complex_equal(&temp2, &(Q[i + 1][j]));
      }
    }
    else {
      for (j = 0; j < sites_on_node * Ndat; j++) {
        CDIV(Q[i + 1][j], temp, temp2);
        set_complex_equal(&temp2, &(Q[i + 1][j]));
      }
    }

    // Cycle over ALL subsequent columns
    for (k = i + 2; k < volume * Ndat; k++) {
      // q_k --> q_k - q_i <q_{i + 1} | D | q_k> - q_{i + 1} <q_i | D | q_k>
      matvec(Q[k], MonC);
      temp = dot(Q[i + 1], MonC);
      if (i + 1 < sites_on_node * Ndat) {
        for (j = 0; j < i + 2; j++) {
          CMUL(Q[i][j], temp, temp2);
          CDIF(Q[k][j], temp2);
        }
      }
      else {
        for (j = 0; j < sites_on_node * Ndat; j++) {
          CMUL(Q[i][j], temp, temp2);
          CDIF(Q[k][j], temp2);
        }
      }

      temp = dot(Q[i], MonC);
      if (i + 1 < sites_on_node * Ndat) {
        for (j = 0; j < i + 2; j++) {
          CMUL(Q[i + 1][j], temp, temp2);
          CSUM(Q[k][j], temp2);
        }
      }
      else {
        for (j = 0; j < sites_on_node * Ndat; j++) {
          CMUL(Q[i + 1][j], temp, temp2);
          CSUM(Q[k][j], temp2);
        }
      }
    }
    // Print some timing information
    // and make sure diagonal elements are still sane
    dtime += dclock();
    node0_printf("%d matvecs in %.4g seconds ", Nmatvecs, dtime);

    // Save diagonal element and free columns i and i+1
    // We use shift to single out the appropriate node, then sum
    // Note that Q[i][i - shift] = 1
    if (0 <= i + 1 - shift && i + 1 - shift < sites_on_node * Ndat)
      set_complex_equal(&(Q[i + 1][i + 1 - shift]), &(diag[i + 1]));
    g_dcomplexsum(&(diag[i + 1]));
    node0_printf("%.8g %.8g\n", diag[i + 1].real, diag[i + 1].imag);
    fflush(stdout);

    // Done with these
    free(Q[i]);
    free(Q[i + 1]);
  }
  free(MonC);
  free(Q);

  // Q is triangular by construction
  // Compute its determinant from its non-trivial diagonal elements diag[i + 1]
  // expressed as the phase and (the log of) the magnitude
  // The diagonal elements are distributed across different nodes
  // according to shift = this_node * sites_on_node * Ndat defined above
  phase = 0.0;
  log_mag = 0.0;
  for (i = 1; i < sites_on_node * Ndat; i += 2) {   // Start at 1, use diag[i]
#ifdef DEBUG_CHECK
    // Check: print out all diagonal elements
    printf("%d+%d: (%.4g, %.4g) --> exp[%.4g + %.4gi]\n",
           this_node, i, diag[i].real, diag[i].imag,
           log(cabs_sq(&(diag[i]))) / 2.0, carg(&(diag[i])));
#endif

    phase += carg(&(diag[i]));
    log_mag += log(cabs_sq(&(diag[i]))) / 2.0;
  }

  // Accumulate phase and (log of) magnitude across all nodes
  g_doublesum(&phase);
  g_doublesum(&log_mag);

  // pf(M) = (det Q)^{-1}
  // Negate log of magnitude and phase
  // Following the C++ code, keep phase in [0, 2pi)
  tr = fmod(-1.0 * phase, TWOPI);
  if (tr < 0)
    phase = tr + TWOPI;
  else
    phase = tr;
  node0_printf("PFAFF %.8g %.8g %.8g %.8g\n", -1.0 * log_mag, phase,
               fabs(cos(phase)), fabs(sin(phase)));
}
// -----------------------------------------------------------------
