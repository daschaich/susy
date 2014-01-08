// -----------------------------------------------------------------
// Measure the phase of the pfaffian
// Based on the algorithm in the appendix of hep-lat/0305002
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the dot product of two complex vectors, adag b
double_complex dot(complex *a, complex *b) {
  int i, Ndat = 16 * DIMF;
  complex tc;
  double_complex dot = cmplx(0.0, 0.0);

  for (i = 0; i < sites_on_node * Ndat; i++) {
    CMUL(a[i], b[i], tc);   // CMULJ_ leads to divergences
    CSUM(dot, tc);
  }
  return dot;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wrapper for the fermion_op matvec
// Convert between complex vectors and Twist_Fermions
void matvec(complex *in, complex *out) {
  register site *s;
  int i, j, mu, nu, iter, Ndat = 16 * DIMF;
  complex *xx;
  Twist_Fermion src[sites_on_node], res[sites_on_node];

  // Copy complex vector into Twist_Fermion src
  // Each Twist_Fermion has Ndat = 16DIMF non-trivial complex components
  FORALLSITES(i, s) {
    clear_TF(&(src[i]));    // Clear diagonal plaquette fermions
    xx = in + Ndat * i;     // This particular site in the vector
    iter = 0;
    for (j = 0; j < DIMF; j++) {
      set_complex_equal(&(xx[iter]), &(src[i].Fsite.c[j]));
      iter++;
    }
    for (mu = 0; mu < NUMLINK; mu++) {   // Do all links before any plaqs
      for (j = 0; j < DIMF; j++) {
        set_complex_equal(&(xx[iter]), &(src[i].Flink[mu].c[j]));
        iter++;
      }
    }
    for (mu = 0; mu < NUMLINK; mu++) {
      for (nu = mu + 1; nu < NUMLINK; nu++) {
        for (j = 0; j < DIMF; j++) {
          set_complex_equal(&(xx[iter]), &(src[i].Fplaq[mu][nu].c[j]));
          CNEGATE(src[i].Fplaq[mu][nu].c[j], src[i].Fplaq[nu][mu].c[j]);
          iter++;
        }
      }
    }
  }

  fermion_op(src, res, 1);    // D

  // Copy the resulting Twist_Fermion res back to complex vector y
  // Each Twist_Fermion has Ndat = 16DIMF non-trivial complex components
  FORALLSITES(i, s) {
    xx = out + Ndat * i;    // This particular site in the vector
    iter = 0;
    for (j = 0; j < DIMF; j++) {
      set_complex_equal(&(res[i].Fsite.c[j]), &(xx[iter]));
      iter++;
    }
    for (mu = 0; mu < NUMLINK; mu++) {   // Do all links before any plaqs
      for (j = 0; j < DIMF; j++) {
        set_complex_equal(&(res[i].Flink[mu].c[j]), &(xx[iter]));
        iter++;
      }
    }
    for (mu = 0; mu < NUMLINK; mu++) {
      for (nu = mu + 1; nu < NUMLINK; nu++) {
        for (j = 0; j < DIMF; j++) {
          set_complex_equal(&(res[i].Fplaq[mu][nu].c[j]), &(xx[iter]));
          iter++;
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void d_phase() {
  register int i, j, k;
  int Ndat = 16 * DIMF;
  double phase, log_mag, tr;
  complex temp, temp2;
  complex *MonC = malloc(sites_on_node * Ndat * sizeof(*MonC));
  complex **Q = malloc(sites_on_node * Ndat * sizeof(complex*));

  // Allocate and initialize Q to unit matrix
  // Each Twist_Fermion has Ndat = 16DIMF non-trivial complex components
  for (i = 0; i < sites_on_node * Ndat; i++) {
    Q[i] = malloc(sites_on_node * Ndat * sizeof(complex));
    for (j = 0; j < sites_on_node * Ndat; j++)
      Q[i][j] = cmplx(0.0, 0.0);
    Q[i][i] = cmplx(1.0, 0.0);
  }

  // Cycle over pairs of columns
  for (i = 0; i < sites_on_node * Ndat - 1; i += 2) {
    // e_{i + 1} --> e_{i + 1} / <e_{i + 1} | D | e_i>
    // But only if <e_{i + 1} | D | e_i> is non-zero!
    matvec(Q[i], MonC);
    temp = dot(Q[i + 1], MonC);
    if (cabs_sq(&temp) > IMAG_TOL) {
      // Check complex number by which we're dividing
//      printf("(%.4g, %.4g)\n", temp.real, temp.imag);
      for (j = 0; j < sites_on_node * Ndat - 1; j++) {
        CDIV(Q[i + 1][j], temp, temp2);
        set_complex_equal(&temp2, &(Q[i + 1][j]));
      }
    }

    // Cycle over all subsequent columns
    for (k = i + 2; k < sites_on_node * Ndat - 1; k++) {
      // e_k --> e_k - e_i <e_{i + 1} | D | e_k> - e_{i + 1} <e_i | D | e_k>
      matvec(Q[k], MonC);
      temp = dot(Q[i + 1], MonC);
      for (j = 0; j < sites_on_node * Ndat - 1; j++) {
        CMUL(Q[i][j], temp, temp2);
        CDIF(Q[k][j], temp2);
      }

      temp = dot(Q[i], MonC);
      for (j = 0; j < sites_on_node * Ndat - 1; j++) {
        CMUL(Q[i + 1][j], temp, temp2);
        CDIF(Q[k][j], temp2);
      }
    }
  }

  // Q is triangular by construction
  // Compute its determinant from its diagonal elements
  // expressed as the phase and (the log of) the magnitude
  phase = 0.0;
  log_mag = 0.0;
  for (i = 0; i < sites_on_node * Ndat; i++) {
    // Check: print out all diagonal elements
//    printf("%d: (%.4g, %.4g) --> %.4g exp[%.4g]\n",
//           i, Q[i][i].real, Q[i][i].imag,
//           cabs_sq(&(Q[i][i])), carg(&(Q[i][i])));

    phase += carg(&(Q[i][i]));
    log_mag += log(cabs_sq(&(Q[i][i])));    // Extra factor of 2
    free(Q[i]);
  }
  free(Q);
  free(MonC);

  // Accumulate phase and (log of) magnitude across all nodes
  g_doublesum(&phase);
  g_doublesum(&log_mag);

  // pf(M) = (det Q)^{-1}
  // Invert (log of) magnitude, accounting for factor of 2 above
  // Negate phase and mod out 2pi
  tr = -1.0 * phase;
  phase = fmod(tr, 2.0 * PI);
  node0_printf("PFAFF %.8g %.8g %.8g %.8g\n", 2.0 / log_mag, phase,
               fabs(cos(phase)), fabs(sin(phase)));
}
// -----------------------------------------------------------------
