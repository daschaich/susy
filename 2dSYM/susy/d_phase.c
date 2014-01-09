// -----------------------------------------------------------------
// Measure the phase of the pfaffian
// Based on the algorithm in the appendix of hep-lat/0305002
// hep-lat/9903014 provides important definitions and caveats (!!!)
#include "susy_includes.h"
#define PFA_TOL 1e-8
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the inner product of two complex vectors, (a b)
double_complex dot(complex *a, complex *b) {
  int i, Ndat = 4 * DIMF;
  complex tc;
  double_complex dot = cmplx(0.0, 0.0);

  for (i = 0; i < sites_on_node * Ndat; i++) {
    CMUL(a[i], b[i], tc);   // !!! Note no conjugation!
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
  int i, j, mu, nu, iter, Ndat = 4 * DIMF;

  // Copy complex vector into Twist_Fermion src
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  // !!! Need adjacent site/link/plaq fermions for functioning basis!
  FORALLSITES(i, s) {
    clear_TF(&(src[i]));    // Clear diagonal plaquette fermions
    iter = i * Ndat;
    for (j = 0; j < DIMF; j++) {
      set_complex_equal(&(in[iter]), &(src[i].Fsite.c[j]));
      iter++;
      for (mu = 0; mu < NUMLINK; mu++) {   // Do all links before any plaqs
        set_complex_equal(&(in[iter]), &(src[i].Flink[mu].c[j]));
        iter++;
        for (nu = mu + 1; nu < NUMLINK; nu++) {
          set_complex_equal(&(in[iter]), &(src[i].Fplaq[mu][nu].c[j]));
          CNEGATE(src[i].Fplaq[mu][nu].c[j], src[i].Fplaq[nu][mu].c[j]);
          iter++;
        }
      }
    }
  }

#ifdef DEBUG_CHECK
  // Check that in and src have the same magnitude
  complex tc = dot(in, in);
  printf(" in magsq = %.4g\n", tc.real);

  Real tr = 0.0;
  FORALLSITES(i, s)
    tr += magsq_TF(&(src[i]));
  printf("src magsq = %.4g\n", tr);
#endif

  fermion_op(src, res, 1);    // D

  // Copy the resulting Twist_Fermion res back to complex vector y
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
  // !!! Need adjacent site/link/plaq fermions for functioning basis!
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

#ifdef DEBUG_CHECK
  // Check that in and src have the same magnitude
  tr = 0.0;
  FORALLSITES(i, s)
    tr += magsq_TF(&(res[i]));
  printf(" res magsq = %.4g\n", tr);

  tc = dot(out, out);
  printf(" out magsq = %.4g\n", tc.real);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void d_phase() {
  register int i, j, k;
  int Ndat = 4 * DIMF;
  double phase, log_mag, tr;
  complex temp, temp2;
  complex *MonC = malloc(sites_on_node * Ndat * sizeof(*MonC));
  complex **Q = malloc(sites_on_node * Ndat * sizeof(complex*));

  // Allocate and initialize Q to unit matrix
  // Each Twist_Fermion has Ndat = 4DIMF non-trivial complex components
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
#ifdef DEBUG_CHECK
    // Check magnitude of MonC, compare with magnitude of res above
    temp = dot(MonC, MonC);
    printf("MonC magsq = %.4g\n", temp.real);
#endif

    temp = dot(Q[i + 1], MonC);
    // !!! Must require non-vanishing matrix element!
    if (cabs_sq(&temp) < PFA_TOL) {
      printf("d_phase: <i+1 | D | i> = (%.4g, %.4g) too small for i = %d\n",
             temp.real, temp.imag, i);
      terminate(1);
    }
#ifdef DEBUG_CHECK
    // Check complex number by which we're dividing
    printf("(%.4g, %.4g)\n", temp.real, temp.imag);
#endif
    for (j = 0; j < sites_on_node * Ndat; j++) {
      CDIV(Q[i + 1][j], temp, temp2);
      set_complex_equal(&temp2, &(Q[i + 1][j]));
    }

    // Cycle over all subsequent columns
    for (k = i + 2; k < sites_on_node * Ndat; k++) {
      // e_k --> e_k - e_i <e_{i + 1} | D | e_k> - e_{i + 1} <e_i | D | e_k>
      matvec(Q[k], MonC);
#ifdef DEBUG_CHECK
      // Check magnitude of MonC, compare with magnitude of res above
      temp = dot(MonC, MonC);
      printf("MonC magsq = %.4g\n", temp.real);
#endif

      temp = dot(Q[i + 1], MonC);
      for (j = 0; j < sites_on_node * Ndat; j++) {
        CMUL(Q[i][j], temp, temp2);
        CDIF(Q[k][j], temp2);
      }

      temp = dot(Q[i], MonC);
      for (j = 0; j < sites_on_node * Ndat; j++) {
        CMUL(Q[i + 1][j], temp, temp2);
        CSUM(Q[k][j], temp2);
      }
    }
  }

//#ifdef DEBUG_CHECK
  // Checks on Q from hep-lat/0305002:
  // <e_j | D | e_i> = 0 except
  // <e_{i + 1} | D | e_i> = 1 for i even
  // <e_{i - 1} | D | e_i> = -1 for i odd
  complex one = cmplx(1.0, 0.0);
  for (i = 0; i < sites_on_node * Ndat; i++) {
    matvec(Q[i], MonC);
    for (j = 0; j < sites_on_node * Ndat; j++) {
      temp = dot(Q[j], MonC);
      if (j == i + 1 && i % 2 == 0) {   // Braces suppress compiler warning
        CDIF(temp, one);    // Should vanish
      }
      else if (j == i - 1 && i % 2 == 1)
        CSUM(temp, one);    // Should vanish

      if (temp.real > PFA_TOL || temp.imag > PFA_TOL) {
        temp = dot(Q[j], MonC);
        printf("PROBLEM: <%d | D | %d> = (%.4g, %.4g)\n",
               j, i, temp.real, temp.imag);
      }
    }
  }
//#endif

  // Q is triangular by construction
  // Compute its determinant from its diagonal elements
  // expressed as the phase and (the log of) the magnitude
  phase = 0.0;
  log_mag = 0.0;
  for (i = 0; i < sites_on_node * Ndat; i++) {
#ifdef DEBUG_CHECK
    // Check: print out all diagonal elements
    printf("%d: (%.4g, %.4g) --> exp[%.4g + %.4gi]\n",
           i, Q[i][i].real, Q[i][i].imag,
           log(cabs_sq(&(Q[i][i]))) / 2.0, carg(&(Q[i][i])));
#endif

    phase += carg(&(Q[i][i]));
    log_mag += log(cabs_sq(&(Q[i][i]))) / 2.0;
    free(Q[i]);
  }
  free(Q);
  free(MonC);

  // Accumulate phase and (log of) magnitude across all nodes
  g_doublesum(&phase);
  g_doublesum(&log_mag);

  // pf(M) = (det Q)^{-1}
  // Negate log of magnitude and phase, modding out 2pi
  tr = -1.0 * phase;
  phase = fmod(tr, 2.0 * PI);
  node0_printf("PFAFF %.8g %.8g %.8g %.8g\n", -1.0 * log_mag, phase,
               fabs(cos(phase)), fabs(sin(phase)));
}
// -----------------------------------------------------------------
