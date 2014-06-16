// -----------------------------------------------------------------
// Measure the Konishi and SUGRA scalar field correlation functions
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA correlators, projected to zero spatial momentum
void d_correlator() {
  register int i;
  register site *s;
  int a, b, t, my_t, tt;
  Real norm, dum = 0.0, corr;
  Real *tK, *tS[NUMLINK][NUMLINK];
  complex ctmp;
  su3_matrix_f *B[NUMLINK], tmat;

  // Allocate tK and tS_{ab}
  tK = malloc(nt * sizeof(*tK));
  for (a = 0; a < NUMLINK; a++) {
    for (b = 0; b < NUMLINK; b++)
      tS[a][b] = malloc(nt * sizeof(*tS[a][b]));
  }

  // Initialize tK and tS_{ab} to zero
  for (t = 0; t < nt; t++) {
    tK[t] = 0.0;
    for (a = 0; a < NUMLINK; a++) {
      for (b = 0; b < NUMLINK; b++)
        tS[a][b][t] = 0.0;
    }
  }

  // Finally allocate B_a and start calculating
  for (a = 0; a < NUMLINK; a++)
    B[a] = malloc(sites_on_node * sizeof(su3_matrix_f));

  FORALLSITES(i, s) {
    for (a = 0; a < NUMLINK; a++)
      dum += realtrace_su3_f(&(s->linkf[a]), &(s->linkf[a]));
  }       // realtrace_su3_f(A, B) returns ReTr(Adag B)
  g_doublesum(&dum);
  dum /= (Real)(NCOL * NUMLINK * volume);

  // B_a = U_a Udag_a - Sum_{a, x} ReTr(U_a Udag_a) / (5Nc vol)
  FORALLSITES(i, s) {
    for (a = 0; a < NUMLINK; a++) {
      mult_su3_na_f(&(s->linkf[a]), &(s->linkf[a]), &(B[a][i]));
      for (t = 0; t < NCOL; t++)
        B[a][i].e[t][t].real -= dum;    // Subtract volume average
    }
  }

  // Now form the correlators, first summing over each timeslice
  // To get average, would normalize correlators by (Nt / vol)^2 below
  FORALLSITES(i, s) {
    my_t = s->t;
    for (a = 0; a < NUMLINK; a++) {
      for (b = 0; b < NUMLINK; b++) {
        mult_su3_nn_f(&(B[a][i]), &(B[b][i]), &tmat);
        ctmp = trace_su3_f(&tmat);

        // Make sure Tr(B_a * B_b) really is real
        if (abs(ctmp.imag) > IMAG_TOL) {
          printf("node%d WARNING: Im(X[%d][%d][%d]) = %.4g > %.4g)\n",
                 this_node, a, b, i, ctmp.imag, IMAG_TOL);
        }

        tS[a][b][my_t] += ctmp.real;
        if (a == b)
          tK[my_t] += ctmp.real;
      }
    }
  }
  // Done with B[a]
  for (a = 0; a < NUMLINK; a++)
    free(B[a]);

  // Sum tK and tS_{ab} across nodes before forming correlators
  for (t = 0; t < nt; t++) {
    g_doublesum(&tK[t]);
    for (a = 0; a < NUMLINK; a++) {
      for (b = 0; b < NUMLINK; b++)
        g_doublesum(&tS[a][b][t]);
    }
  }

  // Subtract trace from SUGRA operator
  for (t = 0; t < nt; t++) {
    for (a = 0; a < NUMLINK; a++)
      tS[a][a][t] -= tK[t] / (Real)NUMLINK;
  }

  // Form and print out correlators, normalized by Nt / vol^2
  // (To average over tt, remove one factor of Nt from normalization)
  // Konishi
  norm = (Real)(nx * ny * nz * volume);
  for (t = 0; t <= (int)(nt / 2); t++) {
    corr = tK[0] * tK[t];
    for (tt = 1; tt < nt; tt++)
      corr += tK[tt] * tK[(t + tt) % nt];
    corr /= norm;
    node0_printf("KONISHI %d %.8g\n", t, corr);
  }

  // SUGRA -- print out all 25
  for (a = 0; a < NUMLINK; a++) {
    for (b = 0; b < NUMLINK; b++) {
      for (t = 0; t <= (int)(nt / 2); t++) {
        corr = tS[a][b][0] * tS[a][b][t];
        for (tt = 1; tt < nt; tt++)
          corr += tS[a][b][tt] * tS[a][b][(t + tt) % nt];
        corr /= norm;
        node0_printf("SUGRA %d %d %d %.8g\n", a, b, t, corr);
      }
    }
  }

  free(tK);
  for (a = 0; a < NUMLINK; a++) {
    for (b = 0; b < NUMLINK; b++)
      free(tS[a][b]);
  }
}
// -----------------------------------------------------------------
