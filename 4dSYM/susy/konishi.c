// -----------------------------------------------------------------
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute traces of three bilinears using two scalar field interpolating ops
// Op 0 ("P") is log of hermitian matrix P from polar decomposition
// Op 1 ("U") is traceless part of U_a Udag_a
// Bilin 0 is Tr[PP]
// Bilin 1 is (Tr[UP] + Tr[PU]) / 2
// Bilin 2 is Tr[UU]
void compute_Ba() {
  register int i;
  register site *s;
  int a, b, j, k;
  Real tr;
  complex tc;
  matrix_f tmat;

  FORALLSITES(i, s) {
    // Construct scalar fields
    for (a = XUP; a < NUMLINK; a++) {
      // Log of hermitian part of polar decomposition
      polar(&(s->linkf[a]), &(Ba[0][a][i]), &tmat);
      matrix_log(&tmat, &(Ba[0][a][i]));

      // Traceless part of U.Udag (hermitian so trace is real)
      mult_na_f(&(s->linkf[a]), &(s->linkf[a]), &(Ba[1][a][i]));
      tc = trace_f(&(Ba[1][a][i]));
      tr = tc.real / (Real)NCOL;
      for (k = 0; k < NCOL; k++)
        Ba[1][a][i].e[k][k].real -= tr;

      // Check reality of resulting scalar fields
      for (j = 0; j < N_B; j++){
        tc = trace_f(&(Ba[j][a][i]));
        if (fabs(tc.imag) > IMAG_TOL) {
          printf("WARNING: Tr[Ba[%d][%d][%d]] = (%.4g, %.4g) is not real\n",
                 j, a, i, tc.real, tc.imag);
        }
      }
    }
  }

  // Compute traces of bilinears
  // Symmetric in a <--> b but store all to simplify SUGRA computation
  // Have checked that all are purely real and gauge invariant
  for (a = XUP; a < NUMLINK; a++) {
    for (b = XUP; b < NUMLINK; b++) {
      FORALLSITES(i, s) {
        mult_nn_f(&(Ba[0][a][i]), &(Ba[0][b][i]), &tmat);
        tc = trace_f(&tmat);
        traceBB[0][a][b][i] = tc.real;    // Tr[PP]

        mult_nn_f(&(Ba[1][a][i]), &(Ba[1][b][i]), &tmat);
        tc = trace_f(&tmat);
        traceBB[2][a][b][i] = tc.real;    // Tr[UU]

        // (Tr[UP] + Tr[PU]) / 2
        mult_nn_f(&(Ba[0][a][i]), &(Ba[1][b][i]), &tmat);
        tc = trace_f(&tmat);
        traceBB[1][a][b][i] = 0.5 * tc.real;
        mult_nn_f(&(Ba[1][a][i]), &(Ba[0][b][i]), &tmat);
        tc = trace_f(&tmat);
        traceBB[1][a][b][i] += 0.5 * tc.real;
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure and print the Konishi and SUGRA operators on each timeslice
void konishi() {
  register int i;
  register site *s;
  int a, b, t, j;
  double norm, *OK[N_K], *OS[N_K];

  // Allocate and initialize Konishi and SUGRA operators on each time slice
  for (j = 0; j < N_K; j++) {
    OK[j] = malloc(nt * sizeof(double));
    OS[j] = malloc(nt * sizeof(double));
    for (t = 0; t < nt; t++) {
      OK[j][t] = 0.0;
      OS[j][t] = 0.0;
    }
  }

  // Compute traces of bilinears of scalar field interpolating ops
  compute_Ba();

  // Now form the zero momentum projected operators (summing across nodes)
  FORALLSITES(i, s) {
    t = s->t;
    for (a = XUP; a < NUMLINK; a++) {
      for (j = 0; j < N_K; j++) {
        OK[j][t] += traceBB[j][a][a][i];    // Konishi

        // Now SUGRA, averaged over 20 off-diagonal components
        for (b = XUP; b < NUMLINK; b++) {
          if (a == b)
            continue;
          OS[j][t] += 0.05 * traceBB[j][a][b][i];
        }
      }
    }
  }
  // Normalization removed from site loop, followed by sum over nodes
  norm = (Real)(nx * ny * nz);
  for (t = 0; t < nt; t++) {
    for (j = 0; j < N_K; j++) {
      OK[j][t] /= norm;
      g_doublesum(&(OK[j][t]));

      OS[j][t] /= norm;
      g_doublesum(&(OS[j][t]));
    }
  }

  // Print each operator on each time slice
  // Subtract either ensemble average or volume average (respectively)
  // Format: TAG  t  a  op[a]  subtracted[a]
  for (j = 0; j < N_K; j++) {
    volK[j] = OK[j][0];
    volS[j] = OS[j][0];
    for (t = 1; t < nt; t++) {
      volK[j] += OK[j][t];
      volS[j] += OS[j][t];
    }
    volK[j] /= (double)nt;
    volS[j] /= (double)nt;
  }

  for (t = 0; t < nt; t++) {
    for (j = 0; j < N_K; j++)
      node0_printf("KONISHI %d %d %.16g %.16g\n", t, j,
                   OK[j][t] - vevK[j], OK[j][t] - volK[j]);
  }
  for (t = 0; t < nt; t++) {
    for (j = 0; j < N_K; j++)     // Assume vanishing vev
      node0_printf("SUGRA %d %d %.16g %.16g\n", t, j,
                   OS[j][t], OS[j][t] - volS[j]);
  }

  for (j = 0; j < N_K; j++) {
    free(OK[j]);
    free(OS[j]);
  }
}
// -----------------------------------------------------------------
