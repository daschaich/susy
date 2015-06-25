// -----------------------------------------------------------------
// Helper functions for smearing and measurements
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate newU = exp(Q).U, overwriting s->linkf
// Here Q is the traceless anti-hermitian lattice field from stout_smear
// Go to eighth order in the exponential:
//   exp(Q) * U = (1 + Q + Q^2/2 + Q^3/6 ...) * U
//              = U + Q*(U + (Q/2)*(U + (Q/3)*( ... )))
void exp_mult() {
  register int i, dir;
  register site *s;
  register Real t2, t3, t4, t5, t6, t7, t8;
  su3_matrix_f *link, tmat, tmat2, htemp;

  // Take divisions out of site loop (can't be done by compiler)
  t2 = 1.0 / 2.0;
  t3 = 1.0 / 3.0;
  t4 = 1.0 / 4.0;
  t5 = 1.0 / 5.0;
  t6 = 1.0 / 6.0;
  t7 = 1.0 / 7.0;
  t8 = 1.0 / 8.0;

  for (dir = XUP; dir < NUMLINK; dir++) {
    FORALLSITES(i, s) {
      uncompress_anti_hermitian(&(Q[dir][i]), &htemp);
      link = &(s->linkf[dir]);

      mult_su3_nn_f(&htemp, link, &tmat);
      scalar_mult_add_su3_matrix_f(link, &tmat, t8, &tmat2);

      mult_su3_nn_f(&htemp, &tmat2, &tmat);
      scalar_mult_add_su3_matrix_f(link, &tmat, t7, &tmat2);

      mult_su3_nn_f(&htemp, &tmat2, &tmat);
      scalar_mult_add_su3_matrix_f(link, &tmat, t6, &tmat2);

      mult_su3_nn_f(&htemp, &tmat2, &tmat);
      scalar_mult_add_su3_matrix_f(link, &tmat, t5, &tmat2);

      mult_su3_nn_f(&htemp, &tmat2, &tmat);
      scalar_mult_add_su3_matrix_f(link, &tmat, t4, &tmat2);

      mult_su3_nn_f(&htemp, &tmat2, &tmat);
      scalar_mult_add_su3_matrix_f(link, &tmat, t3, &tmat2);

      mult_su3_nn_f(&htemp, &tmat2, &tmat);
      scalar_mult_add_su3_matrix_f(link, &tmat, t2, &tmat2);

      mult_su3_nn_f(&htemp, &tmat2, &tmat);
      add_su3_matrix_f(link, &tmat, &(s->linkf[dir]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute traces of bilinears of scalar field interpolating ops:
//   traceless part of the hermitian matrix returned by the polar projection
//   traceless part of U_a Udag_a
#ifdef CORR
void compute_Ba() {
  register int i;
  register site *s;
  int a, b, j;
  complex tc;
  Real tr;
  su3_matrix_f tmat;

  // B_a is traceless part of the hermitian matrix from the polar projection
  // C_a is traceless part of U_a Udag_a
  FORALLSITES(i, s) {
    for (a = XUP; a < NUMLINK; a++) {
      polar(&(s->linkf[a]), &tmat, &(Ba[a][i]));
      tc = trace_su3_f(&(Ba[a][i]));
      CDIVREAL(tc, (Real)NCOL, tc);
      for (j = 0; j < NCOL; j++)
        CDIF(Ba[a][i].e[j][j], tc);

      mult_su3_na_f(&(s->linkf[a]), &(s->linkf[a]), &(Ca[a][i]));
      tc = trace_su3_f(&(Ca[a][i]));
      CDIVREAL(tc, (Real)NCOL, tc);
      for (j = 0; j < NCOL; j++)
        CDIF(Ca[a][i].e[j][j], tc);
    }

    // Bilinear traces are symmetric in a <--> b,
    // but store all to simplify SUGRA computation
    // Make sure all are purely real
    for (a = XUP; a < NUMLINK; a++) {
      for (b = a; b < NUMLINK; b++) {
        mult_su3_nn_f(&(Ba[a][i]), &(Ba[b][i]), &tmat);
        tc = trace_su3_f(&tmat);
        traceBB[0][a][b][i] = tc.real;
        traceBB[0][a][b][i] = traceBB[0][b][a][i];
        if (fabs(tc.imag) > IMAG_TOL) {
          printf("WARNING: Tr(BB[0][%d][%d]) = (%.4g, %.4g) at site %d\n",
                 a, b, tc.real, tc.imag, i);
        }


        mult_su3_nn_f(&(Ba[a][i]), &(Ca[b][i]), &tmat);
        tc = trace_su3_f(&tmat);
        traceBB[1][a][b][i] = tc.real;
        traceBB[1][a][b][i] = traceBB[1][b][a][i];
        if (fabs(tc.imag) > IMAG_TOL) {
          printf("WARNING: Tr(BB[1][%d][%d]) = (%.4g, %.4g) at site %d\n",
                 a, b, tc.real, tc.imag, i);
        }

        // Will probably be unchanged by stout smearing
        mult_su3_nn_f(&(Ca[a][i]), &(Ca[b][i]), &tmat);
        tc = trace_su3_f(&tmat);
        traceBB[2][a][b][i] = tc.real;
        traceBB[2][a][b][i] = traceBB[2][b][a][i];
        if (fabs(tc.imag) > IMAG_TOL) {
          printf("WARNING: Tr(BB[2][%d][%d]) = (%.4g, %.4g) at site %d\n",
                 a, b, tc.real, tc.imag, i);
        }
      }
    }
  }
}
#endif
// -----------------------------------------------------------------
