// -----------------------------------------------------------------
// Measure Tr[X^2] / N, where X is the scalar field
// from the log of hermitian part of the polar decomposition
#include "susy_includes.h"

void measure_tr_xsq() {
  register int i, dir;
  register site *s;
  double norm = 1.0 / ((double)(NUMLINK * volume * NCOL));
  complex tc, tr_xsq = cmplx(0.0, 0.0), tr_x_comm = cmplx(0.0, 0.0);
  matrix X[NDIMS], tmat;

  FORALLSITES(i, s) {
    FORALLDIR(dir) {
      polar(&(s->link[dir]), &(X[dir]), &tmat);
      matrix_log(&tmat, &(X[dir]));

      // Note that we sum over both links
      // and don't normalize by 2 below...
      tc = complextrace_nn(&(X[dir]), &(X[dir]));
      CSUM(tr_xsq, tc);
    }
    mult_nn(&(X[0]), &(X[1]), &tmat);
    mult_nn_dif(&(X[1]), &(X[0]), &tmat);
    tc = complextrace_nn(&tmat, &tmat);
    CSUM(tr_x_comm, tc);
  }
  CMULREAL(tr_xsq, norm, tr_xsq);
  g_complexsum(&tr_xsq);
  g_complexsum(&tr_x_comm);

  if (fabs(tr_xsq.imag) > IMAG_TOL)
    node0_printf("WARNING Im(tr_xsq) = %.4g\n", tr_xsq.imag);
  if (fabs(tr_x_comm.imag) > IMAG_TOL)
    node0_printf("WARNING Im(tr_x_comm) = %.4g\n", tr_x_comm.imag);

  node0_printf("TR_XSQ %.8g\n", tr_xsq.real);
  node0_printf("TR_X_COMM %.8g\n", tr_x_comm.real);
}
// -----------------------------------------------------------------
