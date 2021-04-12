// -----------------------------------------------------------------
// Measure Tr[X^2] / N, where X is the scalar field
// from the log of hermitian part of the polar decomposition
#include "susy_includes.h"

void measure_tr_xsq() {
  register int i, dir;
  register site *s;
  double norm = 1.0 / ((double)(NUMLINK * volume * NCOL));
  complex tc, tr_xsq = cmplx(0.0, 0.0);
  matrix X, tmat;

  FORALLSITES(i, s) {
    FORALLDIR(dir) {
      polar(&(s->link[dir]), &X, &tmat);
      matrix_log(&tmat, &X);

      // Note that we sum over both links
      // and don't normalize by 2 below...
      tc = complextrace_nn(&X, &X);
      CSUM(tr_xsq, tc);
    }
  }
  CMULREAL(tr_xsq, norm, tr_xsq);
  g_complexsum(&tr_xsq);

  if (fabs(tr_xsq.imag) > IMAG_TOL)
    node0_printf("WARNING Im(tr_xsq) = %.4g\n", tr_xsq.imag);

  node0_printf("TR_XSQ %.8g\n", tr_xsq.real);
}
// -----------------------------------------------------------------
