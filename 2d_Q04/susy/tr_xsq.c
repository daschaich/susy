// -----------------------------------------------------------------
// Measure Tr[X^2] / N, where X is the scalar field
// from the polar decomposition
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
      // Take log
      matrix_log(&tmat, &X);
      
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
