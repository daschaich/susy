// -----------------------------------------------------------------
// Measure Tr[X^2] / N, where X is the scalar field
// from the log of hermitian part of the polar decomposition
#include "susy_includes.h"

void measure_tr_xsq() {
  register int i, dir;
  register site *s;
  double norm = 1.0 / ((double)(NUMLINK * volume * NCOL));
  complex tc, tc_comm, tr_xsq = cmplx(0.0, 0.0), tr_x_comm = cmplx(0.0, 0.0);
  matrix X, tmat, tmat2, tmat3,  tmat4, X_one, X_two;

  FORALLSITES(i, s) {
    FORALLDIR(dir) {
      polar(&(s->link[dir]), &X, &tmat);
      matrix_log(&tmat, &X);
      if(dir==TUP){
      mat_copy(&X,&X_one);
      }
      if(dir==XUP){
      mat_copy(&X,&X_two);
      }

      // Note that we sum over both links
      // and don't normalize by 2 below...
      tc = complextrace_nn(&X, &X);
      CSUM(tr_xsq, tc);
    }
   //mult_nn_diff(&X_one,&X_two,&tmat4);
   mult_nn(&X_one,&X_two,&tmat2);
   mult_nn(&X_two,&X_one,&tmat3);
   sub_matrix(&tmat2,&tmat3,&tmat4);
   tc_comm = complextrace_nn(&tmat4,&tmat4);
   CSUM(tr_x_comm,tc_comm);
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

