// -----------------------------------------------------------------
// Measure the average value of Tr[X[i] X[i]] / N
// as well as the width sqrt(<tr^2> - <tr>^2) of its distribution
#include "susy_includes.h"

double scalar_trace(double *Xtr, double *Xwidth) {
  
  register int i, j;
  register site *s;
  double Xtr_ave = 0.0, XtrSq = 0.0, td;
  
  FORALLSITES(i, s) {
    for(j=0;j<NSCALAR;j++) {
      td = realtrace(&(s->X[j]), &(s->X[j]));
      Xtr[j] += td;
      XtrSq += td * td;
    }
  }
  
  for(j=0;j<NSCALAR;j++) {
    Xtr[j] *= (one_ov_N / (double) nt);
    g_doublesum(&(Xtr[j]));
    Xtr_ave += Xtr[j];
  }

  Xtr_ave /= (double) NSCALAR;
  XtrSq /= ((double) nt * NSCALAR);
  g_doublesum(&XtrSq);
  td = Xtr_ave;
  *Xwidth = sqrt(XtrSq - td * td);

  return Xtr_ave;
}
// -----------------------------------------------------------------
