// -----------------------------------------------------------------
// Measure the average value of the link Tr[Udag U] / N
// as well as the width sqrt(<tr^2> - <tr>^2) of its distribution
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Tr[Udag U] / N
double d_link(double *linktr, double *width) {
  register int i, dir;
  register site *s;
  double ave = 0.0, linktrSq = 0.0, td;

  for (dir = XUP; dir < NUMLINK; dir++) {
    linktr[dir] = 0.0;
    FORALLSITES(i, s) {
      td = realtrace_su3_f(&(s->linkf[dir]), &(s->linkf[dir]));
      linktr[dir] += td;
      linktrSq += td * td;
    }
    linktr[dir] /= ((double)volume * NCOL);
    g_doublesum(&(linktr[dir]));
    ave += linktr[dir];
  }
  ave /= ((double)NUMLINK);
  linktrSq /= ((double)volume * NCOL * NCOL * NUMLINK);
  g_doublesum(&linktrSq);
  *width = sqrt(linktrSq - ave * ave);

  return ave;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Tr[Udag U] / N for links in the fermion irrep
double d_link_frep(double *linktr, double *width) {
  register int i, dir;
  register site *s;
  double ave = 0.0, linktrSq = 0.0, td;

  for (dir = XUP; dir < NUMLINK; dir++) {
    linktr[dir] = 0.0;
    FORALLSITES(i, s) {
      td = realtrace_su3(&(s->link[dir]), &(s->link[dir]));
      linktr[dir] += td;
      linktrSq += td * td;
    }
    linktr[dir] /= ((double)volume * NCOL);
    g_doublesum(&(linktr[dir]));
    ave += linktr[dir];
  }
  ave /= ((double)NUMLINK);
  linktrSq /= ((double)volume * NCOL * NCOL * NUMLINK);
  g_doublesum(&linktrSq);
  *width = sqrt(linktrSq - ave * ave);

  return ave;
}
// -----------------------------------------------------------------
