// -----------------------------------------------------------------
// Measure the average value of the link Tr[Udag U] / N
// as well as the width sqrt(<tr^2> - <tr>^2) of its distribution
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Tr[Udag U] / N
void d_link(int bl) {
  register int i, dir;
  register site *s;
  double linktr[NUMLINK], sum = 0.0, linktrSq = 0.0, td;

  for (dir = XUP; dir < NUMLINK; dir++) {
    linktr[dir] = 0.0;
    FORALLSITES(i, s) {
      td = realtrace_su3_f(&(s->linkf[dir]), &(s->linkf[dir]));
      linktr[dir] += td;
      linktrSq += td * td;
    }
    g_doublesum(&(linktr[dir]));
  }
  g_doublesum(&(linktrSq));

  if (bl == 0) {            // Braces suppress compiler complaint
    node0_printf("FLINK");
  }
  else
    node0_printf("BFLINK %d", bl);
  for (dir = XUP; dir < NUMLINK; dir++) {
    linktr[dir] /= ((double)volume * NCOL);
    sum += linktr[dir];
    node0_printf(" %.6g", linktr[dir]);
  }
  sum /= ((double)NUMLINK);
  linktrSq /= ((double)volume * NCOL * NCOL * NUMLINK);
  node0_printf(" %.6g %.6g\n", sum, sqrt(linktrSq - sum * sum));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Tr[Udag U] / N for links in the fermion irrep
void d_link_frep(int bl) {
  register int i, dir;
  register site *s;
  double linktr[NUMLINK], sum = 0.0, linktrSq = 0.0, td;

  for (dir = XUP; dir < NUMLINK; dir++) {
    linktr[dir] = 0.0;
    FORALLSITES(i, s) {
      td = realtrace_su3(&(s->link[dir]), &(s->link[dir]));
      linktr[dir] += td;
      linktrSq += td * td;
    }
    g_doublesum(&(linktr[dir]));
  }

  if (bl == 0) {            // Braces suppress compiler complaint
    node0_printf("ALINK");
  }
  else
    node0_printf("BALINK %d", bl);
  for (dir = XUP; dir < NUMLINK; dir++) {
    linktr[dir] /= ((double)volume * NCOL);
    sum += linktr[dir];
    node0_printf(" %.6g", linktr[dir]);
  }
  sum /= ((double)NUMLINK);
  linktrSq /= ((double)volume * NCOL * NCOL * NUMLINK);
  node0_printf(" %.6g %.6g\n", sum, sqrt(linktrSq - sum * sum));
}
// -----------------------------------------------------------------
