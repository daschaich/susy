// -----------------------------------------------------------------
// Measure the average value of the link Tr[Udag U] / N
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Tr[Udag U] / N
void d_link() {
  register int i, dir;
  register site *s;
  double g_action[NUMLINK], g_sum = 0;

  for (dir = XUP; dir < NUMLINK; dir++) {
    g_action[dir] = 0;
    FORALLSITES(i, s)
      g_action[dir] += realtrace_su3_f(&(s->linkf[dir]), &(s->linkf[dir]))
                       / ((double)(NCOL));
    g_doublesum(&(g_action[dir]));
  }

  node0_printf("FLINK");
  for (dir = XUP; dir < NUMLINK; dir++) {
    g_action[dir] /= ((double)volume);
    g_sum += g_action[dir];
    node0_printf(" %.6g", g_action[dir]);
  }
  node0_printf(" %.6g\n", g_sum / ((double)(NUMLINK)));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Tr[Udag U] / N for links in the fermion irrep
void d_link_frep() {
  register int i, dir;
  register site *s;
  double g_action[NUMLINK], g_sum = 0.0;

  for (dir = XUP; dir < NUMLINK; dir++) {
    g_action[dir] = 0.0;
    FORALLSITES(i, s)
      g_action[dir] += realtrace_su3(&(s->link[dir]), &(s->link[dir]))
                       / ((double)(DIMF));
    g_doublesum(&(g_action[dir]));
  }

  node0_printf("ALINK");
  for (dir = XUP; dir < NUMLINK; dir++) {
    g_action[dir] /= ((double)volume);
    g_sum += g_action[dir];
    node0_printf(" %.6g", g_action[dir]);
  }
  node0_printf(" %.6g\n", g_sum / ((double)NUMLINK));
}
// -----------------------------------------------------------------
