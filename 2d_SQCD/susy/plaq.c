// -----------------------------------------------------------------
// Measure average space--space and space--time plaquettes
// Use tempmat for temporary storage
#include "susy_includes.h"

void plaquette(double *plaq) {
  register int i;
  register site *s;
  double sum = 0.0;
  msg_tag *mtag0, *mtag1;
  matrix tmat;

  // We can exploit a symmetry under TUP<-->XUP
  // gen_pt[0] is U_b(x+a), gen_pt[1] is U_a(x+b)
  mtag0 = start_gather_site(F_OFFSET(link[XUP]), sizeof(matrix),
                            goffset[TUP], EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_site(F_OFFSET(link[TUP]), sizeof(matrix),
                            goffset[XUP], EVENANDODD, gen_pt[1]);

  // tempmat = Udag_b(x) U_a(x)
  FORALLSITES(i, s)
    mult_an(&(s->link[XUP]), &(s->link[TUP]), &(tempmat)[i]);

  // Compute tr[Udag_a(x+b) Udag_b(x) U_a(x) U_b(x+a)]
  wait_gather(mtag0);
  wait_gather(mtag1);
  FORALLSITES(i, s) {
    mult_nn(&(tempmat[i]), (matrix *)(gen_pt[0][i]), &tmat);
    sum += (double)realtrace((matrix *)(gen_pt[1][i]), &tmat);
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);
  g_doublesum(&sum);

  // Average over volume
  *plaq = sum / ((double)volume);
}
// -----------------------------------------------------------------
