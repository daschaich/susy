// -----------------------------------------------------------------
// Measure average space--space and space--time plaquettes
// Use tempmat for temporary storage
#include "susy_includes.h"

void plaquette(double *ss_plaq, double *st_plaq) {
  register int i, dir, dir2;
  register site *s;
  double ss_sum = 0.0, st_sum = 0.0;
  msg_tag *mtag0, *mtag1;
  matrix tmat;

  // We can exploit a symmetry under dir<-->dir2
  for (dir = TUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      // gen_pt[0] is U_b(x+a), gen_pt[1] is U_a(x+b)
      mtag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                goffset[dir], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                                goffset[dir2], EVENANDODD, gen_pt[1]);

      // tempmat = Udag_b(x) U_a(x)
      FORALLSITES(i, s)
        mult_an(&(s->link[dir2]), &(s->link[dir]), &(tempmat)[i]);
      wait_gather(mtag0);
      wait_gather(mtag1);

      // Compute tr[Udag_a(x+b) Udag_b(x) U_a(x) U_b(x+a)]
      if (dir == TUP || dir2 == TUP) {
        FORALLSITES(i, s) {
          mult_nn(&(tempmat[i]), (matrix *)(gen_pt[0][i]), &tmat);
          st_sum += (double)realtrace((matrix *)(gen_pt[1][i]), &tmat);
        }
      }
      else {
        FORALLSITES(i, s) {
          mult_nn(&(tempmat[i]), (matrix *)(gen_pt[0][i]), &tmat);
          ss_sum += (double)realtrace((matrix *)(gen_pt[1][i]), &tmat);
        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);

  // Average over two plaquettes that involve the temporal link
  // and the one does not
  *ss_plaq = ss_sum / ((double)volume);
  *st_plaq = st_sum / ((double)(2.0 * volume));
}
// -----------------------------------------------------------------
