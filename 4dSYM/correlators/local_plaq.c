// -----------------------------------------------------------------
// Measure average space--space and space--time plaquettes
// Use tempmat for temporary storage
// #define MIN_PLAQ turns on measurement of minimum plaquette per config
#define MIN_PLAQ
#include "corr_includes.h"

double local_plaquette(double *ss_plaq, double *st_plaq) {
  register int i, dir, dir2;
  register site *s;
  double ss_sum = 0.0, st_sum = 0.0, cur_plaq, max_plaq = 0.0;
#ifdef MIN_PLAQ
  double min_plaq = 200.0 * NCOL;
#endif
  msg_tag *mtag0, *mtag1;
  matrix tmat;

  // We can exploit a symmetry under dir<-->dir2
  for (dir = YUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      mtag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                goffset[dir], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                                goffset[dir2], EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s)
        mult_an(&(s->link[dir2]), &(s->link[dir]), &(tempmat[i]));
      wait_gather(mtag0);
      wait_gather(mtag1);

      if (dir == TUP || dir2 == TUP) {
        FORALLSITES(i, s) {
          mult_nn(&(tempmat[i]), (matrix *)(gen_pt[0][i]), &tmat);
          cur_plaq = (double)realtrace((matrix *)(gen_pt[1][i]), &tmat);
          if (cur_plaq > max_plaq)
            max_plaq = cur_plaq;
#ifdef MIN_PLAQ
          if (cur_plaq < min_plaq)
            min_plaq = cur_plaq;
#endif
          st_sum += cur_plaq;
        }
      }
      else {
        FORALLSITES(i, s) {
          mult_nn(&(tempmat[i]), (matrix *)(gen_pt[0][i]), &tmat);
          cur_plaq = (double)realtrace((matrix *)(gen_pt[1][i]), &tmat);
          if (cur_plaq > max_plaq)
            max_plaq = cur_plaq;
#ifdef MIN_PLAQ
          if (cur_plaq < min_plaq)
            min_plaq = cur_plaq;
#endif
          ss_sum += cur_plaq;
        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);

  // Average over four plaquettes that involve the temporal link
  // and six that do not
  *ss_plaq = ss_sum / ((double)(6.0 * volume));
  *st_plaq = st_sum / ((double)(4.0 * volume));

  g_doublemax(&max_plaq);
#ifdef MIN_PLAQ
  // Somewhat hacky since we don't have g_doublemin...
  min_plaq = -min_plaq;
  g_doublemax(&min_plaq);
  min_plaq = -min_plaq;
  node0_printf("MIN_PLAQ %.8g", min_plaq);
#endif
  return max_plaq;
}
// -----------------------------------------------------------------
