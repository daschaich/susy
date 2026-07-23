// -----------------------------------------------------------------
// Measure average space--space and space--time plaquettes
// Use tempmat for temporary storage

// #define MIN_PLAQ turns on measurement of minimum plaquette per config
// #define PLAQ_DIST prints out all plaquettes for plotting distribution
// CAUTION: Do not run PLAQ_DIST with MPI!

#define MIN_PLAQ
//#define PLAQ_DIST
#include "susy_includes.h"

double local_plaquette(double *plaq) {
  register int i;
  register site *s;
  double sum = 0.0, cur_plaq, max_plaq = 0.0;
#ifdef MIN_PLAQ
  double min_plaq = 200.0 * NCOL;
#endif
  msg_tag *mtag0, *mtag1;
  matrix tmat;

#ifdef PLAQ_DIST
  if (this_node != 0) {
    printf("plaquette: don't run PLAQ_DIST in parallel\n");
    fflush(stdout);
    terminate(1);
  }
#endif

  mtag0 = start_gather_site(F_OFFSET(link[XUP]), sizeof(matrix),
                            goffset[TUP], EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_site(F_OFFSET(link[TUP]), sizeof(matrix),
                            goffset[XUP], EVENANDODD, gen_pt[1]);

  FORALLSITES(i, s)
    mult_an(&(s->link[XUP]), &(s->link[TUP]), &(tempmat[i]));
  wait_gather(mtag0);
  wait_gather(mtag1);

  FORALLSITES(i, s) {
    mult_nn(&(tempmat[i]), (matrix *)(gen_pt[0][i]), &tmat);
    cur_plaq = (double)realtrace((matrix *)(gen_pt[1][i]), &tmat);
    if (cur_plaq > max_plaq)
      max_plaq = cur_plaq;
#ifdef MIN_PLAQ
    if (cur_plaq < min_plaq)
      min_plaq = cur_plaq;
#endif
#ifdef PLAQ_DIST
    printf("PLAQ_DIST %d %d %.4g\n", s->x, s->t, cur_plaq);
#endif
    sum += cur_plaq;
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);
  g_doublesum(&sum);

  // Average over volume
  *plaq = sum / ((double)(volume));

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
