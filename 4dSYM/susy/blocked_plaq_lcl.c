// -----------------------------------------------------------------
// Measure plaquettes after block RG blocking steps
// Use general_gathers; lattice must be divisible by 2^block in all dirs
// Use tempmat1 and tempmat2 for temporary storage

// #define MIN_PLAQ turns on measurement of minimum plaquette per config
// for tuning smearing as in Hasenfratz & Knechtli, hep-lat/0103029
#define MIN_PLAQ
#include "susy_includes.h"

void blocked_plaq_lcl(int Nsmear, int block) {
  register int i, dir, dir2;
  register site *s;
  register su3_matrix_f *m1, *m4;
  int j, bl = 2, d1[4] = {0, 0, 0, 0}, d2[4] = {0, 0, 0, 0};
  double ss_sum = 0.0, st_sum = 0.0, cur_plaq;
#ifdef MIN_PLAQ
  double min_plaq = 200.0 * NCOL;
#endif
  msg_tag *mtag;
  su3_matrix_f mtmp;

  // Set number of links to stride, bl = 2^block
  // Allow sanity check of reproducing d_plaquette() with this routine
  for (j = 1; j < block; j++)
    bl *= 2;
  if (block <= 0)
    bl = 1;

  // Compute the bl-strided plaquette, exploiting a symmetry under dir<-->dir2
  for (dir = YUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      for (j = 0; j < NDIMS; j++) {
        d1[j] = bl * offset[dir][j];
        d2[j] = bl * offset[dir2][j];
      }
      mtag = start_general_gather_site(F_OFFSET(linkf[dir2]),
                                       sizeof(su3_matrix_f), d1,
                                       EVENANDODD, gen_pt[0]);
      wait_general_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(tempmat2[i]));
      cleanup_general_gather(mtag);

      mtag = start_general_gather_site(F_OFFSET(linkf[dir]),
                                       sizeof(su3_matrix_f), d2,
                                       EVENANDODD, gen_pt[0]);
      FORALLSITES(i, s) {
        m1 = &(s->linkf[dir]);
        m4 = &(s->linkf[dir2]);
        mult_su3_an_f(m4, m1, &(tempmat1[i]));
      }
      wait_general_gather(mtag);

      if (dir == TUP || dir2 == TUP) {
        FORALLSITES(i, s) {
          m1 = &(tempmat2[i]);
          m4 = (su3_matrix_f *)(gen_pt[0][i]);
          mult_su3_nn_f(&(tempmat1[i]), m1, &mtmp);
          cur_plaq = (double)realtrace_su3_f(m4, &mtmp);
#ifdef MIN_PLAQ
          if (cur_plaq < min_plaq)
            min_plaq = cur_plaq;
#endif
          st_sum += cur_plaq;
        }
      }
      else {
        FORALLSITES(i, s) {
          m1 = &(tempmat2[i]);
          m4 = (su3_matrix_f *)(gen_pt[0][i]);
          mult_su3_nn_f(&(tempmat1[i]), m1, &mtmp);
          cur_plaq = (double)realtrace_su3_f(m4, &mtmp);
#ifdef MIN_PLAQ
          if (cur_plaq < min_plaq)
            min_plaq = cur_plaq;
#endif
          ss_sum += cur_plaq;
        }
      }
      cleanup_general_gather(mtag);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);

  // Average over four plaquettes that involve the temporal link
  // and six that do not
#ifdef MIN_PLAQ
  min_plaq = -min_plaq;
  g_doublemax(&min_plaq);
  min_plaq = -min_plaq;
  node0_printf("BMIN_PLAQ %d %d %.8g", Nsmear, block, min_plaq);
#endif
  ss_sum /= ((double)(6.0 * volume));
  st_sum /= ((double)(4.0 * volume));
  node0_printf(" %.8g %.8g\n", ss_sum, st_sum);
}
// -----------------------------------------------------------------
