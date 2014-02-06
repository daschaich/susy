// -----------------------------------------------------------------
// Measure average space--space and space--time plaquettes,
// mallocing the temporary su3_matrix

// #define LOCAL_PLAQ turns on plaquette measurement on hypersurfaces
// plaq_prll[xx] gives the average of the in a fixed MY_X hypersurface,
// as a function of the MY_X coordinate
// Similarly the other plaquettes are given by plaq_perp[xx]

// #define MIN_PLAQ turns on measurement of minimum plaquette per config
// for tuning nHYP as in Hasenfratz & Knechtli, hep-lat/0103029

// #define ALL_PLAQ prints out all irrep plaquettes
// for plotting distribution
// CAUTION: Do not run ALL_PLAQ with MPI!!

#define MIN_PLAQ
#ifdef LOCAL_PLAQ
#define MY_X x
#define MY_DIR XUP          // Printed out first time the routine is called
#define MY_N nx
static int print_dir = 0;   // Controls printing out MY_DIR
#endif

#include "generic_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void d_plaquette_lcl(double *ss_plaq, double *st_plaq) {
  register int i, dir1, dir2;
  register site *s;
  register su3_matrix_f *m1, *m4;
  su3_matrix_f mtmp;
  double ss_sum = 0.0, st_sum = 0.0, cur_plaq;
#ifdef MIN_PLAQ
  double min_plaq = NCOL;
#endif
  msg_tag *mtag0, *mtag1;
  su3_matrix_f *su3mat = malloc(sites_on_node * sizeof(*su3mat));

#ifdef LOCAL_PLAQ
  int xx;
  double *plaq_perp = malloc(MY_N * sizeof(*plaq_perp));
  double *plaq_prll = malloc(MY_N * sizeof(*plaq_prll));

  for (xx = 0; xx < MY_N; xx++) {
    plaq_perp[xx] = 0.0;
    plaq_prll[xx] = 0.0;
  }
#endif

  if (su3mat == NULL) {
    printf("plaquette: can't malloc su3mat\n");
    fflush(stdout);
    terminate(1);
  }

  for (dir1 = YUP; dir1 <= TUP; dir1++) {
    for (dir2 = XUP; dir2 < dir1; dir2++) {
      mtag0 = start_gather_site(F_OFFSET(linkf[dir2]), sizeof(su3_matrix_f),
                                dir1, EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(linkf[dir1]), sizeof(su3_matrix_f),
                                dir2, EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s) {
        m1 = &(s->linkf[dir1]);
        m4 = &(s->linkf[dir2]);
        mult_su3_an_f(m4,m1,&su3mat[i]);
      }
      wait_gather(mtag0);
      wait_gather(mtag1);

      if (dir1 == TUP) {
        FORALLSITES(i, s) {
          mult_su3_nn_f(&(su3mat[i]), (su3_matrix_f *)(gen_pt[0][i]), &mtmp);
          cur_plaq = (double)realtrace_su3_f((su3_matrix_f *)(gen_pt[1][i]),
                                             &mtmp);
#ifdef MIN_PLAQ
          if (cur_plaq < min_plaq)
            min_plaq = cur_plaq;
#endif
          st_sum += cur_plaq;
#ifdef LOCAL_PLAQ
          if (dir1 == MY_DIR || dir2 == MY_DIR)
            plaq_perp[s->MY_X] += cur_plaq;
          else
            plaq_prll[s->MY_X] += cur_plaq;
#endif
        }
      }
      else {
        FORALLSITES(i, s) {
          mult_su3_nn_f(&(su3mat[i]), (su3_matrix_f *)(gen_pt[0][i]), &mtmp);
          cur_plaq = (double)realtrace_su3_f((su3_matrix_f *)(gen_pt[1][i]),
                                             &mtmp);
#ifdef MIN_PLAQ
          if (cur_plaq < min_plaq)
            min_plaq = cur_plaq;
#endif
          ss_sum += cur_plaq;
#ifdef LOCAL_PLAQ
          if (dir1 == MY_DIR || dir2 == MY_DIR)
            plaq_perp[s->MY_X] += cur_plaq;
          else
            plaq_prll[s->MY_X] += cur_plaq;
#endif
        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);

  // Average over three space-space and three time-space plaquettes
  *ss_plaq = ss_sum / ((double)(3.0 * volume));
  *st_plaq = st_sum / ((double)(3.0 * volume));

  free(su3mat);

#ifdef LOCAL_PLAQ
  for (xx = 0; xx < MY_N; xx++) {
    g_doublesum(&(plaq_perp[xx]));
    g_doublesum(&(plaq_prll[xx]));
  }
  // Normalization
  for (xx = 0; xx < MY_N; xx++) {
    plaq_perp[xx] *= ((double)(MY_N)) / (3.0 * volume):
    plaq_prll[xx] *= ((double)(MY_N)) / (3.0 * volume);
  }

  // Print out
  if (this_node == 0) {
    if (print_dir == 0) {
      printf("LOCAL_PLAQ [0=XUP,..., 3=TUP] dir=%d\n", MY_DIR);
      print_dir = 1;
    }
    printf("THIN_PLAQ_PERP");
    for (xx = 0; xx < MY_N; xx++)
      printf(" %e", (double)plaq_perp[xx]);
    printf("\n");
    printf("THIN_PLAQ_PRLL");
    for (xx = 0; xx < MY_N; xx++)
      printf(" %e", (double)plaq_prll[xx]);
    printf("\n");
  }
  free(plaq_perp);
  free(plaq_prll);
#endif

#ifdef MIN_PLAQ
  min_plaq = -min_plaq;
  g_doublemax(&min_plaq);
  min_plaq = -min_plaq;
  node0_printf("MIN_PLAQ %.8g", min_plaq);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifndef PURE_GAUGE
// Irrep plaquette
void d_plaquette_frep_lcl(double *ss_plaq_frep, double *st_plaq_frep) {
  register int i,dir1,dir2;
  register site *s;
  register su3_matrix *m1,*m4;
  su3_matrix mtmp;
  double ss_sum = 0.0, st_sum = 0.0, cur_plaq;
#ifdef MIN_PLAQ
  double min_plaq = DIMF;
#endif
  msg_tag *mtag0,*mtag1;
  su3_matrix *su3mat = malloc(sites_on_node * sizeof(*su3mat));

#ifdef LOCAL_PLAQ
  int xx;
  double *plaq_perp = malloc(MY_N * sizeof(*plaq_perp));
  double *plaq_prll = malloc(MY_N * sizeof(*plaq_prll));

  for (xx = 0; xx < MY_N; xx++) {
    plaq_perp[xx] = 0.0;
    plaq_prll[xx] = 0.0;
  }
#endif

  if (su3mat == NULL) {
    printf("plaquette: can't malloc su3mat\n");
    fflush(stdout);
    terminate(1);
  }

  for (dir1 = YUP; dir1 <= TUP; dir1++) {
    for (dir2 = XUP; dir2 < dir1; dir2++) {
      mtag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(su3_matrix),
                                dir1, EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(link[dir1]), sizeof(su3_matrix),
                                dir2, EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s) {
        m1 = &(s->link[dir1]);
        m4 = &(s->link[dir2]);
        mult_su3_an(m4, m1, &su3mat[i]);
      }

      wait_gather(mtag0);
      wait_gather(mtag1);
      if (dir1 == TUP) {
        FORALLSITES(i, s) {
          mult_su3_nn(&(su3mat[i]), (su3_matrix *)(gen_pt[0][i]), &mtmp);
          cur_plaq = (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                           &mtmp);
#ifdef MIN_PLAQ
          if (cur_plaq < min_plaq)
            min_plaq = cur_plaq;
#endif
#ifdef ALL_PLAQ
          printf("ALL_PLAQ %d %d %d %d %d %d %e\n",
                 s->x, s->y, s->z, s->t, dir1, dir2, cur_plaq);
#endif
          st_sum += cur_plaq;
#ifdef LOCAL_PLAQ
          if (dir1 == MY_DIR || dir2 == MY_DIR)
            plaq_perp[s->MY_X] += cur_plaq;
          else
            plaq_prll[s->MY_X] += cur_plaq;
#endif
        }
      }
      else {
        FORALLSITES(i, s) {
          mult_su3_nn(&(su3mat[i]), (su3_matrix *)(gen_pt[0][i]), &mtmp);
          cur_plaq = (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),
                                           &mtmp);
#ifdef MIN_PLAQ
          if (cur_plaq < min_plaq)
            min_plaq = cur_plaq;
#endif
#ifdef ALL_PLAQ
          printf("ALL_PLAQ %d %d %d %d %d %d %e\n",
                 s->x, s->y, s->z, s->t, dir1, dir2, cur_plaq);
#endif
          ss_sum += cur_plaq;
#ifdef LOCAL_PLAQ
          if (dir1 == MY_DIR || dir2 == MY_DIR)
            plaq_perp[s->MY_X] += cur_plaq;
          else
            plaq_prll[s->MY_X] += cur_plaq;
#endif
        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);

  // Average over three space-space and three time-space plaquettes
  *ss_plaq_frep = ss_sum / ((double)(3.0 * volume));
  *st_plaq_frep = st_sum / ((double)(3.0 * volume));

  free(su3mat);

#ifdef LOCAL_PLAQ
  for (xx = 0; xx < MY_N; xx++) {
    g_doublesum(&(plaq_perp[xx]));
    g_doublesum(&(plaq_prll[xx]));
  }
  // Normalization
  for (xx = 0; xx < MY_N; xx++) {
    plaq_perp[xx] *= ((double)(MY_N)) / (3 * volume);
    plaq_prll[xx] *= ((double)(MY_N)) / (3 * volume);
  }

  // Print out
  if (this_node == 0) {
//    printf("LOCAL_PLAQ dir=%d\n", MY_DIR);
    printf("FAT_PLAQ_PERP");
    for (xx = 0; xx < MY_N; xx++)
      printf(" %.4g", (double)plaq_perp[xx]);
    printf("\n");
    printf("FAT_PLAQ_PRLL");
    for (xx = 0; xx < MY_N; xx++)
      printf(" %.4g", (double)plaq_prll[xx]);
    printf("\n");
  }
  free(plaq_perp);
  free(plaq_prll);
#endif

#ifdef MIN_PLAQ
  min_plaq = -min_plaq;
  g_doublemax(&min_plaq);
  min_plaq = -min_plaq;
  node0_printf("MIN_PLAQ_FERM %e\n",min_plaq);
#endif
}
#endif // Not PURE_GAUGE
// -----------------------------------------------------------------
