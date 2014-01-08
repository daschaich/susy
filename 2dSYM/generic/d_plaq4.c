// -----------------------------------------------------------------
// Measure average plaquette, mallocing the temporary su3_matrix
#include "generic_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void d_plaquette(double *plaq) {
  register int i;
  register site *s;
  register su3_matrix_f *m1, *m4;
  su3_matrix_f mtmp;
  double sum = 0.0;
  msg_tag *mtag0, *mtag1;
  su3_matrix_f *su3mat = malloc(sites_on_node * sizeof(*su3mat));

  if (su3mat == NULL) {
    printf("plaquette: can't malloc su3mat\n");
    fflush(stdout);
    terminate(1);
  }

  mtag0 = start_gather_site(F_OFFSET(linkf[XUP]), sizeof(su3_matrix_f),
                            TUP, EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_site(F_OFFSET(linkf[TUP]), sizeof(su3_matrix_f),
                            XUP, EVENANDODD, gen_pt[1]);

  FORALLSITES(i, s) {
    m1 = &(s->linkf[TUP]);
    m4 = &(s->linkf[XUP]);
    mult_su3_an_f(m4, m1, &su3mat[i]);
  }
  wait_gather(mtag0);
  wait_gather(mtag1);

  FORALLSITES(i, s) {
    mult_su3_nn_f(&(su3mat[i]), (su3_matrix_f *)(gen_pt[0][i]), &mtmp);
    sum += (double)realtrace_su3_f((su3_matrix_f *)(gen_pt[1][i]), &mtmp);
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);
  g_doublesum(&sum);

  // Average over volume
  *plaq = sum / ((double)volume);

  free(su3mat);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifndef PURE_GAUGE
// Irrep plaquette
void d_plaquette_frep(double *plaq_frep) {
  register int i;
  register site *s;
  register su3_matrix *m1, *m4;
  su3_matrix mtmp;
  double sum = 0.0;
  msg_tag *mtag0, *mtag1;
  su3_matrix *su3mat = malloc(sites_on_node * sizeof(*su3mat));

  if (su3mat == NULL) {
    printf("plaquette: can't malloc su3mat\n");
    fflush(stdout);
    terminate(1);
  }

  mtag0 = start_gather_site(F_OFFSET(link[XUP]), sizeof(su3_matrix),
                            TUP, EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_site(F_OFFSET(link[TUP]), sizeof(su3_matrix),
                            XUP, EVENANDODD, gen_pt[1]);

  FORALLSITES(i, s) {
    m1 = &(s->link[TUP]);
    m4 = &(s->link[XUP]);
    mult_su3_an(m4, m1, &su3mat[i]);
  }

  wait_gather(mtag0);
  wait_gather(mtag1);
  FORALLSITES(i, s) {
    mult_su3_nn(&(su3mat[i]), (su3_matrix *)(gen_pt[0][i]), &mtmp);
    sum += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]), &mtmp);
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);
  g_doublesum(&sum);

  // Average over volmue
  *plaq_frep = sum / ((double)volume);

  free(su3mat);
}
#endif // Not PURE_GAUGE
// -----------------------------------------------------------------
