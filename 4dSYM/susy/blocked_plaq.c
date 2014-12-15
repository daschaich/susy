// -----------------------------------------------------------------
// Evaluate the plaquette after block RG blocking steps
// Use general_gathers; lattice must be divisible by 2^block in all dirs
#include "susy_includes.h"

void blocked_plaq(int Nstout, int block) {
  register int i, dir, dir2;
  register site *s;
  register su3_matrix_f *m1, *m4;
  int j, bl = 2, d1[4] = {0, 0, 0, 0}, d2[4] = {0, 0, 0, 0};
  double ss_sum = 0.0, st_sum = 0.0, tr;
  msg_tag *mtag;
  su3_matrix_f tmat;
  su3_matrix_f *su3mat = malloc(sites_on_node * sizeof(*su3mat));

  if (su3mat == NULL) {
    printf("blocked_plaq: can't malloc su3mat\n");
    fflush(stdout);
    terminate(1);
  }

  // Set number of links to stride, bl = 2^block
  // Allow sanity check of reproducing ploop() with this routine
  for (j = 1; j < block; j++)
    bl *= 2;
  if (block <= 0)
    bl = 1;

  // Copy temporal links to tempmat1
  FORALLSITES(i, s)
    su3mat_copy_f(&(s->linkf[TUP]), &(s->tempmat1));

  // Compute the bl-strided plaquette, exploiting a symmetry under dir<-->dir2
  for (dir = YUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      for (j = 0; j < NDIMS; j++) {
        d1[j] = bl * offset[dir][j];
        d2[j] = bl * offset[dir2][j];
      }
      // Can only have one general gather at once...
      mtag = start_general_gather_site(F_OFFSET(linkf[dir2]),
                                        sizeof(su3_matrix_f), d1,
                                        EVENANDODD, gen_pt[0]);

      // su3mat = Udag_b(x) U_a(x)
      FORALLSITES(i, s) {
        m1 = &(s->linkf[dir]);
        m4 = &(s->linkf[dir2]);
        mult_su3_an_f(m4, m1, &su3mat[i]);
      }

      // Copy first gather to tempmat1
      wait_general_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(s->tempmat1));
      cleanup_general_gather(mtag);

      mtag = start_general_gather_site(F_OFFSET(linkf[dir]),
                                        sizeof(su3_matrix_f), d2,
                                        EVENANDODD, gen_pt[0]);
      wait_general_gather(mtag);

      // Compute tr[Udag_a(x+d2) Udag_b(x) U_a(x) U_b(x+d1)]
      if (dir == TUP || dir2 == TUP) {
        FORALLSITES(i, s) {
          mult_su3_nn_f(&(su3mat[i]), &(s->tempmat1), &tmat);
          st_sum += (double)realtrace_su3_f((su3_matrix_f *)(gen_pt[0][i]),
                                            &tmat);
        }
      }
      else {
        FORALLSITES(i, s) {
          mult_su3_nn_f(&(su3mat[i]), &(s->tempmat1), &tmat);
          ss_sum += (double)realtrace_su3_f((su3_matrix_f *)(gen_pt[0][i]),
                                            &tmat);
        }
      }
      cleanup_general_gather(mtag);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);

  // Average over four plaquettes that involve the temporal link
  // and six that do not
  ss_sum /= ((double)(6.0 * volume));
  st_sum /= ((double)(4.0 * volume));
  tr = (ss_sum + st_sum) / 2.0;
  node0_printf("BPLAQ %d %d %.8g %.8g %.8g\n",
               Nstout, block, ss_sum, st_sum, tr);

  free(su3mat);
}
// -----------------------------------------------------------------
