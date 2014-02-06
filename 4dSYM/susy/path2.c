// -----------------------------------------------------------------
// Walk around specified path of fundamental links
// This version includes a "kind" array for modified Wilson loops:
// kind[i] = 0 --> linkf, kind[i] = 1 --> mom = (linkf^{-1})^dag
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dir is a list of the directions in the path, with the given length
// sign is the corresponding list of which way to go in the given dir
// That is, the negative sign means take the adjoint
// kind tells us whether to use linkf or mom = (linkf^{-1})^dag
// Uses tempmat1 to accumulate linkf product along path
void path2(int *dir, int *sign, int *kind, int length) {
  register int i;
  register site *s;
  int j;
  msg_tag *mtag0 = NULL;

  // Initialize tempmat1 with first link in path
  if (sign[0] > 0) {    // Gather forwards, no adjoint
    if (kind[0] == 0) {
      mtag0 = start_gather_site(F_OFFSET(linkf[dir[0]]), sizeof(su3_matrix_f),
                                goffset[dir[0]]+1, EVENANDODD, gen_pt[0]);
    }
    else if (kind[0] == 1) {
      mtag0 = start_gather_site(F_OFFSET(mom[dir[0]]), sizeof(su3_matrix_f),
                                goffset[dir[0]]+1, EVENANDODD, gen_pt[0]);
    }
    else {
      node0_printf("path2: unrecognized kind[0] = %d\n", kind[0]);
      terminate(1);
    }

    wait_gather(mtag0);
    FORALLSITES(i, s)
      su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(s->tempmat1));

    cleanup_gather(mtag0);
  }

  if (sign[0] < 0) {    // Take adjoint, no gather
    FORALLSITES(i, s) {
      if (kind[0] == 0)
        su3_adjoint_f(&(s->linkf[dir[0]]), &(s->tempmat1));
      else if (kind[0] == 1)
        su3_adjoint_f(&(s->mom[dir[0]]), &(s->tempmat1));
      else {
        node0_printf("path2: unrecognized kind[0] = %d\n", kind[0]);
        terminate(1);
      }
    }
  }

  // Accumulate subsequent links in product in tempmat1
  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {    // mult_su3_nn_f then gather backwards
      FORALLSITES(i, s) {
        if (kind[j] == 0)
          mult_su3_nn_f(&(s->tempmat1), &(s->linkf[dir[j]]), &(s->tempmat2));
        else if (kind[j] == 1)
          mult_su3_nn_f(&(s->tempmat1), &(s->mom[dir[j]]), &(s->tempmat2));
        else {
          node0_printf("path2: unrecognized kind[%d] = %d\n", j, kind[j]);
          terminate(1);
        }
      }
      mtag0 = start_gather_site(F_OFFSET(tempmat2), sizeof(su3_matrix_f),
                                goffset[dir[j]] + 1, EVENANDODD, gen_pt[0]);

      wait_gather(mtag0);
      FORALLSITES(i, s)
        su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(s->tempmat1));

      cleanup_gather(mtag0);
    }

    if (sign[j] < 0) {    // Gather forwards then mult_su3_na_f
      mtag0 = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix_f),
                                goffset[dir[j]], EVENANDODD, gen_pt[1]);

      wait_gather(mtag0);
      FORALLSITES(i, s) {
        if (kind[j] == 0) {
          mult_su3_na_f((su3_matrix_f *)(gen_pt[1][i]), &(s->linkf[dir[j]]),
                        &(s->tempmat2));
        }
        else if (kind[j] == 1) {
          mult_su3_na_f((su3_matrix_f *)(gen_pt[1][i]), &(s->mom[dir[j]]),
                        &(s->tempmat2));
        }
        else {
          node0_printf("path2: unrecognized kind[%d] = %d\n", j, kind[j]);
          terminate(1);
        }
      }
      FORALLSITES(i, s)   // Don't want to overwrite tempmat1 too soon
        su3mat_copy_f(&(s->tempmat2), &(s->tempmat1));

      cleanup_gather(mtag0);
    }
  }
}
// -----------------------------------------------------------------
