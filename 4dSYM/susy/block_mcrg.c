// -----------------------------------------------------------------
// Blocking function for MCRG
// Just multiplies together "neighboring" links separated by 2^block sites
// Saves all 2^d products in sub-lattices of original gauge configuration
// Overwrite links, use site->mom for temporary storage
#include "susy_includes.h"

void block_mcrg(int block) {
  register int dir, i;
  register site *s;
  int disp[4], j, bl = 1;
  msg_tag *tag;

  // Set number of links to stride, bl = 2^(block - 1)
  for (j = 1; j < block; j++)
    bl *= 2;
  for (dir = XUP; dir < NUMLINK; dir++) {
    for (j = 0; j < NDIMS; j++)
      disp[j] = bl * offset[dir][j];

    tag = start_general_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[0]);
    wait_general_gather(tag);

    FORALLSITES(i, s)
      mult_nn(&(s->link[dir]), (matrix *)(gen_pt[0][i]), &(s->mom[dir]));
    cleanup_general_gather(tag);
  }

  // Overwrite original links
  FORALLSITES(i, s) {
    FORALLDIR(dir)
      mat_copy_f(&(s->mom[dir]), &(s->link[dir]));
  }
}
// -----------------------------------------------------------------
