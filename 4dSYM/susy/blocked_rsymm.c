// -----------------------------------------------------------------
// Modified rectangular Wilson loops after RG blocking
// Use general_gathers; lattice must be divisible by 2^block in all dirs
// Evaluate in different spatial dirs to check rotational invariance
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Walk around strided path of links specified by dir, sign and kind
// dir is a list of the directions in the path, with the given length
// sign is the corresponding list of which way to go in the given dir
// That is, the negative sign means take the adjoint
// bl = 2^block is number of links to stride
// Use tempmat to accumulate link product along path
void blocked_path(int *dir, int *sign, int length, int bl) {
  register int i;
  register site *s;
  int j, k, d[4];
  msg_tag *tag;
  matrix *mat;

  // Initialize tempmat with first link in path
  if (sign[0] > 0) {    // Gather from site - bl * dir[0], no adjoint
    for (k = 0; k < NDIMS; k++)
      d[k] = -bl * offset[dir[0]][k];
    tag = start_general_gather_site(F_OFFSET(link[dir[0]]), sizeof(matrix),
                                    d, EVENANDODD, gen_pt[0]);

    wait_general_gather(tag);
    FORALLSITES(i, s)
      mat_copy_f((matrix *)(gen_pt[0][i]), &(tempmat[i]));

    cleanup_general_gather(tag);
  }

  if (sign[0] < 0) {    // Take adjoint, no gather
    FORALLSITES(i, s)
      adjoint_f(&(s->link[dir[0]]), &(tempmat[i]));
  }

  // Accumulate subsequent links in product in tempmat
  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {    // mult_nn then gather from site - bl * dir[j]
      FORALLSITES(i, s)
        mult_nn(&(tempmat[i]), &(s->link[dir[j]]), &(tempmat2[i]));

      for (k = 0; k < NDIMS; k++)
        d[k] = -bl * offset[dir[j]][k];
      tag = start_general_gather_field(tempmat2, sizeof(matrix),
                                       d, EVENANDODD, gen_pt[0]);

      wait_general_gather(tag);
      FORALLSITES(i, s)
        mat_copy_f((matrix *)(gen_pt[0][i]), &(tempmat[i]));

      cleanup_general_gather(tag);
    }

    if (sign[j] < 0) {    // Gather from site + bl * dir[j] then mult_na
      for (k = 0; k < NDIMS; k++)
        d[k] = bl * offset[dir[j]][k];
      tag = start_general_gather_field(tempmat, sizeof(matrix),
                                       d, EVENANDODD, gen_pt[1]);

      wait_general_gather(tag);
      FORALLSITES(i, s) {
        mat = (matrix *)(gen_pt[1][i]);
        mult_na(mat, &(s->link[dir[j]]), &(tempmat2[i]));
      }
      FORALLSITES(i, s)   // Don't want to overwrite tempmat too soon
        mat_copy_f(&(tempmat2[i]), &(tempmat[i]));

      cleanup_general_gather(tag);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Walk around strided path of links specified by dir, sign and kind
// dir lists the directions in the path
// sign lists whether to go forward (1) or backwards (-1)
// kind tells us whether to use link (1) or mom = (link^{-1})^dag (-1)
// length is the length of the path, and of each array
// bl = 2^block is number of links to stride
// Use tempmat to accumulate link product along path
void blocked_rsymm_path(int *dir, int *sign, int *kind, int length, int bl) {
  register int i;
  register site *s;
  int j, k, d[4] = {0, 0, 0, 0};
  msg_tag *tag = NULL;
  matrix *mat;

  // Initialize tempmat with first link in path
  if (sign[0] > 0) {    // Gather from site - bl * dir[0], no adjoint
    for (k = 0; k < NDIMS; k++)
      d[k] = -bl * offset[dir[0]][k];
    if (kind[0] > 0) {
      tag = start_general_gather_site(F_OFFSET(link[dir[0]]),
                                      sizeof(matrix), d,
                                      EVENANDODD, gen_pt[0]);
    }
    else if (kind[0] < 0) {
      tag = start_general_gather_site(F_OFFSET(mom[dir[0]]),
                                      sizeof(matrix), d,
                                      EVENANDODD, gen_pt[0]);
    }
    else {
      node0_printf("blocked_rsymm_path: unrecognized kind[0] = %d\n",
                   kind[0]);
      terminate(1);
    }

    wait_general_gather(tag);
    FORALLSITES(i, s)
      mat_copy_f((matrix *)(gen_pt[0][i]), &(tempmat[i]));
    cleanup_general_gather(tag);
  }

  else if (sign[0] < 0) {    // Take adjoint, no gather
    FORALLSITES(i, s) {
      if (kind[0] > 0)
        adjoint_f(&(s->link[dir[0]]), &(tempmat[i]));
      else if (kind[0] < 0)
        adjoint_f(&(s->mom[dir[0]]), &(tempmat[i]));
      else {
        node0_printf("blocked_rsymm_path: unrecognized kind[0] = %d\n",
                     kind[0]);
        terminate(1);
      }
    }
  }
  else {
    node0_printf("blocked_rsymm_path: unrecognized sign[0] = %d\n", sign[0]);
    terminate(1);
  }

  // Accumulate subsequent links in product in tempmat
  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {    // mult_nn then gather from site - bl * dir[j]
      FORALLSITES(i, s) {
        if (kind[j] > 0)
          mult_nn(&(tempmat[i]), &(s->link[dir[j]]), &(tempmat2[i]));
        else if (kind[j] < 0)
          mult_nn(&(tempmat[i]), &(s->mom[dir[j]]), &(tempmat2[i]));
        else {
          node0_printf("blocked_rsymm_path: unrecognized kind[%d] = %d\n",
                       j, kind[j]);
          terminate(1);
        }
      }
      for (k = 0; k < NDIMS; k++)
        d[k] = -bl * offset[dir[j]][k];
      tag = start_general_gather_field(tempmat2, sizeof(matrix),
                                       d, EVENANDODD, gen_pt[0]);

      wait_general_gather(tag);
      FORALLSITES(i, s)
        mat_copy_f((matrix *)(gen_pt[0][i]), &(tempmat[i]));
      cleanup_general_gather(tag);
    }

    // Gather from site + bl * dir[j] then mult_na
    else if (sign[j] < 0) {
      for (k = 0; k < NDIMS; k++)
        d[k] = bl * offset[dir[j]][k];
      tag = start_general_gather_field(tempmat, sizeof(matrix),
                                       d, EVENANDODD, gen_pt[1]);

      wait_general_gather(tag);
      FORALLSITES(i, s) {
        if (kind[j] > 0) {
          mat = (matrix *)(gen_pt[1][i]);
          mult_na(mat, &(s->link[dir[j]]), &(tempmat2[i]));
        }
        else if (kind[j] < 0) {
          mat = (matrix *)(gen_pt[1][i]);
          mult_na(mat, &(s->mom[dir[j]]), &(tempmat2[i]));
        }
        else {
          node0_printf("blocked_rsymm_path: unrecognized kind[%d] = %d\n",
                       j, kind[j]);
          terminate(1);
        }
      }
      FORALLSITES(i, s)   // Don't want to overwrite tempmat too soon
        mat_copy_f(&(tempmat2[i]), &(tempmat[i]));
      cleanup_general_gather(tag);
    }
    else {
      node0_printf("blocked_rsymm_path: unrecognized sign[%d] = %d\n",
                   j, sign[j]);
      terminate(1);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print both usual and transformed Wilson loops
void blocked_rsymm(int Nsmear, int block) {
  register int i;
  register site *s;
  int j, bl = 2, max, dir_normal, dir_inv, dist, dist_inv, mu, length;
  double rsymm_loop, wloop, invlink[NUMLINK], invlink_sum = 0.0, td;
  double invlinkSq = 0.0;
  complex tc;
  matrix tmat;

  // Set number of links to stride, bl = 2^block
  // Allow sanity check of reproducing rsymm() with this routine
  for (j = 1; j < block; j++)
    bl *= 2;
  if (block <= 0)
    bl = 1;

  max = (MAX_X + 1) / bl;
  node0_printf("blocked_rsymm: MAX = %d\n", max);
  int dir[4 * max], sign[4 * max], kind[4 * max];

  // Compute and optionally check inverse matrices
  // Temporarily store the adjoint of the inverse in momentum matrices,
  // since it transforms like the original link
  FORALLDIR(mu) {
    FORALLSITES(i, s) {
      invert(&(s->link[mu]), &tmat);
      adjoint_f(&tmat, &(s->mom[mu]));

#ifdef DEBUG_CHECK
#define INV_TOL 1e-12
#define INV_TOL_SQ 1e-24
      // Check inversion -- tmat should be unit matrix
      int j, k;
      mult_nn(&(s->mom[dir]), &(s->link[dir]), &tmat);
      for (j = 0; j < NCOL; j++) {
        if (fabs(1 - tmat.e[j][j].real) > INV_TOL
            || fabs(tmat.e[j][j].imag) > INV_TOL) {
          printf("Link inversion fails on node%d:\n", this_node);
          dumpmat(&tmat);
        }
        for (k = j + 1; k < NCOL; k++) {
          if (cabs_sq(&(tmat.e[j][k])) > INV_TOL_SQ
              || cabs_sq(&(tmat.e[k][j])) > INV_TOL_SQ) {
            printf("Link inversion fails on node%d:\n", this_node);
            dumpmat(&tmat);
          }
        }
      }
      // Check left multiplication in addition to right multiplication
      mult_nn(&(s->link[dir]), &(s->mom[dir]), &tmat);
      for (j = 0; j < NCOL; j++) {
        if (fabs(1 - tmat.e[j][j].real) > INV_TOL
            || fabs(tmat.e[j][j].imag) > INV_TOL) {
          printf("Link inversion fails on node%d:\n", this_node);
          dumpmat(&tmat);
        }
        for (k = j + 1; k < NCOL; k++) {
          if (cabs_sq(&(tmat.e[j][k])) > INV_TOL_SQ
              || cabs_sq(&(tmat.e[k][j])) > INV_TOL_SQ) {
            printf("Link inversion fails on node%d:\n", this_node);
            dumpmat(&tmat);
          }
        }
      }
#endif
    }
  }

  // First check average value of the inverted link
  // Tr[U^{-1} (U^{-1})^dag] / N and corresponding width
  // Just like link_trace() but use s->mom instead of s->link
  for (dir_inv = XUP; dir_inv < NUMLINK; dir_inv++) {
    invlink[dir_inv] = 0.0;
    FORALLSITES(i, s) {
      td = realtrace(&(s->mom[dir_inv]), &(s->mom[dir_inv]));
      invlink[dir_inv] += td;
      invlinkSq += td * td;
    }
    invlink[dir_inv] *= one_ov_N / ((double)volume);
    g_doublesum(&(invlink[dir_inv]));
    invlink_sum += invlink[dir_inv];
  }
  invlink_sum /= ((double)NUMLINK);
  invlinkSq *= one_ov_N * one_ov_N / ((double)volume * NUMLINK);
  g_doublesum(&(invlinkSq));

  node0_printf("BINVLINK %d %d", Nsmear, block);
  for (dir_inv = XUP; dir_inv < NUMLINK; dir_inv++)
    node0_printf(" %.6g", invlink[dir_inv]);
  td = sqrt(invlinkSq - invlink_sum * invlink_sum);
  node0_printf(" %.6g %.6g\n", invlink_sum, td);

  // Construct and print all loops up to max x max
  // in all NUMLINK * (NUMLINK - 1) directions
  // Invert all links in the second direction in each loop
  for (dir_normal = XUP; dir_normal < NUMLINK; dir_normal++) {
    for (dir_inv = XUP; dir_inv < NUMLINK; dir_inv++) {
      if (dir_inv == dir_normal)
        continue;

      for (dist = 1; dist <= max; dist++) {
        for (dist_inv = 1; dist_inv <= max; dist_inv++) {
          // Set up rectangular Wilson loop path as list of dir, sign * kind
          length = 2 * (dist + dist_inv);
          for (i = 0; i < dist; i++) {
            dir[i] = dir_normal;
            sign[i] = 1;
            kind[i] = 1;
          }
          for (i = dist; i < dist + dist_inv; i++) {
            dir[i] = dir_inv;
            sign[i] = 1;
            kind[i] = -1;
          }
          for (i = dist + dist_inv; i < 2 * dist + dist_inv; i++) {
            dir[i] = dir_normal;
            sign[i] = -1;
            kind[i] = 1;
          }
          for (i = 2 * dist + dist_inv; i < length; i++) {
            dir[i] = dir_inv;
            sign[i] = -1;
            kind[i] = -1;
          }
#ifdef DEBUG_CHECK
          node0_printf("path %d [%d] %d [%d] length %d: ",
                       dir_normal, dist, dir_inv, dist_inv, length);
          for (i = 0; i < length; i++)
            node0_printf(" (%d)*%d(%d) ", sign[i], dir[i], kind[i]);
          node0_printf("\n");
#endif

          // blocked_path and blocked_rsymm_path
          // both accumulate the product in tempmat
          blocked_path(dir, sign, length, bl);
          wloop = 0.0;
          FORALLSITES(i, s) {
            tc = trace_f(&(tempmat[i]));
            wloop += tc.real;
          }
          g_doublesum(&wloop);

          blocked_rsymm_path(dir, sign, kind, length, bl);
          rsymm_loop = 0.0;
          FORALLSITES(i, s) {
            tc = trace_f(&(tempmat[i]));
            rsymm_loop += tc.real;
          }
          g_doublesum(&rsymm_loop);
          // Format: normal [dir] inverted [dir] usual transformed
          node0_printf("BRSYMM %d %d %d [%d] %d [%d] %.8g %.8g\n",
                       Nsmear, block, dist, dir_normal, dist_inv, dir_inv,
                       wloop / volume, rsymm_loop / volume);
        } // dist_inv
      } // dist
    } // dir_inv
  } // dir_normal
}
// -----------------------------------------------------------------
