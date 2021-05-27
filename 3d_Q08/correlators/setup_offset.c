// -----------------------------------------------------------------
// label[k] labels the offset of the kth connection between psibar and psi
// by an integer 1 to 4
// The separation of the paths is in New York metric
// It corresponds to our labeling lambda_1 = lambda(1, 0, 0, 0)
// offset[k] is the actual offset of the kth connection
// This program only generates the positive offsets, 40 of the 80 total
// The others can be found by using adjoints of these paths
// (e.g., (-1, -1, -1, -1) is the adjoint of (1, 1, 1, 1))

#include "corr_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void cubic_neighbor(int x, int y, int z, int t, int *arg, int forw_back,
                    int *xpt, int *ypt, int *zpt, int *tpt) {

  if (forw_back == FORWARDS) {
    *xpt = (x + nx + arg[0]) % nx;
    *ypt = (y + ny + arg[1]) % ny;
    *zpt = (z + nz + arg[2]) % nz;
    *tpt = (t + nt + arg[3]) % nt;
  }
  else {
    *xpt = (x + nx - arg[0]) % nx;
    *ypt = (y + ny - arg[1]) % ny;
    *zpt = (z + nz - arg[2]) % nz;
    *tpt = (t + nt - arg[3]) % nt;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_offset() {
  int i, k;

  // Construct the link paths: one in each direction plus back-diagonal
  for (i = 0; i < NUMLINK - 1; i++) {
    for (k = 0; k < NDIMS; k++)
      offset[i][k] = 0;

    offset[i][i] = 1;
  }
  for (k = 0; k < NDIMS; k++)
    offset[DIR_5][k] = -1;

#ifdef DEBUG_CHECK
  node0_printf("There are %d distinct paths:\n", NUMLINK);
#endif
  // goffset holds indices of gather_array in ../generic/com_mpi.c
  // The first eight elements of gather_array are
  //   XUP, YUP, ZUP, TUP, TDOWN, ZDOWN, YDOWN, XDOWN
  // in that order!
  // In order to use XDOWN = XUP + 1, etc., we make the next ten elements
  //   XUP, XDOWN, YUP, YDOWN, ZUP, ZDOWN, TUP, TDOWN, DIR_5, -DIR_5
  // Then goffset[0]=8, goffset[1]=10, ..., goffset[4]=16
  // But we can't use these in EVEN or ODD gathers!
  for (i = 0; i < NUMLINK; i++) {
    goffset[i] = make_gather(cubic_neighbor, offset[i],
                             WANT_INVERSE, NO_EVEN_ODD, SCRAMBLE_PARITY);

#ifdef DEBUG_CHECK
    int dir;
    node0_printf("  %d ahead:", i);
    for (dir = XUP; dir <= TUP; dir++)
      node0_printf(" %d", offset[i][dir]);

    node0_printf(" (offset %d)\n", goffset[i]);
#endif
  }
}
// -----------------------------------------------------------------
