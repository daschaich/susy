// -----------------------------------------------------------------
// We are reproducing the C++ constructors:
// Lattice_Vector::Lattice_Vector(void) {
//   for (int i = 0; i < D; i++)
//     coords[i] = 0;
// }
// Lattice_Vector::Lattice_Vector(int mu) {
//   if (mu < D) {
//     for (int i = 0; i < D; i++)
//       coords[i] = 0;
//
//     coords[mu] = 1;
//     return;
//   }
//
//   for (int i = 0; i < D; i++)
//     coords[i] = -1;
//
//   return;
// }

// label[k] labels the offset of the kth connection between psibar and psi
// by an integer 1 to 2
// The separation of the paths is in New York metric
// It corresponds to our labeling of rho and lambda, lambda_1 = lambda(1, 0)
// offset[k] is the actual offset of the kth connection
// This program only generates the positive offsets
// The others can be found by using adjoints of these paths
// (e.g., (-1, -1) is the adjoint of (1, 1))

#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void cubic_neighbor(int x, int t, int *arg, int forw_back,
                    int *xpt, int *tpt) {

  if (forw_back == FORWARDS) {
    *xpt = (x + nx + arg[0]) % nx;
    *tpt = (t + nt + arg[1]) % nt;
  }
  else {
    *xpt = (x + nx - arg[0]) % nx;
    *tpt = (t + nt - arg[1]) % nt;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_bc() {
  register int i, dir;
  register site *s;

  for (dir = 0; dir < NUMLINK; dir++) {
    FORALLSITES(i, s) {
      s->bc[dir] = 1.0;
      s->bc[OPP_LDIR(dir)] = 1.0;

      if (s->t + offset[dir][TUP] < 0)
        s->bc[dir] = PBC;
      else if (s->t + offset[dir][TUP] > nt - 1)
        s->bc[dir] = PBC;
      if (s->t - offset[dir][TUP] < 0)
        s->bc[OPP_LDIR(dir)] = PBC;
      else if (s->t - offset[dir][TUP] > nt - 1)
        s->bc[OPP_LDIR(dir)] = PBC;
    }
  }

  // BC test
//  FORALLSITES(i, s) {
//    if (s->x == 1) {
//      printf("%d %d:", s->x, s->t);
//      for (dir = 0; dir < 2 * NUMLINK; dir++)
//        printf("  %4.2g", s->bc[dir]);
//      printf("\n");
//    }
//  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_offset() {
  int i, k;

  // Construct the link paths: one in each direction plus back-diagonal
  for (i = 0; i < NDIMS; i++) {
    for (k = 0; k < NDIMS; k++)
      offset[i][k] = 0;

    offset[i][i] = 1;
  }

#ifdef DEBUG_CHECK
  node0_printf("There are %d distinct paths:\n", NUMLINK);
#endif
  // Now make the gather tables
  for (i = 0; i < NUMLINK; i++) {
    goffset[i] = make_gather(cubic_neighbor, offset[i],
                             WANT_INVERSE, NO_EVEN_ODD, SCRAMBLE_PARITY);

#ifdef DEBUG_CHECK
    int mu;
    node0_printf("  %d ahead:", i);
    for (mu = 0; mu < NDIMS; mu++)
      node0_printf(" %d", offset[i][mu]);

    node0_printf(" (offset %d)\n", goffset[i]);
#endif
  }

  // Make boundary condition tables
  if (PBC > 0.0)
    node0_printf("Periodic temporal boundary conditions\n");
  if (PBC < 0.0)
    node0_printf("Antiperiodic temporal boundary conditions\n");
  setup_bc();
}
// -----------------------------------------------------------------
