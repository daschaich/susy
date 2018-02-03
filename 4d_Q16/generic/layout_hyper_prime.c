// -----------------------------------------------------------------
// Routines to determine the distribution of sites on nodes
// Divide the lattice by prime factors in any of the four directions
// Prefers to divide the longest dimensions to minimize surface area
// Prefers to divide dimensions which have already been divided,
// to avoid introducing more off-node directions

// Start trying to divide with the largest prime factor, work down to 2
// The current maximum prime is 53
// The array prime[] may be extended if necessary

// Requires that the lattice volume be divisible by the number of nodes
// Each dimension must be divisible by a suitable factor,
// such that the product of the four factors is the number of nodes
// Also need even number of sites per hypercube

// setup_layout() does any initial setup
//   When it is called the lattice dimensions {nx, ny, nz, nt} have been set
//   This routine sets the global variables "sites_on_node",
//   "even_sites_on_node" and "odd_sites_on_node"
// node_number(x, y, z, t) returns the node number on which a site lives
// node_index(x, y, z, t) returns the index of the site on the node
//   i.e., the site is lattice[node_index(x, y, z, t)]
// num_sites(node) returns the (constant) number of sites on a node
// get_logical_dimensions() returns the machine dimensions
// get_logical_coordinates() returns the mesh coordinates of this node
#include "generic_includes.h"
#ifdef HAVE_QMP
#include <qmp.h>
#endif

static int squaresize[4];           // Dimensions of hypercubes
static int nsquares[4];             // Number of hypercubes in each direction
static int machine_coordinates[4];  // Logical machine coordinates

int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
#define MAXPRIMES (sizeof(prime) / sizeof(int))
// -----------------------------------------------------------------



// -----------------------------------------------------------------
static void setup_hyper_prime() {
  int i, j, k, dir;

  node0_printf("hyper_prime\n");

  // Figure out dimensions of rectangle
  squaresize[XUP] = nx;
  squaresize[YUP] = ny;
  squaresize[ZUP] = nz;
  squaresize[TUP] = nt;
  nsquares[XUP] = 1;
  nsquares[YUP] = 1;
  nsquares[ZUP] = 1;
  nsquares[TUP] = 1;

  i = 1;  // Current number of hypercubes
  while (i < numnodes()) {
    // Figure out which prime to divide by starting with largest
    k = MAXPRIMES - 1;
    while ((numnodes() / i) % prime[k] != 0 && k > 0)
      --k;

    // Figure out which direction to divide
    // Find largest even dimension of h-cubes
    for (j = 1, dir = XUP; dir <= TUP; dir++) {
      if (squaresize[dir] > j && squaresize[dir] % prime[k] == 0)
        j = squaresize[dir];
    }

    // If direction with largest extent has already been divided,
    // divide it again
    // Otherwise divide first direction with largest extent
    FORALLUPDIR(dir) {
      if (squaresize[dir] == j && nsquares[dir] > 1)
        break;
    }
    if (dir > TUP) {
      FORALLUPDIR(dir) {
        if (squaresize[dir] == j)
          break;
      }
    }
    // This can fail if I run out of prime factors in the dimensions
    if (dir > TUP) {
      node0_printf("ERROR: Not enough factors of %d to lay out lattice\n",
                   prime[k]);
      g_sync();
      terminate(1);
    }

    // Do the surgery
    i *= prime[k];
    squaresize[dir] /= prime[k];
    nsquares[dir] *= prime[k];
  }

  // Dividing all-odd lattice among multiple nodes can give squaresize[i]=0
  // Check for that and exit gracefully if encountered
  FORALLUPDIR(dir) {
    if (squaresize[dir] < 1) {
      node0_printf("ERROR: Can't lay out lattice ");
      node0_printf("(bad squaresize in layout_hyper_prime)\n");
      g_sync();
      terminate(1);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_layout() {
  int k = mynode();

  node0_printf("LAYOUT = Hypercubes, options = ");
  setup_hyper_prime();

  // Compute machine coordinates
  machine_coordinates[XUP] = k % squaresize[XUP];
  k /= squaresize[XUP];
  machine_coordinates[YUP] = k % squaresize[YUP];
  k /= squaresize[YUP];
  machine_coordinates[ZUP] = k % squaresize[ZUP];
  k /= squaresize[ZUP];
  machine_coordinates[TUP] = k % squaresize[TUP];

  // Number of sites on node
  sites_on_node = squaresize[XUP] * squaresize[YUP]
                * squaresize[ZUP] * squaresize[TUP];

  // Need even number of sites per hypercube
  if (sites_on_node % 2 != 0) {
    node0_printf("ERROR: Can't lay out lattice ");
    node0_printf("so that all nodes have an even number of sites\n");
    g_sync();
    terminate(1);
  }

  node0_printf("ON EACH NODE %d x %d x %d x %d\n",
               squaresize[XUP], squaresize[YUP],
               squaresize[ZUP], squaresize[TUP]);

  even_sites_on_node = sites_on_node / 2;
  odd_sites_on_node = sites_on_node / 2;
}

int node_number(int x, int y, int z, int t) {
  register int i;
  x /= squaresize[XUP];
  y /= squaresize[YUP];
  z /= squaresize[ZUP];
  t /= squaresize[TUP];
  i = x + nsquares[XUP] * (y + nsquares[YUP] * (z + nsquares[ZUP] * t));
  return i;
}

int node_index(int x, int y, int z, int t) {
  register int i, xr, yr, zr, tr;
  xr = x % squaresize[XUP];
  yr = y % squaresize[YUP];
  zr = z % squaresize[ZUP];
  tr = t % squaresize[TUP];
  i = xr;
  i += squaresize[XUP] * (yr + squaresize[YUP] * (zr + squaresize[ZUP] * tr));
  if ((x + y + z + t) % 2 == 0)   // Even site
    return (i / 2);
  else
    return ((i + sites_on_node) / 2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
size_t num_sites(int node) {
  return sites_on_node;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
const int *get_logical_dimensions() {
  return nsquares;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Coordinates simulate a mesh architecture
// Must correspond to the node_number result
const int *get_logical_coordinate() {
  return machine_coordinates;
}
// -----------------------------------------------------------------
