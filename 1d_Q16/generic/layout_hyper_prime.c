// -----------------------------------------------------------------
// Routines to determine the distribution of sites on nodes
// Divide the lattice by prime factors

// Start trying to divide with the largest prime factor, work down to 2
// The current maximum prime is 53
// The array prime[] may be extended if necessary

// Requires that nt be divisible by the number of cores
// Also need even number of sites per core

// setup_layout() does any initial setup
//   When it is called nt has been set
//   This routine sets the global variables "sites_on_node",
//   "even_sites_on_node" and "odd_sites_on_node"
// node_number(t) returns the node number on which a site lives
// node_index(t) returns the index of the site on the node
//   i.e., the site is lattice[node_index(t)]
// num_sites(node) returns the (constant) number of sites on a node
// get_logical_dimensions() returns the machine dimensions
// get_logical_coordinates() returns the mesh coordinates of this node
#include "generic_includes.h"

static int squaresize;           // Dimensions of hypercubes
static int nsquares;             // Number of hypercubes
static int machine_coordinates;  // Logical machine coordinates

int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
#define MAXPRIMES (sizeof(prime) / sizeof(int))
// -----------------------------------------------------------------



// -----------------------------------------------------------------
static void setup_hyper_prime() {
  int i, k;

  node0_printf("hyper_prime\n");

  // Figure out dimensions of rectangle
  squaresize = nt;
  nsquares = 1;

  i = 1;  // Current number of hypercubes
  while (i < numnodes()) {
    // Figure out which prime to divide by starting with largest
    k = MAXPRIMES - 1;
    while ((numnodes() / i) % prime[k] != 0 && k > 0)
      --k;

    // This can fail if I run out of prime factors
    if (k < 0) {
      node0_printf("ERROR: Not enough factors of %d to lay out lattice\n",
                   prime[k]);
      g_sync();
      terminate(1);
    }

    // Do the surgery
    i *= prime[k];
    squaresize /= prime[k];
    nsquares *= prime[k];
  }

  // Dividing odd-nt lattice among multiple nodes can give squaresize=0
  // Check for that and exit gracefully if encountered
  if (squaresize < 1) {
    node0_printf("ERROR: Can't lay out lattice ");
    node0_printf("(bad squaresize in layout_hyper_prime)\n");
    g_sync();
    terminate(1);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_layout() {
  int k = mynode();

  node0_printf("LAYOUT = Hypercubes, options = ");
  setup_hyper_prime();

  // Compute machine coordinates
  machine_coordinates = k % squaresize;

  // Number of sites on node
  sites_on_node = squaresize;

  // Need even number of sites per hypercube
  if (sites_on_node % 2 != 0) {
    node0_printf("ERROR: Can't lay out lattice ");
    node0_printf("so that all nodes have an even number of sites\n");
    g_sync();
    terminate(1);
  }

  node0_printf("ON EACH NODE %d\n", squaresize);

  even_sites_on_node = sites_on_node / 2;
  odd_sites_on_node = even_sites_on_node;
}

int node_number(int t) {
  register int i;
  t /= squaresize;
  i = t;
  return i;
}

int node_index(int t) {
  register int i, tr;
  tr = t % squaresize;
  i = tr;
  if (t % 2 == 0)   // Even site
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
  return &nsquares;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Coordinates simulate a mesh architecture
// Must correspond to the node_number result
const int *get_logical_coordinate() {
  return &machine_coordinates;
}
// -----------------------------------------------------------------
