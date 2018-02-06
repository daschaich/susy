// -----------------------------------------------------------------
// Allocate space for lattice fields in the site struct
// Fill in coordinates, parity, index
// Allocate gen_pt pointers for gather results
// Initialize site-based random number generator if SITERAND defined
#include "generic_includes.h"
#include <defines.h>                  // For SITERAND

void make_lattice() {
  register int i;
  int t;

  // Allocate space for lattice
  node0_printf("Mallocing %.1f MBytes per core for lattice\n",
               (double)sites_on_node * sizeof(site) / 1e6);
  lattice = malloc(sizeof *lattice * sites_on_node);
  if (lattice == NULL) {
    printf("node%d: no room for lattice\n", this_node);
    terminate(1);
  }

  // Allocate address vectors
  for (i = 0; i < N_POINTERS; i++) {
    gen_pt[i] = malloc(sizeof(char*) * sites_on_node);
    if (gen_pt[i] == NULL) {
      printf("node%d: no room for pointer array\n", this_node);
      terminate(1);
    }
  }

  // Fill in parity, coordinates and index
  for (t = 0; t < nt; t++) {
    if (node_number(t) == mynode()) {
      i = node_index(t);
      lattice[i].t = t;
      lattice[i].index = t;
      if (t % 2 == 0)
        lattice[i].parity = EVEN;
      else
        lattice[i].parity = ODD;
#ifdef SITERAND
      initialize_prn(&(lattice[i].site_prn), iseed, lattice[i].index);
#endif
    }
  }
}

void free_lattice() {
  int i;

  for (i = 0; i < N_POINTERS; i++)
    free(gen_pt[i]);

  free(lattice);
}
// -----------------------------------------------------------------
