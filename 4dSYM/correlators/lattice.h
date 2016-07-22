// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/macros.h"    // For MAXFILENAME
#include "../include/io_lat.h"    // For gauge_file
#include "../include/susy.h"
#include "../include/random.h"    // For double_prn
#include "../include/dirs.h"      // For NDIMS
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, z, t;   // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

  matrix link[NUMLINK];       // Gauge links

#ifdef SMEAR          // For directional staple routine
  matrix mom[NUMLINK];
#endif
} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

// Global variables
EXTERN int nx, ny, nz, nt;  // Lattice dimensions
EXTERN int volume;          // Volume of lattice
EXTERN int Nblock;          // Number of blocks used in analyses
EXTERN int Nmeas;           // Number of lattices to load per block
EXTERN int tot_meas;        // Nblock * Nmeas

EXTERN double g_ssplaq, g_stplaq;         // Global plaqs for I/O
EXTERN double_complex linktrsum;
EXTERN u_int32type nersc_checksum;
EXTERN char cfg[MAX_CFG][MAXFILENAME];    // Lattice file names

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

// Stuff for gathers -- both forwards and backwards
EXTERN int offset[NUMLINK][NDIMS];    // Path along each link
EXTERN int goffset[2 * NUMLINK];

// Temporary matrices
EXTERN matrix *tempmat, *tempmat2, *staple;

// Arrays to be used by LAPACK in unit.c
EXTERN int *ipiv;
EXTERN double *store, *work, *Rwork, *eigs;

EXTERN gauge_file *startlat_p;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
// Need two for plaquettes, three for smearing
#define N_POINTERS 3
EXTERN char **gen_pt[N_POINTERS];

// Scalar fields and their bilinear traces
// Just need one set of N_B / N_K for current measurement
#define N_B 1
#define N_K 1     // N_B * (N_B + 1) / 2
EXTERN matrix *Ba[N_B][NUMLINK];
EXTERN double *traceBB[N_K][NUMLINK][NUMLINK];
EXTERN Real one_ov_N;                 // For trace subtraction

// Structs for operators and correlators
// Will need Nblock of each
// Need to correlate every operator with every other
#define N_op 2    // 2 N_K
typedef struct {
  double O[N_op];
} Kops;
typedef struct {
  double C[N_op][N_op];
} Kcorrs;
EXTERN Kops *tempops, *tempops2;      // For shifting
EXTERN Kops **ops;                    // Will be Nblock of these

// Correlator stuff that only depends on volume
EXTERN Real *lookup;                  // Unsorted list of scalar distances r
EXTERN Real *norm;                    // By number of four-vectors for each r
EXTERN Real MAX_r, total_r;

#ifdef SMEAR
// APE or stout smearing stuff
EXTERN int smearflag;                 // NONE, STOUT, APE
EXTERN int Nsmear;
EXTERN double alpha;
EXTERN anti_hermitmat *Q[NUMLINK];    // To be exponentiated
#endif
#endif // _LATTICE_H
// -----------------------------------------------------------------
