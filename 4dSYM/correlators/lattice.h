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

EXTERN double g_ssplaq, g_stplaq;   // Global plaqs for I/O
EXTERN double_complex linktrsum;
EXTERN u_int32type nersc_checksum;
EXTERN char startfile[MAXFILENAME];
EXTERN int startflag;     // Beginning lattice: CONTINUE, RELOAD, FRESH

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

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
#define N_POINTERS 1
EXTERN char **gen_pt[N_POINTERS];

#define N_B 2
#define N_K 3    // N_B * (N_B + 1) / 2
// Multiple scalar fields and their bilinear traces
EXTERN matrix *Ba[N_B][NUMLINK];
EXTERN double *traceBB[N_K][NUMLINK][NUMLINK];

// Ensemble averages and volume averages for subtracting
EXTERN double vevK[N_K], vevS[N_K], volK[N_K], volS[N_K];

// Structs for operators and correlators with either subtraction
typedef struct {
  double OK[N_K][2];
  double OS[N_K][2];
} Kops;
typedef struct {
  double C[N_K][N_K][2];
} Kcorrs;
EXTERN Kops *tempops, *tempops2;

#ifdef SMEAR
// APE or stout smearing stuff
EXTERN int smearflag;                 // NONE, STOUT, APE
EXTERN int Nsmear;
EXTERN double alpha;
EXTERN anti_hermitmat *Q[NUMLINK];    // To be exponentiated
#endif
#endif // _LATTICE_H
// -----------------------------------------------------------------
