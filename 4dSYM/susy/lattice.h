// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/macros.h"    // For MAXFILENAME
#include "../include/io_lat.h"    // For gauge_file
#include "../include/su3.h"
#include "../include/random.h"    // For double_prn
#include "../include/dirs.h"      // For NDIMS
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Twist_Fermion struct
typedef struct {
  su3_vector Fsite;
  su3_vector Flink[NUMLINK];
  su3_vector Fplaq[NPLAQ];
} Twist_Fermion;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, z, t;   // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

#ifdef SITERAND
  // The state information for a random number generator
  double_prn site_prn;
#endif

  su3_matrix_f linkf[NUMLINK];    // Gauge links
  su3_matrix link[NUMLINK];       // Adjoint links

#ifdef HMC_ALGORITHM
  su3_matrix_f old_linkf[NUMLINK];  // For accept/reject
#endif

  // Momentum matrices in each direction are just U(N) matrices
  // as opposed to anti-hermitian matrices
  su3_matrix_f mom[NUMLINK];
  su3_matrix_f f_U[NUMLINK];        // Force matrices

  // Boundary conditions -- many unused
  Real bc1[2 * NUMLINK], bc2[2 * NUMLINK][2 * NUMLINK];
  Real bc3[2 * NUMLINK][2 * NUMLINK][2 * NUMLINK];

#ifdef PL_CORR
  complex print_var, ploop_corr, fft1, fft2;
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
EXTERN int PBC;             // Temporal fermion boundary condition
EXTERN int volume;          // Volume of lattice
EXTERN int iseed;           // Random number seed
EXTERN int warms, trajecs, niter, propinterval;
EXTERN Real traj_length;

// U(N) generators, epsilon tensor
EXTERN su3_matrix_f Lambda[DIMF], Lambda_prod[DIMF][DIMF];
EXTERN Real perm[NUMLINK][NUMLINK][NUMLINK][NUMLINK][NUMLINK];

// Submatrix to convert from 5- to 4-component notation
// Fifth row is just all 1 / sqrt(5)
EXTERN Real P[NDIMS][NUMLINK];

// Translate (mu, nu) to linear index of anti-symmetric matrix
EXTERN int plaq_index[NUMLINK][NUMLINK];

EXTERN Real rsqmin, rsqprop;
EXTERN Real lambda, kappa, bmass, fmass, kappa_u1, G, B;
EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktrsum;
EXTERN u_int32type nersc_checksum;
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN int startflag; // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int fixflag;   // Gauge fixing: COULOMB_GAUGE_FIX, NO_GAUGE_FIX
EXTERN int saveflag;  // 1 if we will save the lattice;
EXTERN int total_iters;   // To be incremented by the multi-mass CG

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

// Stuff for multi-mass CG and RHMC
EXTERN int nsteps[2];           // Fermion and gauge steps
EXTERN Real ampdeg, *amp, *shift;
EXTERN Real ampdeg4, *amp4, *shift4;
EXTERN Real ampdeg8, *amp8, *shift8;
EXTERN int Nroot, Norder;
EXTERN Real gnorm, *fnorm, max_gf, *max_ff;

// Each node maintains a structure with the pseudorandom number
// generator state
EXTERN double_prn node_prn;

// Stuff for derivative and link terms
EXTERN int offset[NUMLINK][NDIMS];    // Path along each link
EXTERN int label[NUMLINK];
EXTERN int q_off_max, q_offset[NQLINK][4];

// Stuff for gathers -- both forwards and backwards
EXTERN int goffset[2 * NUMLINK];
EXTERN int gq_offset[2 * NQLINK];

EXTERN int DbmP_d1[NTERMS], DbmP_d2[NTERMS];
EXTERN int DbpP_d1[NTERMS], DbpP_d2[NTERMS];
EXTERN int DbplusPtoP_lookup[NTERMS][NUMLINK];
EXTERN int DbminusPtoP_lookup[NTERMS][NUMLINK];

EXTERN int F1Q_d1[NTERMS], F1Q_d2[NTERMS];
EXTERN int F2Q_d1[NTERMS], F2Q_d2[NTERMS];
EXTERN int FQ_lookup[NTERMS][NUMLINK];

// Persistent site, link and plaq fermions for matrix--vector operation
// Used in fermion_op and assemble_fermion_force
EXTERN su3_vector *site_src, *link_src[NUMLINK], *plaq_src[NPLAQ];
EXTERN su3_vector *site_dest, *link_dest[NUMLINK], *plaq_dest[NPLAQ];

// For convenience in calculating action and force
// May be wasteful of space
EXTERN complex minus1;
EXTERN complex *tr_dest, *Tr_Uinv[NUMLINK], *plaqdet[NUMLINK][NUMLINK];
EXTERN complex *ZWstar[NUMLINK][NUMLINK];
EXTERN su3_matrix_f *DmuUmu, *Fmunu[NPLAQ];
EXTERN su3_matrix_f *Uinv[NUMLINK], *Udag_inv[NUMLINK], *UpsiU[NUMLINK];

// Temporary vectors, matrices and Twist_Fermion
EXTERN su3_vector *tempvec[NUMLINK];
EXTERN su3_matrix_f *tempmat1, *tempmat2, *staple;
EXTERN Twist_Fermion *tempTF;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
// Need 10 for gauge-fixing, 9 for Q-invariant determinant force
#define N_POINTERS 10
EXTERN char **gen_pt[N_POINTERS];

#ifdef CORR
EXTERN su3_matrix_f *Ba[NUMLINK];         // Scalar field interpolating op
EXTERN Real *traceBB[NUMLINK][NUMLINK];   // Tr[B_a B_b]
EXTERN double vevK;                       // Konishi vacuum subtraction
#endif

#ifdef BILIN
EXTERN int nsrc;
#endif

#ifdef SMEAR
// APE or stout smearing stuff
EXTERN int Nsmear;
EXTERN double alpha;
#endif
// This is needed for smearing a `hot-start' random configuration
EXTERN anti_hermitmat *Q[NUMLINK];    // To be exponentiated

#ifdef EIG
// Eigenvalue stuff
EXTERN int Nvec;
EXTERN double *eigVal;
EXTERN Twist_Fermion *src, *res;      // For av_ov matvec
EXTERN Twist_Fermion **eigVec;
EXTERN Real eig_tol;          // Tolerance for the eigenvalue computation
EXTERN int maxIter;           // Maximum iterations
#endif

#ifdef PHASE
// Pfaffian phase stuff
EXTERN Twist_Fermion *src, *res;    // For fieldwise matvec
EXTERN int Nmatvecs;                // For timing/counting
EXTERN int ckpt_load, ckpt_save;    // For checkpointing
#endif

#ifdef MODE
// Mode number (and associated step function) stuff
EXTERN int Npts;
EXTERN Real M, spacing;
EXTERN Twist_Fermion **source;

EXTERN int step_order;  // Selects between options hard-coded in coeffs.c
EXTERN double eps;
EXTERN double delta;    // Unused, but may be useful to record in the output
EXTERN double starSq, star;     // Ratio (Omega / Omega_*)^2 and its sqrt
EXTERN double *coeffs;
#endif

#endif // _LATTICE_H
// -----------------------------------------------------------------
