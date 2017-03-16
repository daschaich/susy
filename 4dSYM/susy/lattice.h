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
// Twist_Fermion struct
typedef struct {
  matrix Fsite;
  matrix Flink[NUMLINK];
  matrix Fplaq[NPLAQ];
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

  matrix link[NUMLINK];       // Gauge links

#ifdef HMC_ALGORITHM
  matrix old_link[NUMLINK];   // For accept/reject
#endif

  // Momentum matrices in each direction are just U(N) matrices
  // as opposed to anti-hermitian matrices
  matrix mom[NUMLINK], f_U[NUMLINK];        // Force matrices

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

EXTERN int nx, ny, nz, nt;  // Lattice dimensions
EXTERN int PBC;             // Temporal fermion boundary condition
EXTERN int volume;          // Volume of lattice
EXTERN int iseed;           // Random number seed
EXTERN int warms, trajecs, niter, propinterval;
EXTERN Real traj_length;

// U(N) generators, epsilon tensor
EXTERN matrix Lambda[DIMF];
EXTERN Real perm[NUMLINK][NUMLINK][NUMLINK][NUMLINK][NUMLINK];

// Translate (mu, nu) to linear index of anti-symmetric matrix
EXTERN int plaq_index[NUMLINK][NUMLINK];

EXTERN Real rsqmin, lambda, kappa, bmass, fmass, kappa_u1, G, B;
EXTERN int doG, doB;
EXTERN double g_ssplaq, g_stplaq;   // Global plaqs for I/O
EXTERN double_complex linktrsum;
EXTERN u_int32type nersc_checksum;
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN int startflag;     // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int fixflag;       // Gauge fixing: COULOMB_GAUGE_FIX, NO_GAUGE_FIX
EXTERN int saveflag;      // 1 if we will save the lattice;
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
EXTERN matrix *site_src, *link_src[NUMLINK], *plaq_src[NPLAQ];
EXTERN matrix *site_dest, *link_dest[NUMLINK], *plaq_dest[NPLAQ];

// For convenience in calculating action and force
// May be wasteful of space
EXTERN Real one_ov_N;
EXTERN complex minus1, *tr_eta;
EXTERN complex *tr_dest, *Tr_Uinv[NUMLINK], *plaqdet[NUMLINK][NUMLINK];
EXTERN complex *ZWstar[NUMLINK][NUMLINK], *tempdet[NUMLINK][NUMLINK];
#ifndef LINEAR_DET
EXTERN complex *tempZW[NUMLINK][NUMLINK];
#endif
EXTERN matrix *DmuUmu, *Fmunu[NPLAQ];
EXTERN matrix *Uinv[NUMLINK], *Udag_inv[NUMLINK], *UpsiU[NUMLINK];

// CG Twist_Fermions
EXTERN Twist_Fermion *mpm, *pm0, *rm;

// Temporary matrices and Twist_Fermion
EXTERN matrix *tempmat, *tempmat2, *staple;
EXTERN Twist_Fermion *tempTF;

// Arrays to be used by LAPACK in determinant.c
EXTERN int *ipiv;
EXTERN double *store, *work;

// Allocate some more arrays to be used by LAPACK in unit.c
EXTERN double *Rwork, *eigs;

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
#define N_B 2
#define N_K 3    // N_B * (N_B + 1) / 2
// Multiple scalar fields and their bilinear traces
EXTERN matrix *Ba[N_B][NUMLINK];
EXTERN double *traceBB[N_K][NUMLINK][NUMLINK];

// Structs for operators and correlators
typedef struct {
  double OK[N_K];
  double OS[N_K];
} Kops;
typedef struct {
  double C[N_K][N_K];
} Kcorrs;
EXTERN Kops *tempops, *tempops2;
#endif

#ifdef SMEAR
// APE or stout smearing stuff
EXTERN int smearflag;                 // NONE, STOUT, APE
EXTERN int Nsmear;
EXTERN double alpha;
EXTERN anti_hermitmat *Q[NUMLINK];    // To be exponentiated
#endif

#ifdef BILIN
EXTERN int nsrc;
#endif

#ifdef EIG
// Eigenvalue stuff
EXTERN int Nvec;
EXTERN double *eigVal;
EXTERN Twist_Fermion *src, *res;      // For av_ov matvec
EXTERN Twist_Fermion **eigVec;
EXTERN Real eig_tol;          // Tolerance for the eigenvalue computation
EXTERN int maxIter;           // Maximum iterations
#endif

#if defined(CHEB) || defined(MODE)
// Z2 random source and stochastic source stuff for both
// Chebyshev spectral density and Giusti--Luescher mode number
EXTERN int Nstoch;                    // Number of stochastic sources
EXTERN Real sqrt1_ov_Nm1;             // 1 / sqrt(N - 1) for averaging
EXTERN Twist_Fermion *z_rand;         // Z2 random Twist_Fermion
#endif

#ifdef CHEB
// Chebyshev spectral density stuff
EXTERN int cheb_order;                // Number of coefficients to compute
EXTERN Real *cheb_coeff, *cheb_err;   // Results for coefficients
EXTERN Real lambda_min, lambda_max;   // Bounds on spectral range
#endif

#ifdef MODE
// Giusti--Luescher mode number and step function stuff
EXTERN int step_order;            // Selects between options in mode_coeffs.c
EXTERN int numOmega;              // Number of Omega at which to evaluate nu
EXTERN Real *Omega;               // List of Omega at which to evaluate nu
EXTERN Real OmStar;               // Scaled Omega / star given to CG

EXTERN Real step_eps;             // Epsilon of step function approximation
EXTERN Real delta;                // Just recorded in the output
EXTERN Real starSq, star;         // Ratio (Omega / Omega_*)^2 and its sqrt
EXTERN double *step_coeff;        // Options hard-coded in mode_coeffs.c
EXTERN double *mode, *err;        // Results for mode number

// Z2 stochastic sources and temporary storage
EXTERN Twist_Fermion **source, *XPXSq, *hX, *dest;
EXTERN Twist_Fermion *bj, *bjp1;  // More temporaries for Clenshaw algorithm
#endif

#ifdef PHASE
// Pfaffian phase stuff
EXTERN long int Nmatvecs;           // For timing/counting
EXTERN int ckpt_load, ckpt_save;    // For checkpointing
EXTERN Twist_Fermion *src, *res;    // For fieldwise matvec
#endif

#endif // _LATTICE_H
// -----------------------------------------------------------------
