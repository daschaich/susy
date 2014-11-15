// -----------------------------------------------------------------
// Include files for supersymmetric evolution
#include "../include/config.h"  // Keep this first
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>             // For print_var.c, setup.c, gauge_info.c
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);
int update();
void update_h(Real eps);
void update_u(Real eps);
void gauge_action(double *result);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Susy routines
// Lots of things to initialize and set up
void compute_Fmunu();
void compute_DmuUmu();
void setup_lambda();
void setup_PtoP();
void setup_FQ();
void setup_offset();
void setup_qclosed_offset();
void setup_rhmc();
void fermion_rep();

// Gaussian random source
int grsource(Twist_Fermion *source);

// Action calculations
double d_action(Twist_Fermion *source, Twist_Fermion **sol);
double d_gauge_action();
double d_fermion_action();

double gauge_force(Real eps);
double fermion_force(Real eps, Twist_Fermion *source, Twist_Fermion **psim);
double det_force(Real eps);

// Link-to-plaq term in action
void Dplus(su3_vector *src[NUMLINK], su3_vector *dest[NPLAQ]);

// Plaq-to-link term in action
void Dminus(su3_vector *src[NPLAQ], su3_vector *dest[NUMLINK]);

// First plaq-to-plaq term in action
void DbplusPtoP(su3_vector *src[NPLAQ], su3_vector *dest[NPLAQ]);

// Site-to-link term in action
void DbplusStoL(su3_vector *src, su3_vector *dest[NUMLINK]);

// Second plaq-to-plaq term in action
void DbminusPtoP(su3_vector *src[NPLAQ], su3_vector *dest[NPLAQ]);

// Link-to-site term in action
void DbminusLtoS(su3_vector *src[NUMLINK], su3_vector *dest);

// Fermion matrix--vector operators (D & D^2) and multi-mass CG
void fermion_op(Twist_Fermion *src, Twist_Fermion *dest, int sign);
void hdelta0_field(Twist_Fermion *src, Twist_Fermion *dest);
int congrad_multi_field(Twist_Fermion *src, Twist_Fermion **psim,
                        int MaxCG, Real RsdCG, Real *size_r);

// Compute average link Tr[Udag U] / N_c
void d_link();
void d_link_frep();

Real order(int i, int j, int k, int l, int m);
void epsilon();

// Basic Twist_Fermion and gauge field manipulations
// May eventually move to libraries
void dump_TF(Twist_Fermion *vec);
void conjTF(Twist_Fermion *src, Twist_Fermion *dest);
void copy_TF(Twist_Fermion *src, Twist_Fermion *dest);
void clear_TF(Twist_Fermion *dest);
Real magsq_TF(Twist_Fermion *vec);
complex TF_dot(Twist_Fermion *a, Twist_Fermion *b);
void scalar_mult_add_TF(Twist_Fermion *src1, Twist_Fermion *src2, Real s,
                        Twist_Fermion *dest);
void scalar_mult_TF(Twist_Fermion *src, Real s, Twist_Fermion *dest);
void gauge_field_copy_f(field_offset src, field_offset dest);
void shiftmat(field_offset src, field_offset dest, int dir);

// Determinant-related routines
void measure_det();
complex find_det(su3_matrix_f *Q);

// Adjugate matrix needed by det_force
void adjugate(su3_matrix_f *in, su3_matrix_f *out);

// Matrix invert is just adjugate divided by determinant
void invert(su3_matrix_f *in, su3_matrix_f *out);

// Modified Wilson loops use invert and path
void path(int *dir, int *sign, int length);
void rsymm();
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// More measurements
#ifdef CORR
// Konishi and SUGRA correlators
void setup_P();
void compute_Bmu();
void d_correlator();    // Projected to zero spatial momentum
void d_correlator_r();  // Functions of (x, y, z, t)
#endif
#ifdef BILIN
// vevs to explore susy breaking
int d_bilinear();       // Fermion bilinear
int d_susyTrans();      // Supersymmetry transformation
#endif
#ifdef PL_CORR
// Polyakov loop correlator -- NOT CURRENTLY IN USE
void ploop_c();
void print_var3(char *label);
#endif
#ifdef WLOOP
// Wilson loops -- including determinant division and polar projection
// These look at correlators of products of temporal links,
// which requires gauge fixing
void hvy_pot();
void hvy_pot_polar();

// These construct explicit paths along lattice principal axes, for checking
void hvy_pot_loop();
void hvy_pot_polar_loop();

// Use LAPACK in the polar projection
// http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/zheev.html
// First argument turns on eigenvector computations
// Second argument chooses between storing upper or lower triangle
// Third and fifth arguments are the dimensions of the matrix
// Fourth argument is that matrix, overwritten by the eigenvectors
// Sixth argument holds the computed eigenvalues
// Seventh argument is real workspace of size given by the eighth argument
// Ninth argument is real workspace of size 3 * NCOL - 2
// Final argument reports success or information about failure
void zheev_(char *doV, char *uplo, int *N1, double *store, int *N2,
            double *eigs, double *work, int *Nwork, double *Rwork, int *stat);
void polar(su3_matrix_f *a, su3_matrix_f *b);
#endif

// Monopole computation uses find_det
void monopole();

#ifdef MCRG
void block_mcrg(int bl);
void blocked_ops(int bl);
void blocked_ploop(int bl);
void blocked_rsymm(int bl);
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Eigenvalue routines
#ifdef EIG
#include "primme.h"
int make_evs(int Nvec, Twist_Fermion **eigVec, double *eigVal, int flag);
void check_Dmat(int Nvec, Twist_Fermion **eigVec);

// Use LAPACK to diagonalize <psi_j | D | psi_i>
// on the subspace of Ddag.D eigenvalues psi
// http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/zgeev.html
// First two arguments turn off eigenvector computations
// Third and fifth arguments are the dimensions of the matrix
// Fourth argument is that matrix, which will be overwritten
// Sixth argument holds the computed eigenvalues
// Seventh and ninth arguments are eigenvectors
// Eighth and tenth arguments are the dimensions of the eigenvectors
// Eleventh argument is real workspace, of size given by the twelfth argument
// Thirteenth argument is real workspace, of size given by the third argument
// Final argument reports success or information about failure
void zgeev_(char *doL, char *doR, int *N1, double *store, int *N2, double *eigs,
            double *dumL, int *NL, double *dumR, int *NR,
            double *work, int *Nwork, double *Rwork, int *stat);
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Pfaffian phase
#ifdef PHASE
void d_phase();
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Stochastic mode number
#ifdef MODE
void coefficients();
void step(Twist_Fermion *src, Twist_Fermion *res);
#endif
// -----------------------------------------------------------------
