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
double update_gauge_step(Real eps);
void update_h(Real eps);
void update_u(Real eps);

void gauge_action(double *result);
void udadu_mu_nu(field_offset lsrc, field_offset rsrc, field_offset mat,
                 int mu, int nu, int parity);
void udadu_mat_mu_nu(field_offset matsrc, field_offset matdest,
                     int mu, int nu);
void chain_rule(su3_matrix_f *sigmaf, su3_matrix *sigma,
                su3_matrix_f *gaugelinkf);
void apply_bc(su3_matrix *sigma, int dir, int t);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Susy routines
// Lots of things to initialize and set up
#ifdef CATTERALL_ALG
void initialize_TF();
void initialize_p_F();
#endif
void compute_Fmunu();
void compute_DmuUmu();
void setup_lambda();
void setup_PtoP();
void setup_FQ();
void setup_offset();
void setup_qclosed_offset();
void setup_rhmc();
void fermion_rep();
void make_fermion_rep_matrix(su3_matrix_f *a, su3_matrix *b);

// Gaussian random source
int grsource(Twist_Fermion *source);

// Action calculations
double d_action(Twist_Fermion *source, Twist_Fermion **sol);
double d_gauge_action();
double d_fermion_action();

double gauge_force(Real eps);
double fermion_force(Real eps, Twist_Fermion *source, Twist_Fermion **psim);
#ifdef DET
double det_force(Real eps);
#endif

// Reproduces void Dplus(Link_Field *src, Plaq_Field *dest);
void Dplus(su3_vector *src[NUMLINK], su3_vector *dest[NUMLINK][NUMLINK]);

// Reproduces void Dminus(Plaq_Field *src, Link_Field *dest);
void Dminus(su3_vector *src[NUMLINK][NUMLINK], su3_vector *dest[NUMLINK]);

// Reproduces void DbplusPtoP(Plaq_Field *src, Plaq_Field *dest);
void DbplusPtoP(su3_vector *src[NUMLINK][NUMLINK],
                su3_vector *dest[NUMLINK][NUMLINK]);

// Reproduces void DbplusStoL(Site_Field * src, Link_Field *dest);
void DbplusStoL(su3_vector * src, su3_vector *dest[NUMLINK]);

// Reproduces void DbminusPtoP(Plaq_Field *src, Plaq_Field *dest);
void DbminusPtoP(su3_vector *src[NUMLINK][NUMLINK],
                 su3_vector *dest[NUMLINK][NUMLINK]);

// Reproduces void DbminusLtoS(Link_Field *src, Site_Field *dest);
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

// More determinant routines
#ifdef DET
void measure_det();
complex find_det(su3_matrix_f *Q);

void monopole();

void invert(su3_matrix_f *in, su3_matrix_f *out);
void d_flavor();    // Uses find_det
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Correlators and Wilson loops
#ifdef CORR
void d_correlator();    // Konishi and SUGRA correlators
#endif
#ifdef BILIN            // vevs to explore susy breaking
int d_bilinear();       // Fermion bilinear
int d_susyTrans();      // Supersymmetry transformation
#endif
#ifdef PL_CORR
void ploop_c();         // Polyakov loop correlator
void print_var3(char *label);
#endif
#ifdef WLOOP
void hvy_pot();
void polar(su3_matrix_f *a, su3_matrix_f *b);
void hvy_pot_polar();

void path(int *dir, int *sign, int length);
void hvy_pot_loop();
void hvy_pot_polar_loop();

#ifdef DET
// Modified Wilson loops
void path2(int *dir, int *sign, int *kind, int length);
void hvy_pot_rest_loop();    // Uses find_det
#endif
#endif

// Stout smearing
#ifdef STOUT
void block_stout();
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Eigenvalue routines
#ifdef EIG
#include "primme.h"
int make_evs(int Nvec, Twist_Fermion **eigVec, double *eigVal, int flag);
void check_Dmat(int Nvec, Twist_Fermion **eigVec);

// LAPACK is only used to diagonalize <psi_j | D | psi_i>
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
