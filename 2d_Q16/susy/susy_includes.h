// -----------------------------------------------------------------
// Include files for supersymmetric evolution
#include "../include/config.h"  // Keep this first
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>             // For setup.c, gauge_info.c
#include "../include/complex.h"
#include "../include/susy.h"
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
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Susy routines
// Lots of things to initialize and set up
void setup_lambda();
void setup_offset();
void setup_rhmc();

// Helper routines for action and force computations
void compute_plaqdet();
void compute_Uinv();
void compute_DmuUmu();
//--edited 13-9-23-----
void compute_Dmuphi();
void compute_Dmuvarphi();
void compute_varphi_phi_commutator();
void compute_phi_commutator();
void compute_varphi_commutator();
//void compute_Trm1();
//void compute_Trm2();

//---end----
void compute_Fmunu();

// Gaussian random momentum matrices and pseudofermions
void ranmom();
void ranmom_phi(); //--edited---
int grsource(Twist_Fermion *source);

// Basic observables
// #define PLAQ_DIST in local_plaquette to print all plaquettes in serial
void plaquette(double *ss_plaq, double *st_plaq);
double local_plaquette(double *ss_plaq, double *st_plaq); // Return max plaq
complex ploop(int dir, int project, double *plpMod);

// Scalar eigenvalues: averages, extrema and width
// #define SCALAR_EIG_DIST to print all eigenvalues in serial
void scalar_eig(int project, double *ave_eigs, double *eig_widths,
                double *min_eigs, double *max_eigs);

// Action routines
//---edited----------
double action(Twist_Fermion **source, Twist_Fermion ***sol);
double gauge_action(int do_det);
double sa(int do_det);
double bmass_action();

// Force routines
//---edited-----
double gauge_force(Real eps);
double scalar_force(Real eps); // edited
double fermion_force(Real eps, Twist_Fermion *source, Twist_Fermion **psim);
double det_force(Real eps);

// Fermion matrix--vector operators (D & D^2) and multi-mass CG
void fermion_op(Twist_Fermion *src, Twist_Fermion *dest, int sign);
void DSq(Twist_Fermion *src, Twist_Fermion *dest);
int congrad_multi(Twist_Fermion *src, Twist_Fermion **psim,
                  int MaxCG, Real RsdCG, Real *size_r);

// Compute average Tr[Udag U] / N_c
double link_trace(double *linktr, double *linktr_width,
                  double *dets, double *det_ave, double *det_width);

// Basic Twist_Fermion utilities
// May eventually move to libraries
void dump_TF(Twist_Fermion *in);
void copy_TF(Twist_Fermion *src, Twist_Fermion *dest);
void clear_TF(Twist_Fermion *dest);
Real magsq_TF(Twist_Fermion *in);
complex TF_dot(Twist_Fermion *a, Twist_Fermion *b);
void TF_rdot_sum(Twist_Fermion *a, Twist_Fermion *b, Real *c);
void sum_TF(Twist_Fermion *b, Twist_Fermion *c);
void scalar_mult_sum_TF(Twist_Fermion *b, Real s, Twist_Fermion *c);
void scalar_mult_add_TF(Twist_Fermion *a, Twist_Fermion *b, Real s,
                        Twist_Fermion *c);
void scalar_mult_dif_TF(Twist_Fermion *b, Real s, Twist_Fermion *c);
void scalar_mult_TF(Twist_Fermion *src, Real s, Twist_Fermion *dest);

// Other routines in library_util.c that loop over all sites
void gauge_field_copy(field_offset src, field_offset dest);
void scalar_field_copy(field_offset src, field_offset dest); //--edited------
void shiftmat(matrix *dat, matrix *temp, int dir);

// Random gauge transformation for testing gauge invariance
void random_gauge_trans(Twist_Fermion *TF);

// Determinant-related routines
void measure_det();
void widths();        // Widths of plaquette and det distributions
complex find_det(matrix *Q);
void det_project(matrix *in, matrix *out);

// Use LAPACK in determinant and matrix calculations
// Compute LU decomposition of a complex matrix
// http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/zgetrf.html
// First and second arguments are the dimensions of the matrix
// Third argument is the LU-decomposed matrix
// Fourth argument is the
// Fifth argument is the LU decomposition pivot matrix
// Final argument reports success or information about failure
void zgetrf_(int *N1, int *N2, double *store, int *lda, int *ipiv, int *stat);

// Invert a complex matrix given its LU decomposition
// http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/zgetri.html
// First four and last arguments are defined above
// Fifth argument is real workspace of size given by the sixth argument
void zgetri_(int *N, double *store, int *lda, int *ipiv,
             double *work, int *Nwork, int* stat);

// Matrix inverse via LAPACK
void invert(matrix *in, matrix *out);

// Modified Wilson loops use invert and path
void path(int *dir, int *sign, int length);
void rsymm();
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// More measurements
#ifdef CORR
// Konishi and SUGRA correlators
void compute_Ba();
void konishi();       // Operators averaged over each timeslice

// Map (x, t) to scalar displacements r
Real r_map(int x_in, int t_in);
void correlator_r();            // Functions of r
#endif

#ifdef BILIN
// Ward identity involving eta.psi_a fermion bilinear
int bilinearWard();
#endif

#ifdef PL_CORR
// Polyakov loop correlator -- NOT CURRENTLY IN USE
void ploop_c();         // Computes Polyakov loop at each spatial site
#endif

#ifdef WLOOP
// Wilson loops -- including determinant division and polar projection
// These look at correlators of products of temporal links,
// which requires gauge fixing
void hvy_pot(int do_det);
void hvy_pot_polar();

// These construct explicit paths along lattice principal axes, for checking
void hvy_pot_loop(int do_det);
void hvy_pot_polar_loop();
#endif

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
void polar(matrix *in, matrix *u, matrix *P);
void matrix_log(matrix *in, matrix *out);

// Monopole computation uses find_det
void monopole();

#ifdef SMEAR
// Stout and APE-like smearing, the latter with no final SU(N) projections
// APE-like smearing optionally builds staples from projected links
void exp_mult();
void stout_smear(int Nsmear, double alpha);
void APE_smear(int Nsmear, double alpha, int project);
#endif

#ifdef MCRG
void block_mcrg(int bl);
void blocked_plaq(int Nsmear, int bl);    // Also monitors det and widths
void blocked_local_plaq(int Nsmear, int bl);    // Print max plaq
void blocked_ops(int Nsmear, int bl);
void blocked_ploop(int Nsmear, int bl);
void blocked_rsymm(int Nsmear, int bl);
#ifdef SMEAR            // Blocked stout and APE-like smearing, as above
void blocked_stout(int Nsmear, double alpha, int block);
void blocked_APE(int Nsmear, double alpha, int project, int block);
#endif
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Eigenvalue routines
#ifdef EIG
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
void phase();
#endif
// -----------------------------------------------------------------
