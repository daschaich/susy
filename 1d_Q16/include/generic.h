// -----------------------------------------------------------------
// Macros and declarations for miscellaneous generic routines
#ifndef _GENERIC_H
#define _GENERIC_H

// Other generic directory declarations are elsewhere:
//   See comdefs.h for communications
//   See io_lat.h for I/O
#include <stdio.h>
#include "../include/int32type.h"
#include "../include/complex.h"
#include "../include/susy.h"
#include "../include/macros.h"
#include "../include/random.h"
#include "../include/file_types.h"
#include "../include/io_lat.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// check_unitarity.c and check_antihermity.c
Real check_unit(matrix *c);
Real check_unitarity();
Real check_ah();
Real check_antihermity();

// io_helpers.c
gauge_file *save_lattice(int flag, char *filename);
gauge_file *reload_lattice(int flag, char *filename);
int ask_starting_lattice(FILE *fp, int prompt, int *flag, char *filename);
int ask_ending_lattice(FILE *fp, int prompt, int *flag, char *filename);
int ask_gauge_fix(FILE *fp, int prompt, int *flag);
void coldlat();
void randomlat();
void funnylat();
int get_f(FILE *fp, int prompt, char *variable_name_string, Real *value);
int get_i(FILE *fp, int prompt, char *variable_name_string, int *value);
int get_vi(FILE *fp, int prompt, char *variable_name_string,
           int *value, int nvalues);
int get_vf(FILE *fp, int prompt, char *variable_name_string,
           Real *value, int nvalues);
int get_s(FILE *fp, int prompt, char *variable_name_string, char *value);
int get_prompt(FILE *fp, int *value);

// layout_hyper_prime.c
void setup_layout();
int node_number(int t);
int node_index(int t);
size_t num_sites(int node);
const int *get_logical_dimensions();
const int *get_logical_coordinate();

// make_lattice.c
void make_lattice();
void free_lattice();

// nersc_cksum.c
void sum_linktr(double_complex *linktrsum);
u_int32type nersc_cksum();

// remap_stdio_from_args.c
int remap_stdio_from_args(int argc, char *argv[]);

// ranstuff.c
void initialize_prn(double_prn *prn_pt, int seed, int index);
Real myrand(double_prn *prn_pt);

// reunitarize.c and reantihermize.c
int check_deviation();
void reunitarize();
void reantihermize();
// Use LAPACK singular value decomposition for reunitarization
// http://www.netlib.org/lapack/explore-3.1.1-html/zgesvd.f.html
// First and second arguments tell LAPACK to compute all singular values
// Third and fourth arguments are the dimensions of the matrix (both NCOL)
// Fifth argument is the input matrix (lost)
// Sixth argument is the leading dimension (NCOL)
// Seventh argument is the array of singular values (discarded)
// Eight argument is the matrix of left singular vectors (left)
// Ninth argument is the dimension of left (NCOL)
// Tenth argument is the matrix of right singular vectors (right^dag)
// Eleventh argument is the dimension of right^dag (NCOL)
// Twelfth argument is complex workspace of size given by the 13th argument
// Fourteenth argument is real workspace of size 5 * NCOL
// Final argument reports success or information about failure
void zgesvd_(char *A1, char *A2, int *N1, int *N2, double *store,
             int *lda, double *junk, double *left, int *Nl,
             double *right, int *Nr, double *work, int *Nwork,
             double *Rwork, int *stat);
#endif
// -----------------------------------------------------------------
