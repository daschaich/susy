// -----------------------------------------------------------------
// Include files for supersymmetric evolution
#include "../include/config.h"  // Keep this first
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>             // For print_var.c, setup.c, gauge_info.c
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
void setup_offset();

// Basic observables
void plaquette(double *ss_plaq, double *st_plaq);
double local_plaquette(double *ss_plaq, double *st_plaq); // Return max plaq

// Map (x, y, z, t) to scalar displacements r
Real A4map(int x_in, int y_in, int z_in, int t_in);

// Count and return total number of unique scalar distances
int count_points(Real MAX_r);

// Accumulate Konishi and SUGRA operators within current block
void add_konishi();

// Compute correlators C(r) of vacuum-subtracted Konishi and SUGRA operators
// Print both average and standard error
void correlator_r();

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
// -----------------------------------------------------------------
