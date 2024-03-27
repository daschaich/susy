// -----------------------------------------------------------------
// Defines and subroutine declarations for N=4 SYM and its dim'l reductions
#ifndef _SUSY_H
#define _SUSY_H

#include "../include/complex.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// N=4 SYM fermions always in adjoint rep, as are dim-reduced scalars
// This file is only included for d-dimensional systems with d>=2
// The gauge group is U(NCOL) and size of adjoint rep DIMF = NCOL^2
#define NCOL 2
#define DIMF (NCOL * NCOL)
#define N_OFFDIAG (NCOL * (NCOL - 1) / 2)

// Fundamental representation
// will attempt to keep everything in memory as NCOL * NCOLF matrices
// and add comments stating what the actual representation of variables
// this way only a single set of extra functions will be needed
// phi, lambda, zeta_a and xi are NCOL * NCOLF while
// varphi, lambda_a, zeta are NCOLF * NCOL
// NCOLF is the  number of flavours
#define NCOLF 3
#define DIMFUND (NCOL * NCOLF)

typedef struct { fcomplex e[NCOL][NCOL]; } fmatrix;
typedef struct { fcomplex c[NCOL]; } fvector;

// Anti-hermitian matrices for general NCOL
typedef struct {
  fcomplex m[N_OFFDIAG];
  float im_diag[NCOL];
} fanti_hermitmat;

typedef struct { fcomplex e[NCOL][NCOLF]; } ffunmatrix;

typedef struct { dcomplex e[NCOL][NCOL]; } dmatrix;
typedef struct { dcomplex c[NCOL]; } dvector;
typedef struct {
  dcomplex m[N_OFFDIAG];
  double im_diag[NCOL];
} danti_hermitmat;

typedef struct { dcomplex e[NCOL][NCOLF]; } dfunmatrix;

#if (PRECISION == 1)
#define matrix    fmatrix
#define vector    fvector
#define anti_hermitmat  fanti_hermitmat
#define funmatrix ffunmatrix
#else
#define matrix    dmatrix
#define vector    dvector
#define anti_hermitmat  danti_hermitmat
#define funmatrix dfunmatrix
#endif

// Need SU(2) matrices for any U(N), e.g. for gauge hits when gauge-fixing
typedef struct { complex e[2][2]; } su2_matrix;

// Flags for selecting M or M_adjoint
#define PLUS 1
#define MINUS -1
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Subroutine definitions
// Vector operations
// In file clearvec.c
void clearvec(vector *v);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Matrix operations
// In file dumpmat.c
void dumpmat(matrix *a);

// In file clear_mat.c
void clear_mat(matrix *c);

// In file trace.c
complex trace(matrix *a);

// In file realtr.c
Real realtrace_nn(matrix *a, matrix *b);
Real realtrace(matrix *a, matrix *b);

// In file complextr.c
complex complextrace_nn(matrix *a, matrix *b);
complex complextrace_an(matrix *a, matrix *b);
complex complextrace_na(matrix *a, matrix *b);

// In file mat_copy.c
void mat_copy(matrix *a, matrix *b);

// In file adjoint.c
void adjoint(matrix *a, matrix *b);
void neg_adjoint(matrix *a, matrix *b);

// In file addmat.c
void sum_matrix(matrix *b, matrix *c);
void add_matrix(matrix *a, matrix *b, matrix *c);

// In file addamat.c
void sum_adj_matrix(matrix *b, matrix *c);
void add_adj_matrix(matrix *a, matrix *b, matrix *c);

// In file submat.c
void dif_matrix(matrix *b, matrix *c);
void sub_matrix(matrix *a, matrix *b, matrix *c);

// In file subamat.c
void sub_adj_matrix(matrix *a, matrix *b, matrix *c);

// In file s_a_d_mat.c
void scalar_add_diag(matrix *a, Real s);

// In file s_m_mat.c
void scalar_mult_matrix(matrix *src, Real s, matrix *c);

// In file s_m_amat.c
void scalar_mult_adj_matrix(matrix *src, Real s, matrix *c);

// In file s_m_a_mat.c
void scalar_mult_sum_matrix(matrix *b, Real s, matrix *c);
void scalar_mult_add_matrix(matrix *a, matrix *b, Real s, matrix *c);

// In file s_m_a_amat.c
void scalar_mult_sum_adj_matrix(matrix *b, Real s, matrix *c);

// In file s_m_s_mat.c
void scalar_mult_dif_matrix(matrix *b, Real s, matrix *c);

// In file s_m_s_amat.c
void scalar_mult_dif_adj_matrix(matrix *b, Real s, matrix *c);

// In file cs_a_d_mat.c
void c_scalar_add_diag(matrix *a, complex *f);

// In file cs_m_mat.c
void c_scalar_mult_mat(matrix *b, complex *s, matrix *c);

// In file cs_m_a_mat.c
void c_scalar_mult_sum_mat(matrix *b, complex *s, matrix *c);

// In file cs_m_a_amat.c
void c_scalar_mult_sum_adj_mat(matrix *b, complex *s, matrix *c);
void c_scalar_mult_dif_adj_mat(matrix *b, complex *s, matrix *c);

// In file cs_m_a_mata.c
void c_scalar_mult_sum_mat_adj(matrix *b, complex *s, matrix *c);

// In file cs_m_s_mat.c
void c_scalar_mult_dif_mat(matrix *b, complex *s, matrix *c);

// In file m_mat_nn.c
void mult_nn_sum(matrix *a, matrix *b, matrix *c);
void mult_nn_dif(matrix *a, matrix *b, matrix *c);
void mult_nn(matrix *a, matrix *b, matrix *c);

// In file m_mat_na.c
void mult_na_sum(matrix *a, matrix *b, matrix *c);
void mult_na_dif(matrix *a, matrix *b, matrix *c);
void mult_na(matrix *a, matrix *b, matrix *c);

// In file m_mat_an.c
void mult_an_sum(matrix *a, matrix *b, matrix *c);
void mult_an_dif(matrix *a, matrix *b, matrix *c);
void mult_an(matrix *a, matrix *b, matrix *c);

// In file s_m_mat_nn.c
void scalar_mult_nn_sum(matrix *a, matrix *b, Real s, matrix *c);
void scalar_mult_nn_dif(matrix *a, matrix *b, Real s, matrix *c);
void scalar_mult_nn(matrix *a, matrix *b, Real s, matrix *c);

// In file s_m_mat_na.c
void scalar_mult_na_sum(matrix *a, matrix *b, Real s, matrix *c);
void scalar_mult_na_dif(matrix *a, matrix *b, Real s, matrix *c);
void scalar_mult_na(matrix *a, matrix *b, Real s, matrix *c);

// In file s_m_mat_an.c
void scalar_mult_an_sum(matrix *a, matrix *b, Real s, matrix *c);
void scalar_mult_an_dif(matrix *a, matrix *b, Real s, matrix *c);
void scalar_mult_an(matrix *a, matrix *b, Real s, matrix *c);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Anti-hermitian matrix routines
// In file make_ahmat.c
void make_anti_hermitian(matrix *m, anti_hermitmat *ah);

// In file rand_ahmat.c
void random_anti_hermitian(anti_hermitmat *ah, double_prn *prn_pt);

// In file uncmp_ahmat.c
void uncompress_anti_hermitian(anti_hermitmat *ah, matrix *m);

// In file cmp_ahmat.c
void compress_anti_hermitian(matrix *m, anti_hermitmat *ah);

// In file dump_ahmat.c
void dump_ahmat(anti_hermitmat *ahm);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Routines mixing SU(2) and U(N)
// In file l_su2_hit_n.c
void left_su2_hit_n(su2_matrix *u, int p, int q, matrix *link);

// In file r_su2_hit_a.c
void right_su2_hit_a(su2_matrix *u, int p, int q, matrix *link);

// In file dumpsu2.c
void dumpsu2(su2_matrix *u);

// In file m_su2_mat_vec_n.c
void mult_su2_mat_vec_elem_n(su2_matrix *u, complex *x0, complex *x1);

// In file m_su2_mat_vec_a.c
void mult_su2_mat_vec_elem_a(su2_matrix *u, complex *x0, complex *x1);
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Routines for funmatrices
// in file fun_o_functions.c
void fun_o_sub(funmatrix *a, funmatrix *b, funmatrix *c);
void fun_o_sum(funmatrix *a, funmatrix *b);
void fun_o_scalar_mult_sum(funmatrix *a, Real R, funmatrix *b);
void funmat_copy(funmatrix *a, funmatrix *b);

// in file fun_adjoint.c
void fun_scaled_adjoint(funmatrix *a, Real R, funmatrix *b);
void fun_scaled_adjoint_sum(funmatrix *a, Real R, funmatrix *b);

// in file fun_clear_mat.c
void fun_clear_mat(funmatrix *m);

// in file fun_nn_functions.c
Real fun_nn_realtrace(funmatrix *a, funmatrix *b);

// in file fun_ta_functions.c
void fun_ta_mult_nn(funmatrix *a, matrix *b, funmatrix *c);

// in file fun_tt_functions.c
void fun_tt_mult_an(funmatrix *a, funmatrix *b, matrix *c);
void fun_tt_mult_an_sum(funmatrix *a, funmatrix *b, matrix *c);


// -----------------------------------------------------------------
// Miscellaneous routines
// In file gaussrand.c
Real gaussian_rand_no(double_prn *prn_pt);

// In file z2rand.c
Real Z2_rand_no(double_prn *prn_pt);

// In file byterevn.c
#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);

#endif
// -----------------------------------------------------------------
