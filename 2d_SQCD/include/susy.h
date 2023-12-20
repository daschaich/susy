// -----------------------------------------------------------------
// Defines and subroutine declarations for N=4 SYM and its dim'l reductions
#ifndef _SUSY_H
#define _SUSY_H

#include "../include/complex.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// SQCD has fermions in two reps
// Adjoint rep, as are dim-reduced scalars
// This file is only included for d-dimensional systems with d>=2
// The gauge group is U(NCOL) and size of adjoint rep DIMF = NCOL^2
#define NCOL 2
#define DIMF (NCOL * NCOL)
#define N_OFFDIAG (NCOL * (NCOL - 1) / 2)

// Fundamental rep
// for consistency all field will be represented as NCOL * NCOLF
// however the vector field should be understood as the transpose
// NCOLF is the number of flavours
// funmatrix are the NCOL x NCOLF matrices
// funamatrix are the NCOLF x NCOL matrices
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
typedef struct { fcomplex e[NCOLF][NCOL]; } ffunamatrix;

typedef struct { dcomplex e[NCOL][NCOL]; } dmatrix;
typedef struct { dcomplex c[NCOL]; } dvector;
typedef struct {
  dcomplex m[N_OFFDIAG];
  double im_diag[NCOL];
} danti_hermitmat;

typedef struct { dcomplex e[NCOL][NCOLF]; } dfunmatrix;
typedef struct { dcomplex e[NCOLF][NCOL]; } dfunamatrix;

#if (PRECISION == 1)
#define matrix    fmatrix
#define vector    fvector
#define anti_hermitmat  fanti_hermitmat
#define funmatrix       ffunmatrix
#define funamatrix       ffunamatrix
#else
#define matrix    dmatrix
#define vector    dvector
#define anti_hermitmat  danti_hermitmat
#define funmatrix       dfunmatrix
#define funamatrix       dfunamatrix
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
void scalar_sub_diag(matrix *a, Real s);

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
// For now all in funmatrix_functions.c but plan to split the file
void fun_sum_matrix(funmatrix *b, funmatrix *c);
void fun_add_matrix(funmatrix *a, funmatrix *b, funmatrix *c);

void funa_sum_matrix(funamatrix *b, funamatrix *c);
void funa_add_matrix(funamatrix *a, funamatrix *b, funamatrix *c);

void fun_adjoint(funmatrix *a, funamatrix *b);
void fun_neg_adjoint(funmatrix *a, funamatrix *b);

void funa_adjoint(funamatrix *a, funmatrix *b);
void funa_neg_adjoint(funamatrix *a, funmatrix *b);

void fun_clear_mat(funmatrix *m);
void funa_clear_mat(funamatrix *m);

void fun_c_scalar_mult_sum_adj_mat(funmatrix *b, complex *s, funamatrix *c);
void fun_c_scalar_mult_dif_adj_mat(funmatrix *b, complex *s, funamatrix *c);

void funa_c_scalar_mult_sum_adj_mat(funamatrix *b, complex *s, funmatrix *c);
void funa_c_scalar_mult_dif_adj_mat(funamatrix *b, complex *s, funmatrix *c);

void fun_c_scalar_mult_sum_mat(funmatrix *b, complex *s, funmatrix *c);
void funa_c_scalar_mult_sum_mat(funamatrix *b, complex *s, funamatrix *c);

void fun_c_scalar_mult_sum_mat_adj(funmatrix *b, complex *s, funamatrix *c);
void funa_c_scalar_mult_sum_mat_adj(funamatrix *b, complex *s, funmatrix *c);

void fun_c_scalar_mult_mat(funmatrix *b, complex *s, funmatrix *c);
void funa_c_scalar_mult_mat(funamatrix *b, complex *s, funamatrix *c);

void fun_c_scalar_mult_dif_mat(funmatrix *b, complex *s, funmatrix *c);
void funa_c_scalar_mult_dif_mat(funamatrix *b, complex *s, funamatrix *c);

void fun_dumpmat(funmatrix *m);
void funa_dumpmat(funamatrix *m);

void funa_mult_an_sum(funamatrix *a, funamatrix *b, matrix *c);
void funa_mult_an_dif(funamatrix *a, funamatrix *b, matrix *c);
void funa_mult_an(funamatrix *a, funamatrix *b, matrix *c);

void fun_mult_na_sum(funmatrix *a, funmatrix *b, matrix *c);
void fun_mult_na_dif(funmatrix *a, funmatrix *b, matrix *c);
void fun_mult_na(funmatrix *a, funmatrix *b, matrix *c);

void fun_mult_nn_sum(funmatrix *a, funamatrix *b, matrix *c);
void fun_mult_nn_dif(funmatrix *a, funamatrix *b, matrix *c);
void fun_mult_nn(funmatrix *a, funamatrix *b, matrix *c);

void fun_mat_copy(funmatrix *src, funmatrix *dest);
void funa_mat_copy(funamatrix *src, funamatrix *dest);

void fun_scalar_mult_sum_adj_matrix(funmatrix *b, Real s, funamatrix *c);
void funa_scalar_mult_sum_adj_matrix(funamatrix *b, Real s, funmatrix *c);

void fun_scalar_mult_sum_matrix(funmatrix *b, Real s, funmatrix *c);
void fun_scalar_mult_add_matrix(funmatrix *a, funmatrix *b, Real s, funmatrix *c);

void funa_scalar_mult_sum_matrix(funamatrix *b, Real s, funamatrix *c);
void funa_scalar_mult_add_matrix(funamatrix *a, funamatrix *b, Real s, funamatrix *c);

void fun_scalar_mult_adj_matrix(funmatrix *a, Real s, funamatrix *b);
void funa_scalar_mult_adj_matrix(funamatrix *a, Real s, funmatrix *b);

void fun_scalar_mult_matrix(funmatrix *a, Real s, funmatrix *b);
void funa_scalar_mult_matrix(funamatrix *a, Real s, funamatrix *b);

void funa_scalar_mult_an_sum(funamatrix *a, funamatrix *b, Real s, matrix *c);
void funa_scalar_mult_an_dif(funamatrix *a, funamatrix *b, Real s, matrix *c);
void funa_scalar_mult_an(funamatrix *a, funamatrix *b, Real s, matrix *c);

void fun_scalar_mult_na_sum(funmatrix *a, funmatrix *b, Real s, matrix *c);
void fun_scalar_mult_na_dif(funmatrix *a, funmatrix *b, Real s, matrix *c);
void fun_scalar_mult_na(funmatrix *a, funmatrix *b, Real s, matrix *c);

void fun_scalar_mult_nn_sum(funmatrix *a, funamatrix *b, Real s, matrix *c);
void fun_scalar_mult_nn_dif(funmatrix *a, funamatrix *b, Real s, matrix *c);
void fun_scalar_mult_nn(funmatrix *a, funamatrix *b, Real s, matrix *c);

void fun_scalar_mult_dif_adj_matrix(funmatrix *b, Real s, funamatrix *c);
void funa_scalar_mult_dif_adj_matrix(funamatrix *b, Real s, funmatrix *c);

void fun_scalar_mult_dif_matrix(funmatrix *b, Real s, funmatrix *c);
void funa_scalar_mult_dif_matrix(funamatrix *b, Real s, funamatrix *c);

void fun_dif_matrix(funmatrix *b, funmatrix *c);
void fun_sub_matrix(funmatrix *a, funmatrix *b, funmatrix *c);

void funa_dif_matrix(funamatrix *b, funamatrix *c);
void funa_sub_matrix(funamatrix *a, funamatrix *b, funamatrix *c);

void fun_mat_prod(funmatrix *a, matrix *b, funmatrix *c);
void fun_mat_prod_dif(funmatrix *a, matrix *b, funmatrix *c);

void funa_mat_prod(funamatrix *a, matrix *b, funamatrix *c);
void funa_mat_prod_sum(funamatrix *a, matrix *b, funamatrix *c);
void funa_mat_prod_dif(funamatrix *a, matrix *b, funamatrix *c);
void funa_mat_mult_an(funmatrix *a, matrix *b, funamatrix *c);
void funa_scalar_mat_mult_an(funmatrix *a, matrix *b, Real s, funamatrix *c);
void funa_mat_prod_an_sum(funmatrix *a, matrix *b, funamatrix *c);
void funa_mat_prod_an_dif(funmatrix *a, matrix *b, funamatrix *c);

void funa_scalar_mat_prod_an_dif(funmatrix *a,matrix *b,Real s,funamatrix *c);

Real fun_realtrace(funmatrix *a, funmatrix *b);
Real funa_realtrace(funamatrix *a, funamatrix *b);

complex fun_complextrace_an(funmatrix *a, funmatrix *b);
complex fun_complextrace_na(funmatrix *a, funmatrix *b);

complex funa_complextrace_an(funamatrix *a, funamatrix *b);
complex funa_complextrace_na(funamatrix *a, funamatrix *b);

void funa_dif_adj_matrix(funmatrix *b, funamatrix *c);
// -----------------------------------------------------------------

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
