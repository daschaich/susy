// -----------------------------------------------------------------
// Defines and subroutine declarations for N=4 SYM with gauge group U(NCOL)
// and fermions in the DIMF-dimensional adjoint rep
// The original names now refer to objects of dimension DIMF or DIMFxDIMF
// New objects with suffix _f have dimension NCOL or NCOLxNCOL
#ifndef _SUN_H
#define _SUN_H

#include "../include/complex.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// N=4 SYM fermions always in adjoint rep
// Gauge group U(NCOL) and size of adjoint rep DIMF = NCOL^2
#define NCOL 2
#define DIMF 4

#if (NCOL!=3 || DIMF!=3)
#ifdef FAST
  #error "FAST only works if NCOL=DIMF=3!"
#endif
#endif

typedef struct { fcomplex e[NCOL][NCOL]; } fmatrix_f;
typedef struct { fcomplex e[DIMF][DIMF]; } fmatrix;
typedef struct { fcomplex c[NCOL]; } fvector_f;
typedef struct { fcomplex c[DIMF]; } fvector;

// Anti-hermitian matrices for general NCOL
typedef struct {
  fcomplex m[NCOL * (NCOL - 1) / 2];
  float im_diag[NCOL];
} fanti_hermitmat;

typedef struct { dcomplex e[NCOL][NCOL]; } dmatrix_f;
typedef struct { dcomplex e[DIMF][DIMF]; } dmatrix;
typedef struct { dcomplex c[NCOL]; } dvector_f;
typedef struct { dcomplex c[DIMF]; } dvector;
typedef struct {
  dcomplex m[NCOL * (NCOL - 1) / 2];
  double im_diag[NCOL];
} danti_hermitmat;

#if (PRECISION == 1)
#define matrix_f    fmatrix_f
#define matrix      fmatrix
#define vector_f    fvector_f
#define vector      fvector
#define anti_hermitmat  fanti_hermitmat
#else
#define matrix_f    dmatrix_f
#define matrix      dmatrix
#define vector_f    dvector_f
#define vector      dvector
#define anti_hermitmat  danti_hermitmat
#endif

// Need SU(2) matrices for any U(N), e.g. for gauge hits when gauge-fixing
typedef struct { complex e[2][2]; } su2_matrix;

#define GAMMAFIVE -1    // Some integer which is not a direction
#define PLUS 1          // Flags for selecting M or M_adjoint
#define MINUS -1
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Subroutine definitions
// Fundamental rep vector operations
// In file clearvec_f.c
void clearvec_f(vector_f *v);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fundamental rep matrix operations (including anti-hermitian matrices)
// In file clear_mat_f
void clear_mat_f(matrix_f *c);

// In file trace_f.c
complex trace_f(matrix_f *a);

// In file realtr_f.c
Real realtrace_f(matrix_f *a, matrix_f *b);

// In file complextr_f.c
complex complextrace_f(matrix_f *a, matrix_f *b);

// b <-- a, in file mat_copy_f.c
void mat_copy_f(matrix_f *a, matrix_f *b);

// b <-- adag, in file adjoint_f.c
void adjoint_f(matrix_f *a, matrix_f *b);

// In file addmat_f.c
void sum_matrix_f(matrix_f *b, matrix_f *c);
void add_matrix_f(matrix_f *a, matrix_f *b, matrix_f *c);

// In file addamat_f.c
void add_adj_matrix_f(matrix_f *a, matrix_f *b, matrix_f *c);

// In file submat_f.c
void dif_matrix_f(matrix_f *b, matrix_f *c);
void sub_matrix_f(matrix_f *a, matrix_f *b, matrix_f *c);

// In file subamat_f.c
void sub_adj_matrix_f(matrix_f *a, matrix_f *b, matrix_f *c);

// In file s_a_d_mat_f.c
void scalar_add_diag_f(matrix_f *a, Real s);

// In file s_m_mat_f.c
void scalar_mult_matrix_f(matrix_f *src, Real s, matrix_f *c);

// In file s_m_amat_f.c
void scalar_mult_adj_matrix_f(matrix_f *src, Real s, matrix_f *c);

// In file s_m_a_mat_f.c
void scalar_mult_sum_matrix_f(matrix_f *b, Real s, matrix_f *c);
void scalar_mult_add_matrix_f(matrix_f *a, matrix_f *b, Real s, matrix_f *c);

// In file s_m_a_amat_f.c
void scalar_mult_sum_adj_matrix_f(matrix_f *b, Real s, matrix_f *c);

// In file s_m_s_mat_f.c
void scalar_mult_dif_matrix_f(matrix_f *b, Real s, matrix_f *c);

// In file s_m_s_amat_f.c
void scalar_mult_dif_adj_matrix_f(matrix_f *b, Real s, matrix_f *c);

// In file cs_a_d_mat_f.c
void c_scalar_add_diag_f(matrix_f *a, complex *f);

// In file cs_m_mat_f.c
void c_scalar_mult_mat_f(matrix_f *b, complex *s, matrix_f *c);

// In file cs_m_a_mat_f.c
void c_scalar_mult_sum_mat_f(matrix_f *b, complex *s, matrix_f *c);

// In file cs_m_a_amat_f.c
void c_scalar_mult_sum_adj_mat_f(matrix_f *b, complex *s, matrix_f *c);

// In file cs_m_a_mata_f.c
void c_scalar_mult_sum_mat_adj_f(matrix_f *b, complex *s, matrix_f *c);

// In file cs_m_s_mat_f.c
void c_scalar_mult_dif_mat_f(matrix_f *b, complex *s, matrix_f *c);

// In file dumpmat_f.c
void dumpmat_f(matrix_f *m);

// In file m_mat_nn_f.c
void mult_nn_f(matrix_f *a, matrix_f *b, matrix_f *c);

// In file m_mat_na_f.c
void mult_na_sum_f(matrix_f *a, matrix_f *b, matrix_f *c);
void mult_na_f(matrix_f *a, matrix_f *b, matrix_f *c);

// In file m_mat_an_f.c
void mult_an_f(matrix_f *a, matrix_f *b, matrix_f *c);

// Relate fundamental matrices and anti-hermitian matrices
// In file make_ahmat.c
void make_anti_hermitian(matrix_f *m, anti_hermitmat *ah);

// In file rand_ahmat.c
void random_anti_hermitian(anti_hermitmat *ah, double_prn *prn_pt);

// In file uncmp_ahmat.c
void uncompress_anti_hermitian(anti_hermitmat *ah, matrix_f *m);

// In file cmp_ahmat.c
void compress_anti_hermitian(matrix_f *m, anti_hermitmat *ah);

// In file dump_ahmat.c
void dump_ahmat(anti_hermitmat *ahm);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion rep vector operations
// In file dumpvec.c
void dumpvec(vector *v);

// In file clearvec.c
void clearvec(vector *v);

// In file msq_vec.c
Real magsq_vec(vector *v);

// In file vec_copy.c
void vec_copy(vector *a, vector *b);

// In file rdot.c
Real rdot(vector *a, vector *b);

// In file dot.c
complex dot(vector *a, vector *b);

// In file addvec.c
void add_vector(vector *a, vector *b, vector *c);

// In file subvec.c
void sub_vector(vector *a, vector *b, vector *c);

// In file s_m_vec.c
void scalar_mult_vector(vector *src, Real s, vector *c);

// In file s_m_sum_vec.c
void scalar_mult_sum_vector(vector *a, vector *b, Real scalar);

// In file s_m_a_vec.c
void scalar_mult_add_vector(vector *a, vector *b,
                            Real s, vector *c);

// In file s_m_s_vec.c
void scalar_mult_sub_vector(vector *a, vector *b,
                            Real s, vector *c);

// In file cs_m_vec.c
void c_scalar_mult_vec(vector *src, complex *phase, vector *c);

// In file cs_m_a_vec.c
void c_scalar_mult_add_vec(vector *v1, complex *phase, vector *v2);

// In file cs_m_s_vec.c
void c_scalar_mult_sub_vec(vector *v1, complex *phase, vector *v2);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion rep matrix operations
// In file trace.c
complex trace(matrix *a);

// In file realtr.c
Real realtrace(matrix *a, matrix *b);

// In file complextr.c
complex complextrace(matrix *a, matrix *b);

// In file addmat.c
void add_matrix(matrix *a, matrix *b, matrix *c);

// In file submat.c
void sub_matrix(matrix *a, matrix *b, matrix *c);

// In file s_m_mat.c
void scalar_mult_matrix(matrix *src, Real s, matrix *c);

// In file s_m_a_mat.c
void scalar_mult_add_matrix(matrix *a, matrix *b, Real s, matrix *c);

// In file s_m_s_mat.c
void scalar_mult_sub_matrix(matrix *a, matrix *b, Real s, matrix *c);

// In file cs_m_mat.c
void c_scalar_mult_mat(matrix *src, complex *s, matrix *c);

// In file cs_m_a_mat.c
void c_scalar_mult_add_mat(matrix *a, matrix *b, complex *s, matrix *c);

// In file cs_m_s_mat.c
void c_scalar_mult_sub_mat(matrix *a, matrix *b, complex *s, matrix *c);

// In file adjoint.c
void adjoint(matrix *a, matrix *b);

// In file clear_mat.c
void clear_mat(matrix *c);

// In file mat_copy.c
void mat_copy(matrix *a, matrix *b);

// In file dumpmat.c
void dumpmat(matrix *m);

// In file m_mat_nn.c
void mult_nn(matrix *a, matrix *b, matrix *c);

// In file m_mat_na.c
void mult_na(matrix *a, matrix *b, matrix *c);

// In file m_mat_an.c
void mult_an(matrix *a, matrix *b, matrix *c);
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Routines involving both fermion-rep matrices and vectors
// c <-- a * b, in file m_matvec.c
void mult_mat_vec(matrix *a, vector *b, vector *c);

// c <-- c + a * b, in file m_matvec_s.c
void mult_mat_vec_sum(matrix *a, vector *b, vector *c);

// c <-- c - a * b, in file m_matvec_ns.c
void mult_mat_vec_dif(matrix *a, vector *b, vector *c);

// c <-- adag * b, in file m_amatvec.c
void mult_adj_mat_vec(matrix *a, vector *b, vector *c);

// c <-- c + adag * b, in file m_amatvec_s.c
void mult_adj_mat_vec_sum(matrix *a, vector *b, vector *c);

// c <-- c - adag * b, in file m_amatvec_ns.c
void mult_adj_mat_vec_dif(matrix *a, vector *b, vector *c);

// c <-- b * a, in file m_vecmat.c
void mult_vec_mat(vector *b, matrix *a, vector *c);

// c <-- c + b * a, in file m_vecmat_s.c
void mult_vec_mat_sum(vector *b, matrix *a, vector *c);

// c <-- c - b * a, in file m_vecmat_ns.c
void mult_vec_mat_dif(vector *b, matrix *a, vector *c);

// c <-- b * adag, in file m_vecamat.c
void mult_vec_adj_mat(vector *b, matrix *a, vector *c);

// c <-- c + b * adag, in file m_vecamat_s.c
void mult_vec_adj_mat_sum(vector *b, matrix *a, vector *c);

// c <-- c - b * adag, in file m_vecamat_ns.c
void mult_vec_adj_mat_dif(vector *b, matrix *a, vector *c);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Routines mixing SU(2) and U(N)
// In file l_su2_hit_n.c
void left_su2_hit_n(su2_matrix *u, int p, int q, matrix *link);

// In file r_su2_hit_a.c
void right_su2_hit_a(su2_matrix *u, int p, int q, matrix *link);

// In file l_su2_hit_n_f.c
void left_su2_hit_n_f(su2_matrix *u, int p, int q, matrix_f *link);

// In file r_su2_hit_a_f.c
void right_su2_hit_a_f(su2_matrix *u, int p, int q, matrix_f *link);

// In file dumpsu2.c
void dumpsu2(su2_matrix *u);

// In file m_su2_mat_vec_n.c
void mult_su2_mat_vec_elem_n(su2_matrix *u, complex *x0, complex *x1);

// In file m_su2_mat_vec_a.c
void mult_su2_mat_vec_elem_a(su2_matrix *u, complex *x0, complex *x1);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Miscellaneous routines
// In file gaussrand.c
Real gaussian_rand_no(double_prn *prn_pt);

#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);

#endif
// -----------------------------------------------------------------
