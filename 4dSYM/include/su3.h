// -----------------------------------------------------------------
// Defines and subroutine declarations for U(N) and SU(N)
// with fermions in the DIMF-dimensional adjoint rep
// The original names now refer to objects of dimension DIMF or DIMFxDIMF
// New objects with suffix _f have dimension NCOL or NCOLxNCOL
#ifndef _SUN_H
#define _SUN_H

#include "../include/complex.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermions in adjoint rep
#define FREP adjoint
// #define REALREP      // Someday use for adjoint

// U(2)
#define NCOL 2
#define DIMF 4
#define NUMGEN 4

// U(3)
//#define NCOL 3
//#define DIMF 9
//#define NUMGEN 9

// U(4)
//#define NCOL 4
//#define DIMF 16
//#define NUMGEN 16

// SU(2)
//#define NCOL 2
//#define DIMF 3
//#define NUMGEN 3

// Restrictions due to anti_hermitian structs
#if (NCOL > 4)
  #error "NCOL must be <5!"
#endif

#if (NCOL!=3 || DIMF!=3)
#ifdef FAST
  #error "FAST only works if NCOL=DIMF=3!"
#endif
#endif

typedef struct { fcomplex e[NCOL][NCOL]; } fsu3_matrix_f;
#ifndef REALREP
typedef struct { fcomplex e[DIMF][DIMF]; } fsu3_matrix;
#else
typedef struct { float    e[DIMF][DIMF]; } fsu3_matrix;
#endif
typedef struct { fcomplex c[NCOL]; } fsu3_vector_f;
typedef struct { fcomplex c[DIMF]; } fsu3_vector;

// Anti-hermitian matrices require NCOL<=4
typedef struct {
#if (NCOL == 2)
  fcomplex m01;
  float m00im, m11im;
#endif
#if (NCOL == 3)
  fcomplex m01, m02, m12;
  float m00im, m11im, m22im;
  float space;
#endif
#if (NCOL == 4)
  fcomplex m01, m02, m03, m12, m13, m23;
  float m00im, m11im, m22im, m33im;
#endif
} fanti_hermitmat;

typedef struct { dcomplex e[NCOL][NCOL]; } dsu3_matrix_f;
#ifndef REALREP
typedef struct { dcomplex e[DIMF][DIMF]; } dsu3_matrix;
#else
typedef struct { double e[DIMF][DIMF]; } dsu3_matrix;
#endif
typedef struct { dcomplex c[NCOL]; } dsu3_vector_f;
typedef struct { dcomplex c[DIMF]; } dsu3_vector;
typedef struct {
#if (NCOL == 2)
  dcomplex m01;
  double m00im, m11im;
#endif
#if (NCOL == 3)
  dcomplex m01, m02, m12;
  double m00im, m11im, m22im;
  double space;
#endif
#if (NCOL == 4)
  dcomplex m01, m02, m03, m12, m13, m23;
  double m00im, m11im, m22im, m33im;
#endif
} danti_hermitmat;

#if (PRECISION==1)
#define su3_matrix_f  fsu3_matrix_f
#define su3_matrix      fsu3_matrix
#define su3_vector_f  fsu3_vector_f
#define su3_vector      fsu3_vector
#define anti_hermitmat  fanti_hermitmat
#else
#define su3_matrix_f  dsu3_matrix_f
#define su3_matrix      dsu3_matrix
#define su3_vector_f  dsu3_vector_f
#define su3_vector      dsu3_vector
#define anti_hermitmat  danti_hermitmat
#endif

// Used if REALREP
typedef struct { complex e[DIMF][DIMF]; } complex_su3_matrix;

// SU(2)
typedef struct { complex e[2][2]; } su2_matrix;

#define GAMMAFIVE -1    // Some integer which is not a direction
#define PLUS 1          // Flags for selecting M or M_adjoint
#define MINUS -1

// Macros to multiply complex numbers by -1 and +-i
#define TIMESMINUSONE(a,b) { (b).real =  -(a).real; (b).imag = -(a).imag; }
#define TIMESPLUSI(a,b) { (b).real = -(a).imag; (b).imag =  (a).real; }
#define TIMESMINUSI(a,b) { (b).real =  (a).imag; (b).imag = -(a).real; }
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Subroutine definitions
// Fundamental rep vector operations
// In file clearvec_f.c
void clearvec_f(su3_vector_f *v);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fundamental rep matrix operations (including anti-hermitian matrices)
// In file clear_mat_f
void clear_su3mat_f(su3_matrix_f *dest);

// In file trace_su3_f.c
complex trace_su3_f(su3_matrix_f *a);

// In file realtr_f.c
Real realtrace_su3_f(su3_matrix_f *a, su3_matrix_f *b);

// In file complextr_f.c
complex complextrace_su3_f(su3_matrix_f *a, su3_matrix_f *b);

// b <-- a, in file su3mat_copy_f.c
void su3mat_copy_f(su3_matrix_f *a, su3_matrix_f *b);

// b <-- adag, in file su3_adjoint_f.c
void su3_adjoint_f(su3_matrix_f *a, su3_matrix_f *b);

// In file addmat_f.c
void add_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);

// In file submat_f.c
void sub_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);

// In file s_m_mat_f.c
void scalar_mult_su3_matrix_f(su3_matrix_f *src, Real scalar,
                              su3_matrix_f *dest);

// In file s_m_s_mat_f.c
void scalar_mult_sub_su3_matrix_f(su3_matrix_f *src1, su3_matrix_f *src2,
                                  Real scalar, su3_matrix_f *dest);

// In file s_a_d_mat_f.c
void scalar_add_diag_su3_f(su3_matrix_f *a, Real s);

// In file cs_a_d_mat_f.c
void c_scalar_add_diag_su3_f(su3_matrix_f *a, complex *f);

// In file cs_m_mat_f.c
void c_scalar_mult_su3mat_f(su3_matrix_f *b, complex *s, su3_matrix_f *c);

// In file cs_m_a_mat_f.c
void c_scalar_mult_add_su3mat_f(su3_matrix_f *m1, su3_matrix_f *m2,
                                complex *phase, su3_matrix_f *m3);

// In file dumpmat_f.c
void dumpmat_f(su3_matrix_f *m);

// In file m_mat_nn_f.c
void mult_su3_nn_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);

// In file m_mat_na_f.c
void mult_su3_na_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);

// In file m_mat_an_f.c
void mult_su3_an_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);

// In file s_m_a_mat_f.c
void scalar_mult_add_su3_matrix_f(su3_matrix_f *src1, su3_matrix_f *src2,
                                  Real scalar, su3_matrix_f *dest);

// Relate fundamental matrices and anti-hermitian matrices
// In file make_ahmat.c
void make_anti_hermitian(su3_matrix_f *m, anti_hermitmat *ah);

// In file rand_ahmat.c
void random_anti_hermitian(anti_hermitmat *ah, double_prn *prn_pt);

// In file uncmp_ahmat.c
void uncompress_anti_hermitian(anti_hermitmat *ah, su3_matrix_f *m);

// In file cmp_ahmat.c
void compress_anti_hermitian(su3_matrix_f *m, anti_hermitmat *ah);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion rep vector operations
// In file dumpvec.c
void dumpvec(su3_vector *v);

// In file clearvec.c
void clearvec(su3_vector *v);

// In file msq_su3vec.c
Real magsq_su3vec(su3_vector *v);

// In file su3vec_copy.c
void su3vec_copy(su3_vector *a, su3_vector *b);

// In file su3_rdot.c
Real su3_rdot(su3_vector *a, su3_vector *b);

// In file su3_dot.c
complex su3_dot(su3_vector *a, su3_vector *b);

// In file addvec.c
void add_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c);

// In file subvec.c
void sub_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c);

// In file s_m_vec.c
void scalar_mult_su3_vector(su3_vector *src, Real scalar,
                            su3_vector *dest);

// In file s_m_sum_vec.c
void scalar_mult_sum_su3_vector(su3_vector *src1, su3_vector *src2,
                                Real scalar);

// In file s_m_a_vec.c
void scalar_mult_add_su3_vector(su3_vector *src1, su3_vector *src2,
                                Real scalar, su3_vector *dest);

// In file s_m_s_vec.c
void scalar_mult_sub_su3_vector(su3_vector *src1, su3_vector *src2,
                                Real scalar, su3_vector *dest);

// In file cs_m_vec.c
void c_scalar_mult_su3vec(su3_vector *src, complex *phase, su3_vector *dest);

// In file cs_m_a_vec.c
void c_scalar_mult_add_su3vec(su3_vector *v1, complex *phase, su3_vector *v2);

// In file cs_m_s_vec.c
void c_scalar_mult_sub_su3vec(su3_vector *v1, complex *phase, su3_vector *v2);

// In file sub4vecs.c
void sub_four_su3_vecs(su3_vector *a, su3_vector *b1, su3_vector *b2,
                       su3_vector *b3, su3_vector *b4);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion rep matrix operations
// In file trace_su3.c
complex trace_su3(su3_matrix *a);

// In file realtr.c
Real realtrace_su3(su3_matrix *a, su3_matrix *b);

// In file complextr.c
complex complextrace_su3(su3_matrix *a, su3_matrix *b);

// In file addmat.c
void add_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file submat.c
void sub_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file s_m_mat.c
void scalar_mult_su3_matrix(su3_matrix *src, Real scalar, su3_matrix *dest);

// In file s_m_a_mat.c
void scalar_mult_add_su3_matrix(su3_matrix *src1, su3_matrix *src2,
                                Real scalar, su3_matrix *dest);

// In file s_m_s_mat.c
void scalar_mult_sub_su3_matrix(su3_matrix *src1, su3_matrix *src2,
                                Real scalar, su3_matrix *dest);

// In file cs_m_mat.c
void c_scalar_mult_su3mat(su3_matrix *src, complex *scalar,
                          su3_matrix *dest);

// In file cs_m_a_mat.c
void c_scalar_mult_add_su3mat(su3_matrix *src1, su3_matrix *src2,
                              complex *scalar, su3_matrix *dest);

// In file cs_m_s_mat.c
void c_scalar_mult_sub_su3mat(su3_matrix *src1, su3_matrix *src2,
                              complex *scalar, su3_matrix *dest);

// In file su3_adjoint.c
void su3_adjoint(su3_matrix *a, su3_matrix *b);

// In file clear_mat.c
void clear_su3mat(su3_matrix *dest);

// In file su3mat_copy.c
void su3mat_copy(su3_matrix *a, su3_matrix *b);

// In file dumpmat.c
void dumpmat(su3_matrix *m);

// In file m_mat_nn.c
void mult_su3_nn(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file m_mat_na.c
void mult_su3_na(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file m_mat_an.c
void mult_su3_an(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// c <-- a * b, in file m_matvec.c
void mult_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- c + a * b, in file m_matvec_s.c
void mult_su3_mat_vec_sum(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- c - a * b, in file m_matvec_ns.c
void mult_su3_mat_vec_nsum(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- adag * b, in file m_amatvec.c
void mult_adj_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- c + adag * b, in file m_amatvec_s.c
void mult_adj_su3_mat_vec_sum(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- c - adag * b, in file m_amatvec_ns.c
void mult_adj_su3_mat_vec_nsum(su3_matrix *a, su3_vector *b, su3_vector *c);

// c <-- b * a, in file m_vecmat.c
void mult_su3_vec_mat(su3_vector *b, su3_matrix *a, su3_vector *c);

// c <-- b * adag, in file m_vecamat.c
void mult_su3_vec_adj_mat(su3_vector *b, su3_matrix *a, su3_vector *c);

// c <-- c + b * adag, in file m_vecamat_s.c
void mult_su3_vec_adj_mat_sum(su3_vector *b, su3_matrix *a, su3_vector *c);

// In file m_amv_4dir.c
void mult_adj_su3_mat_vec_4dir(su3_matrix *a, su3_vector *b, su3_vector *c);

// In file m_amv_4vec.c
void mult_adj_su3_mat_4vec(su3_matrix *mat, su3_vector *src,
                           su3_vector *dest0, su3_vector *dest1,
                           su3_vector *dest2, su3_vector *dest3);

// In file m_mv_s_4dir.c
void mult_su3_mat_vec_sum_4dir(su3_matrix *a, su3_vector *b0,
                               su3_vector *b1, su3_vector *b2,
                               su3_vector *b3, su3_vector *c);

// In file su3_proj.c
void su3_projector(su3_vector *a, su3_vector *b, su3_matrix *c);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Routines mixing SU(2) and SU(N)
// In file l_su2_hit_n.c
void left_su2_hit_n(su2_matrix *u, int p, int q, su3_matrix *link);

// In file r_su2_hit_a.c
void right_su2_hit_a(su2_matrix *u, int p, int q, su3_matrix *link);

// In file l_su2_hit_n_f.c
void left_su2_hit_n_f(su2_matrix *u, int p, int q, su3_matrix_f *link);

// In file r_su2_hit_a_f.c
void right_su2_hit_a_f(su2_matrix *u, int p, int q, su3_matrix_f *link);

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
