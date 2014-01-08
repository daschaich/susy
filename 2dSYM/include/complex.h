// -----------------------------------------------------------------
// Complex Numbers
// Typedefs are included for single- and double-precision complex numbers
// The functions cannot be overloaded, so there are separate routines
// All of the macros, however, will work with both types and mix types freely
#ifndef _SUSY_COMPLEX_H
#define _SUSY_COMPLEX_H

// The following functions are provided in single and double precision
//   complex cmplx(Real r, Real i);   (r, i)
//   complex cexp(complex *a);        exp(*a)
//   complex clog(complex *a);        log(*a)
//   complex csqrt(complex *a);       sqrt(a)
//   complex ce_itheta(Real theta);   exp(i * theta)

// The following macros work for both single and double precision,
// and also for mixed precision
// 1) Macros which appear to return values (Real or double, as appropriate):
//    cabs(*a)        Magnitude of the complex number *a
//    cabs_sq(*a)     Square of the magnitude (faster than cabs)
//    carg(*a)        Phase of the complex number *a

// 2) Macro to convert from single to double or double to single:
//    set_complex_equal(*a, *b)     b = a by components to convert

// 3) Macros for fast in-line operations:
//    (We force use of macros for complex conjugation, addition, subtraction,
//    multiplication and division)
//    CONJG(a, b)         b = conjg(a)
//    CADD(a, b, c)       c = a + b
//    CSUM(a, b)          a += b
//    CDIF(a, b)          a -= b
//    CSUB(a, b, c)       c = a - b
//    CMUL(a, b, c)       c = a * b
//    CDIV(a, b, c)       c = a / b
//    CMUL_J(a, b, c)     c = a * bdag
//    CMULJ_(a, b, c)     c = adag * b
//    CMULJJ(a, b, c)     c = (a * b)dag
//    CNEGATE(a, b)       b = -a
//    CMUL_I(a, b)        b = ia
//    CMUL_MINUS_I(a, b)  b = -ia
//    CMULREAL(a, b, c)   c = b * a, with b real and a & c complex
//    CDIVREAL(a, b, c)   c = a / b, with a & c complex and b real
//    CSUM_TPI(a, b)      a += ib, with a and b complex
//    CSUM_TMI(a, b)      a += -ib, with a and b complex
#include "../include/precision.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Rename to prevent conflict with gcc 3.xx standard functions
#define cexp cexp_single
#define clog clog_single
#define csqrt csqrt_single

// Complex number definitions
typedef struct { float real;  float imag;  } fcomplex;
typedef struct { double real; double imag; } double_complex;

#if (PRECISION==1)
#define complex fcomplex
#else
#define complex dcomplex
#endif

// Alternative name for double complex
typedef double_complex dcomplex;

// Generic precision function prototypes for complex numbers
complex cmplx(Real x, Real y);
complex cexp(complex *a);
complex clog(complex *a);
complex csqrt(complex *z);
complex ce_itheta(Real theta);

double_complex dcmplx(double x, double y);
double_complex dcexp(double_complex *a);
double_complex dclog(double_complex *a);
double_complex dcsqrt(double_complex *z);
double_complex dce_itheta(double theta);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Macros for Complex Numbers
// *b <-- *a
#define set_complex_equal(a, b) { (*b).real = (*a).real; \
                                  (*b).imag = (*a).imag; }

// sqrt(|*a|^2)
#define cabs(a) (sqrt( (*a).real * (*a).real + (*a).imag * (*a).imag))

// |*a|^2
#define cabs_sq(a) ( (*a).real * (*a).real + (*a).imag * (*a).imag)

// phase(*a)
#define carg(a) (atan2((double)(*a).imag, (double)(*a).real))

// b = adag
#define CONJG(a, b) { (b).real = (a).real; (b).imag = -(a).imag; }

// c = a + b
#define CADD(a, b, c) { (c).real = (a).real + (b).real; \
                        (c).imag = (a).imag + (b).imag; }

// a += b
#define CSUM(a, b) { (a).real += (b).real; (a).imag += (b).imag; }

// a -= b
#define CDIF(a, b) { (a).real -= (b).real; (a).imag -= (b).imag; }

// c = a - b
#define CSUB(a, b, c) { (c).real = (a).real - (b).real; \
                        (c).imag = (a).imag - (b).imag; }

// c = a * b
#define CMUL(a, b, c) { \
  (c).real = (a).real * (b).real - (a).imag * (b).imag; \
  (c).imag = (a).real * (b).imag + (a).imag * (b).real; }

// c = a / b
#define CDIV(a, b, c) { \
  double t = (b).real * (b).real + (b).imag * (b).imag; \
  (c).real = ((a).real * (b).real + (a).imag * (b).imag) / t; \
  (c).imag = ((a).imag * (b).real - (a).real * (b).imag) / t; }

// c = a * bdag
#define CMUL_J(a, b, c) { \
  (c).real = (a).real * (b).real + (a).imag * (b).imag; \
  (c).imag = (a).imag * (b).real - (a).real * (b).imag; }

// c = adag * b
#define CMULJ_(a, b, c) { \
  (c).real = (a).real * (b).real + (a).imag * (b).imag; \
  (c).imag = (a).real * (b).imag - (a).imag * (b).real; }

// c = (a * b)dag
#define CMULJJ(a, b, c) { \
  (c).real = (a).real * (b).real - (a).imag * (b).imag; \
  (c).imag = -(a).real * (b).imag - (a).imag * (b).real; }

// b = -a
#define CNEGATE(a, b) { (b).real = -(a).real; (b).imag = -(a).imag; }

// b = ia
#define CMUL_I(a, b) { (b).real = -(a).imag; (b).imag = (a).real; }

// b = -ia
#define CMUL_MINUS_I(a, b) { (b).real = (a).imag; (b).imag = -(a).real; }

// c = b * a
#define CMULREAL(a, b, c) { (c).real = (b) * (a).real; \
                            (c).imag = (b) * (a).imag; }

// c = a / b
#define CDIVREAL(a, b, c) { (c).real = (a).real / (b); \
                            (c).imag = (a).imag / (b); }

// a += ib
#define CSUM_TPI(a, b) { (a).real -= (b).imag; (a).imag += (b).real; }

// a += -ib
#define CSUM_TMI(a, b) { (a).real += (b).imag; (a).imag -= (b).real; }

#endif
// -----------------------------------------------------------------
