// -----------------------------------------------------------------
// Create gaussian random anti-hermitian matrix
// Normalization is <|m01|^2> = 1, <m01.real^2> = 0.5
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"

// prn_pt is a pointer to be passed to gaussian_rand_no()
void random_anti_hermitian(anti_hermitmat *dest, double_prn *prn_pt) {
  Real r3;
#if (NCOL > 2)
  Real r8, sqrt_third;
#endif
#if (NCOL > 3)
  Real r15, sqrt_sixth;
#endif

  // Off-diagonal elements
  r3 = gaussian_rand_no(prn_pt);
  dest->m01.real = gaussian_rand_no(prn_pt);
  dest->m01.imag = gaussian_rand_no(prn_pt);
#if (NCOL > 2)
  sqrt_third = sqrt((double)(1.0 / 3.0));
  r8 = gaussian_rand_no(prn_pt);
  dest->m02.real = gaussian_rand_no(prn_pt);
  dest->m12.real = gaussian_rand_no(prn_pt);
  dest->m02.imag = gaussian_rand_no(prn_pt);
  dest->m12.imag = gaussian_rand_no(prn_pt);
#endif
#if (NCOL > 3)
  sqrt_sixth = sqrt((double)(1.0 / 6.0));
  r15 = gaussian_rand_no(prn_pt);
  dest->m03.real = gaussian_rand_no(prn_pt);
  dest->m13.real = gaussian_rand_no(prn_pt);
  dest->m23.real = gaussian_rand_no(prn_pt);
  dest->m03.imag = gaussian_rand_no(prn_pt);
  dest->m13.imag = gaussian_rand_no(prn_pt);
  dest->m23.imag = gaussian_rand_no(prn_pt);
#endif

  // Diagonal elements
#if (NCOL==2)
  dest->m00im =  r3;
  dest->m11im = -r3;
#endif
#if (NCOL==3)
  dest->m00im =   r3 + sqrt_third * r8;
  dest->m11im =  -r3 + sqrt_third * r8;
  dest->m22im = -2.0 * sqrt_third * r8;
#endif
#if (NCOL==4)
  dest->m00im =   r3 + sqrt_third * r8 + sqrt_sixth * r15;
  dest->m11im =  -r3 + sqrt_third * r8 + sqrt_sixth * r15;
  dest->m22im = -2.0 * sqrt_third * r8 + sqrt_sixth * r15;
  dest->m33im =                   -3.0 * sqrt_sixth * r15;
#endif
#if (NCOL > 4)
  node0_printf("random_anti_hermitian only works for NCOL<=4!");
#endif
}
// -----------------------------------------------------------------
