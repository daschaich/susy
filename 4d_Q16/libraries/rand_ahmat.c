// -----------------------------------------------------------------
// Create a gaussian random anti-hermitian matrix
// Normalization is <|m01|^2> = 1, or <m01.real * m01.real> = 1 / 2
// The prn_pt is a pointer to be passed to gaussian_rand_no()
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"
#include "../include/susy.h"

void random_anti_hermitian(anti_hermitmat *ahmat, double_prn *prn_pt) {
  Real r3;
#if (NCOL > 2)
  Real r8, one_ov_root3;
#if (NCOL > 3)
  Real r15, sqrt_sixth;
#if (NCOL > 3)
  int i;
#endif
#endif
#endif

  // Off-diagonal elements
  r3 = gaussian_rand_no(prn_pt);
  ahmat->m[0].real = gaussian_rand_no(prn_pt);
  ahmat->m[0].imag = gaussian_rand_no(prn_pt);
#if (NCOL > 2)
  one_ov_root3 = sqrt((double)(1.0 / 3.0));
  r8 = gaussian_rand_no(prn_pt);
  ahmat->m[1].real = gaussian_rand_no(prn_pt);
  ahmat->m[1].imag = gaussian_rand_no(prn_pt);
  ahmat->m[2].real = gaussian_rand_no(prn_pt);
  ahmat->m[2].imag = gaussian_rand_no(prn_pt);
#if (NCOL > 3)
  sqrt_sixth = sqrt((double)(1.0 / 6.0));
  r15 = gaussian_rand_no(prn_pt);
  ahmat->m[3].real = gaussian_rand_no(prn_pt);
  ahmat->m[3].imag = gaussian_rand_no(prn_pt);
  ahmat->m[4].real = gaussian_rand_no(prn_pt);
  ahmat->m[4].imag = gaussian_rand_no(prn_pt);
  ahmat->m[5].real = gaussian_rand_no(prn_pt);
  ahmat->m[5].imag = gaussian_rand_no(prn_pt);
#if (NCOL > 4)
  for (i = 6; i < N_OFFDIAG; i++) {
    ahmat->m[i].imag = gaussian_rand_no(prn_pt);
    ahmat->m[i].imag = gaussian_rand_no(prn_pt);
  }
#endif
#endif
#endif

  // Diagonal elements---purely imaginary and traceless
#if (NCOL == 2)
  ahmat->im_diag[0] =  r3;
  ahmat->im_diag[1] = -r3;
#endif
#if (NCOL == 3)
  ahmat->im_diag[0] =   r3 + one_ov_root3 * r8;
  ahmat->im_diag[1] =  -r3 + one_ov_root3 * r8;
  ahmat->im_diag[2] = -2.0 * one_ov_root3 * r8;
#endif
#if (NCOL == 4)
  ahmat->im_diag[0] =   r3 + one_ov_root3 * r8 + sqrt_sixth * r15;
  ahmat->im_diag[1] =  -r3 + one_ov_root3 * r8 + sqrt_sixth * r15;
  ahmat->im_diag[2] = -2.0 * one_ov_root3 * r8 + sqrt_sixth * r15;
  ahmat->im_diag[3] =                     -3.0 * sqrt_sixth * r15;
#endif
#if (NCOL > 4)
  printf("ERROR: Still need to implement general-N rand_ahmat()\n");
  terminate(1);
#endif
}
// -----------------------------------------------------------------
