// -----------------------------------------------------------------
// Create a gaussian random anti-hermitian matrix
// Normalization is <|m01|^2> = 1, or <m01.real * m01.real> = 1 / 2
// The prn_pt is a pointer to be passed to gaussian_rand_no()
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"
#include "../include/susy.h"

#ifdef RANDAH_DEBUG
#include <stdio.h>
#endif

void random_anti_hermitian(anti_hermitmat *ahmat, double_prn *prn_pt) {
  Real r3;
#if (NCOL > 2)
  Real r8;
#if (NCOL > 3)
  Real r15;
#if (NCOL > 4)
  int i, j;
  Real tr;
#endif
#endif
#endif

  // Off-diagonal elements
  r3 = gaussian_rand_no(prn_pt);
  ahmat->m[0].real = gaussian_rand_no(prn_pt);
  ahmat->m[0].imag = gaussian_rand_no(prn_pt);
#if (NCOL > 2)
  r8 = gaussian_rand_no(prn_pt);
  r8 *= sqrt((double)(1.0 / 3.0));
  ahmat->m[1].real = gaussian_rand_no(prn_pt);
  ahmat->m[1].imag = gaussian_rand_no(prn_pt);
  ahmat->m[2].real = gaussian_rand_no(prn_pt);
  ahmat->m[2].imag = gaussian_rand_no(prn_pt);
#if (NCOL > 3)
  r15 = gaussian_rand_no(prn_pt);
  r15 *= sqrt((double)(1.0 / 6.0));
  ahmat->m[3].real = gaussian_rand_no(prn_pt);
  ahmat->m[3].imag = gaussian_rand_no(prn_pt);
  ahmat->m[4].real = gaussian_rand_no(prn_pt);
  ahmat->m[4].imag = gaussian_rand_no(prn_pt);
  ahmat->m[5].real = gaussian_rand_no(prn_pt);
  ahmat->m[5].imag = gaussian_rand_no(prn_pt);
#if (NCOL > 4)
  for (i = 6; i < N_OFFDIAG; i++) {
    ahmat->m[i].real = gaussian_rand_no(prn_pt);
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
  ahmat->im_diag[0] =   r3 + r8;
  ahmat->im_diag[1] =  -r3 + r8;
  ahmat->im_diag[2] = -2.0 * r8;
#endif
#if (NCOL > 3)
  ahmat->im_diag[0] =   r3 + r8 + r15;
  ahmat->im_diag[1] =  -r3 + r8 + r15;
  ahmat->im_diag[2] = -2.0 * r8 + r15;
  ahmat->im_diag[3] =      -3.0 * r15;
#endif
#if (NCOL > 4)
  for (i = 4; i < NCOL; i++) {
    // Reuse r15 for next properly scaled random number
    // Scaling factor is sqrt(2 / (i * (i + 1)))
    // Check: Previously sqrt of 1, 1/3, 1/6 for i=1, 2, 3
    r15 = gaussian_rand_no(prn_pt);
    tr = (Real)i * (i + 1);
    r15 *= sqrt((double)(1.0 / tr));

    // Add new r15 to i existing elements
    for (j = 0; j < i; j++)
      ahmat->im_diag[j] += r15;

    // Initialize next element with -i * r15
    ahmat->im_diag[i] = -1.0 * i * r15;
  }
#endif

#ifdef RANDAH_DEBUG
  // Check tracelessness
  int k;
  r3 = ahmat->im_diag[0];
  for (k = 1; k < NCOL; k++)
    r3 += ahmat->im_diag[k];
  if (fabs(r3) > 1e-6)
    printf("Warning: random_anti_hermitian trace = %.4g\n", r3);
#endif
}
// -----------------------------------------------------------------
