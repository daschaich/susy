// -----------------------------------------------------------------
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <math.h>

#define NMAX 1000

typedef struct {
  double epsilon;
  int order;
  double *c;
  double alpha;
} error2_pars_type;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Chebyshev polynomials
// T[i, z] = 2zT[i - 1, z] - T[i - 2, z]
static void Chebyshev(double *ret, int n, double z) {
  ret[0] = 1.0;
  ret[1] = z;
  for (int i = 2; i <= n; i++)
    ret[i] = 2.0 * z * ret[i - 1] - ret[i - 2];
}

// Derivative of Chebyshev polynomials
// T'[i, z] = 2T[i - 1, z] + 2zT'[i - 1, z] - T'[i - 2, z]
static void DerChebyshev(double *der, double *T, int n, double z) {
  Chebyshev(T, n, z);
  der[0] = 0.0;
  der[1] = 1.0;
  for (int i = 2; i <= n; i++)
    der[i] = 2.0 * T[i - 1] + 2.0 * z * der[i - 1] - der[i - 2];
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Approximating polynomial
static double P(int n, double *c, double y, double epsilon) {
  int i;
  double z = (2. * y - 1. - epsilon) / (1. - epsilon);
  double T[n + 1];
  double ret = 0.;

  Chebyshev(T, n, z);
  for (i = 0; i <= n; i++)
    ret += c[i]*T[i];

  return ret;
}

// Derivative of the approximating polynomial
static double derP(int n, double *c, double y, double epsilon) {
  int i;
  double z = (2.0 * y - 1.0 - epsilon) / (1.0 - epsilon);
  double T[n + 1], derT[n + 1];
  double ret = 0.0;

  DerChebyshev(derT, T, n, z);
  for (i = 0; i <= n; i++)
    ret += c[i] * derT[i];

  return (ret * 2.0 / (1.0 - epsilon));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Error function
static double errh(int n, double *c, double y, double epsilon) {
  return (1.0 - sqrt(y) * P(n, c, y, epsilon));
}


// Derivative of the error function
static double dererrh(int n, double *c, double y, double epsilon) {
  double p = P(n, c, y, epsilon), dp = derP(n, c, y, epsilon);
  return -0.5 * p / sqrt(y) - sqrt(y) * dp;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Solve [ sqrt(y[m]) * P(c, y[m]) + (-1)^m u == 1. , {c, u} ]
static void solvesystem(double *c, double *u, int n, double *y,
                        double epsilon) {

  int i, m;
  double T[n + 1], z, sy, sign = 1.0, X;
  gsl_matrix *A = gsl_matrix_alloc(n + 2, n + 2);
  gsl_vector *b = gsl_vector_alloc(n + 2);
  gsl_vector *x = gsl_vector_alloc(n + 2);

  for (m = 0; m <= n + 1; m++) {
    z = (2.0 * y[m] - 1.0 - epsilon) / (1.0 - epsilon);
    sy = sqrt(y[m]);
    Chebyshev(T, n, z);
    for (i = 0; i <= n; i++)
      gsl_matrix_set(A, m, i, sy * T[i]);

    gsl_matrix_set(A, m, n + 1, sign);
    sign = -sign;
    gsl_vector_set(b, m, 1.0);
  }

  gsl_linalg_HH_solve(A, b, x);
  for (i = 0; i <= n; i++)
    c[i] = gsl_vector_get(x, i);

  *u = gsl_vector_get(x, n + 1);

  for (m = 0; m <= n + 1; m++) {
    z = (2.0 * y[m] - 1.0 - epsilon) / (1.0 - epsilon);
    sy = sqrt(y[m]);
    Chebyshev(T, n, z);
    X = 1.0;
    for (i = 0; i <= n; i++)
      X = X - c[i] * sy * T[i];
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Zero of the (derivative of the) error function
// Assuming that exactly one zero exists
static double zero(double (*fptr)(int, double*, double, double), int n,
                   double *c, double epsilon, double a, double b,
                   double prec) {

  int counter = 0;
  double h1, h2, hmid, mid;
  double prec2 = prec*prec;

  h1 = fptr(n, c, a, epsilon);
  h2 = fptr(n, c, b, epsilon);
  if (h1 * h1 < prec2)
    return a;
  if (h2 * h2 < prec2)
    return b;
  while(1) {
    mid = (a + b) / 2.0;
    counter++;
    hmid = fptr(n, c, mid, epsilon);
    if (hmid * hmid < prec2)
      return mid;
    if (h1 * hmid > 0.0) {
      a = mid;
      h1 = hmid;
    }
    else {
      b = mid;
      h2 = hmid;
    }
  }
  return 0.0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
static int squareroot(double delta, double epsilon, double *c, int *order,
                      double *err) {

  int i, m, n = 4, counter;
  double y[NMAX + 2];
  double zeroes[NMAX + 1];
  double u = 0., error = 1., tmp;

  while (error > delta) {
    n++;
    error = 1.0;
    u = 0.0;
    counter = 0;

    // Initial guess
    for (m = 0; m <= n + 1; m++) {
      y[m] = 1.0 + epsilon - (1.0 - epsilon) * cos(m * M_PI / (n + 1.0));
      y[m] *= 0.5;
    }

    while (error - fabs(u) > fabs(u) * .01) {
      solvesystem(c, &u, n, y, epsilon);

      for (i = 0; i <= n; i++)
        zeroes[i] = zero(&errh, n, c, epsilon, y[i], y[i + 1], u / 1000.0);

      error = fabs(errh(n, c, y[0], epsilon));
      for (i = 0; i < n; i++) {
        y[i + 1] = zero(&dererrh, n, c, epsilon, zeroes[i], zeroes[i + 1],
                        u / (zeroes[i + 1] - zeroes[i]) / 1000.0);
        tmp = fabs(errh(n, c, y[i + 1], epsilon));
        if (tmp > error)
          error = tmp;
      }
      tmp = fabs(errh(n, c, y[n + 1], epsilon));
      if (tmp > error)
        error = tmp;
      counter++;
    }
  }
  *order = n;
  *err = error;
  return counter;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
static double h(double x, double epsilon, int order, double *c) {
  int n;
  double b0 = 0, b1 = 0, b2;
  double z = (2.0 * x * x - 1.0 - epsilon) / (1.0 - epsilon);

  for (n = order; n >= 0; n--) {
    b2 = b1;
    b1 = b0;
    b0 = c[n] + 2.0 * z * b1 - b2;
  }
  return 0.5 - 0.5 * x * (b0 - b1 * z);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
static double error2(double x, void *pars) {
  double al = ((error2_pars_type*)pars)->alpha;
  double hv = h(x, ((error2_pars_type*)pars)->epsilon,
                   ((error2_pars_type*)pars)->order,
                   ((error2_pars_type*)pars)->c);

  // A**B = exp[B log A]
  double num_pow = al + 1.0,         num_arg = 1.0 + x;
  double den_pow = (3.0 + al) / 2.0, den_arg = 1.0 - (x * x);
  double num = gsl_sf_exp(num_pow * gsl_sf_log(num_arg));
  double den = gsl_sf_exp(den_pow * gsl_sf_log(den_arg));
  return (hv * hv * hv * hv * num / den);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Returns ratio (Omega / Omega_*)^(2(alpha + 1))
double star(double delta, double epsilon, int order, double *c,
            double alpha) {

  double pow = alpha + 1.0;
  double intres, interr, ret, temp, seps = sqrt(epsilon);
  error2_pars_type error2_pars;
  error2_pars.epsilon = epsilon;
  error2_pars.order = order;
  error2_pars.c = c;
  error2_pars.alpha = alpha;

  gsl_function gsl_error2 = {&error2, &error2_pars};
  gsl_integration_workspace *gsl_ws_int
                            = gsl_integration_workspace_alloc(1000000);

  gsl_integration_qag(&gsl_error2, -seps, seps, 1.e-2, 0, 1000000,
                      GSL_INTEG_GAUSS15, gsl_ws_int, &intres, &interr);

  temp = gsl_sf_exp(pow * gsl_sf_log((1.0 - seps) / (1.0 + seps)) / 2.0);
  ret = temp + pow * intres;
  temp = gsl_sf_exp(gsl_sf_log(ret) / pow);     // A**B = exp[B log A]
  ret = temp * temp;                            // Squared!
  gsl_integration_workspace_free(gsl_ws_int);

  return ret;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Main program prints out coefficients and other pertinent information
int main(int argc, char* argv[]) {
  int i, order;
  double alpha, epsilon, delta, err, c[NMAX + 1];

  if (argc != 3) {
    printf("Usage: %s epsilon delta\n", argv[0]);
    return -1;
  }

  epsilon = atof(argv[1]);
  delta = atof(argv[2]);
  squareroot(delta, epsilon, c, &order, &err);

//  for(i = 0; i < 5001;i++) {
//    double x = 2. * i / 5000. - 1.;
//    printf("%e %e SIGN\n", x, h(x));
//  }

  printf("order = %d;\n", order);
  printf("epsilon = %.8g;\n", epsilon);
  printf("err = %.8g;\n\n", err);
  for (alpha = 0; alpha <= 3.0; alpha += 0.5) {
    printf("star(%.2g) = %.8g;\n", alpha,
           star(delta, epsilon, order, c, alpha));
  }
  printf("\n");
  for (i = 0; i <= order; i++)
    printf("coeffs[%d] = %.16g;\n", i, c[i]);

  return 0;
}
// -----------------------------------------------------------------
