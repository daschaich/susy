// -----------------------------------------------------------------
// Code for making rational functions for N=4 SYM
// Based on poly4.C by Doug Toussaint

// Usage: poly < [params_in_file] > [params_rhmc_file]

// Where params_in_file contains
//    Nroot  MD_order  action_order  eval_min  eval_max  precision
// Example: 2  15  15  0.0000001  1500  65

// Calculate three rational approximations, for molecular dynamics
// evolution, gaussian random heatbath update, and fermion action computation
// Molecular dynamics: PROD D^(-1 / (4 * Nroot))
// Heatbath ("GR"):    PROD D^(1 / (8 * Nroot))
// Fermion Action:     PROD D^(-1 / (8 * Nroot))

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include"alg_remez.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void print_check_ratfunc(int sign, int Nroot4, double lambda_low, int order,
                         double norm,  double *res, double *pole, char *tag) {

  int i;
  double x, sum, f1;
  if (strcmp(tag, "OMIT") == 0)
    return;

  printf("\n\n# Rational function for %s\n", tag);
  printf("order_%s %d\n\n", tag, order);
  printf("res_%s NORM  %18.16e \n", tag, norm);
  for (i = 0; i < order; i++)
    printf("res_%s %d  %18.16e\n", tag, i, res[i]);
  printf("\n");
  printf("pole_%s 99.9\n", tag); //DUMMY!!
  for (i = 0; i < order; i++) {
    printf("pole_%s %d %18.16e\n", tag, i, pole[i]);
  }
  printf("\n");

  // Check - compute the function at the low endpoint of the interval
  int ii;
  x = lambda_low;
  sum = norm;
  for (i = 0; i < order; i++)
    sum += res[i] / (x + pole[i]);

  f1 = pow(x, ((double)sign / Nroot4));

  printf("# CHECK: f(%e) = %e = %e?\n", x, sum, f1);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int work(int Nroot, int order, double lambda_low,  double lambda_high,
         int precision, char *tag1, char *tag2) {

  int Nroot4 = 4 * Nroot;
  double *res = new double[order];
  double *pole = new double[order];
  double error;   // The error from the approximation
                  // The relative error is minimised
                  // If another error minimisation is required,
                  // change line 398 in alg_remez.C

  // The partial fraction expansion takes the form
  // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
  double norm;
  double bulk = exp(0.5 * (log(lambda_low) + log(lambda_high)));

  // Instantiate the Remez class
  AlgRemez remez(lambda_low, lambda_high, precision);

  // Generate the required approximation
  fprintf(stderr, "Generating a (%d,%d) rational function ", order, order);
  fprintf(stderr, "using %d digit precision\n", precision);
  error = remez.generateApprox(order, order, 1, Nroot4, 0.0, 0, -1, -1.0,
                                             0, -1,    -1.0, 0, -1, -1.0);

  // Find the partial fraction approximation to the function x^{y / z}
  // This only works currently for the special case that n = d
  remez.getPFE(res, pole, &norm);

  print_check_ratfunc(1, Nroot4, lambda_low, order, norm, res, pole, tag1);

  // Find pfe of the inverse function
  remez.getIPFE(res, pole, &norm);

  print_check_ratfunc(-1, Nroot4, lambda_low, order, norm, res, pole, tag2);

  FILE *error_file = fopen("error.dat", "w");
  double x, f, r;
  for (x = lambda_low; x < lambda_high; x *= 1.01) {
    f = remez.evaluateFunc(x);
    r = remez.evaluateApprox(x);
    fprintf(error_file, "%e %e\n", x,  (r - f) / f);
  }
  fclose(error_file);

  delete res;
  delete pole;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main (int argc, char* argv[]) {
  char *tag1, *tag2;
  int i = 0, Nroot;
  int order1;                     // The degree of the numerator
  int order2;                     // The degree of the denominator
  int precision;                  // The precision that gmp uses
  double lambda_low, lambda_high; // The bounds of the approximation

  // Read the number of roots
  scanf("%d", &Nroot);
  printf("Nroot %d\n\n", Nroot);

  // Set the required degrees of approximation
  scanf("%d", &order1);
  scanf("%d", &order2);

  // Set the approximation bounds
  scanf("%le", &lambda_low);
  scanf("%le", &lambda_high);

  // Set the precision of the arithmetic
  scanf("%d", &precision);

  // For the MD term we need only the inverse
  tag1 = (char *)"OMIT";
  tag2 = (char *)"MD";
  work(Nroot, order1, lambda_low, lambda_high, precision, tag1, tag2);

  // The random source term takes the function
  // The action term takes the inverse
  Nroot *= 2;
  tag1 = (char *)"GR";
  tag2 = (char *)"FA";
  work(Nroot, order2, lambda_low, lambda_high, precision, tag1, tag2);
}
// -----------------------------------------------------------------
