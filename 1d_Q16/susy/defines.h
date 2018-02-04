// -----------------------------------------------------------------
// Compiler macros common to all targets in this application
#ifndef _DEFINES_H
#define _DEFINES_H

#define SITERAND              // Use site-based random number generators
//#define TIMING              // Not currently used
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Integrator stuff
// Omelyan lambda, 2lambda and 1 - 2lambda
#define LAMBDA 0.193
#define TWO_LAMBDA 0.386
#define LAMBDA_MID 0.614
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Susy stuff
#define BMN
//#define EIG_POT           // Single rather than double-trace scalar potential
                          // Single-trace case is minimally tested!!!
#define SV                // Site/vector terms in action
#define PUREGAUGE
//#define DEBUG_CHECK       // Print lambdas, offsets, etc.
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measurement stuff
// Threshold to print warning about non-zero imaginary components
// of quantities expected to be real
#define IMAG_TOL 1.0e-8
#define SQ_TOL 1.0e-16

// Maximum time value and spatial distance for Wilson loops
#define MAX_T (nt / 2)
#define MAX_X (nx / 2 - 1)
#endif
// -----------------------------------------------------------------
