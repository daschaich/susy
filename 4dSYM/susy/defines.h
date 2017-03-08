// -----------------------------------------------------------------
// Compiler macros common to all targets in this application
#ifndef _DEFINES_H
#define _DEFINES_H

#define SITERAND              // Use site-based random number generators
#define GAUGE_FIX_TOL 1.0e-7  // For gauge fixing
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
//#define HYBRID            // Truncate the gauge links to live in SL(N,C)
//#define EIG_POT           // Single rather than double-trace scalar potential
                          // Single-trace case is not actively tested!!!
#define LINEAR_DET        // det-1 rather than |det-1|^2
                          // |det-1|^2 is not actively tested!!!
#define SV                // Site/vector terms in action
#define VP                // Vector/plaquette terms in action
#define QCLOSED           // Q-closed terms in action
//#define DEBUG_CHECK       // Print lambdas, offsets, etc.

#define NQLINK 28         // Number of offsets for Q-closed terms
#define NTERMS 30         // Number of DbmP, DbpP F1Q and F2Q terms

// Tunable parameter in gauge action
#define C2 1.0

// Whether or not to project determinant out of given observable
#define NODET 0
#define YESDET 1
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measurement stuff
// Gauge fixing options
#define NO_GAUGE_FIX 30
#define COULOMB_GAUGE_FIX 31

// Smearing options
#define NO_SMEAR 40
#define STOUT_SMEAR 41
#define APE_SMEAR 42

// Threshold to print warning about non-zero imaginary components
// of quantities expected to be real
#define IMAG_TOL 1.0e-8
#define SQ_TOL 1.0e-16

// Maximum time value and spatial distance for Wilson loops
#define MAX_T (nt / 2)
#define MAX_X (nx / 2 - 1)
#endif
// -----------------------------------------------------------------
