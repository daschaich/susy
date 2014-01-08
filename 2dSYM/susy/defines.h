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

// RHMC degree -- set Norder at runtime so it can be reset (e.g., to 1)
#define DEGREE 15
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Temporal boundary condition stuff
#define OPP_LDIR(Ldir) (3 - (Ldir))   // Same as OPP_DIR
#define PBC -1.0    // 1.0 for PBC, -1.0 for APBC
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Susy stuff
#define SV                  // Site/vector terms in action
#define VP                  // Vector/plaquette terms in action
#define DET                 // Determinant term in the action
//#define DEBUG_CHECK         // Print lambdas, offsets, etc.

#define NUMLINK 2           // Number of A2* links per site, and offsets
#define NUMGEN (NCOL*NCOL)  // Number of U(N) generators
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measurement stuff
// Threshold to print warning about non-zero imaginary components
// of quantities expected to be real
#define IMAG_TOL 1e-16

// Maximum time value and spatial distance for Wilson loops
#define MAX_T (nt / 2)
#define MAX_X (nx / 2 - 1)

// Calculate eigenvalues of Ddag.D (as opposed to untested D eigenvalues)
#define DDdag
#endif
// -----------------------------------------------------------------
