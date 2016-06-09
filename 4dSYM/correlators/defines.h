// -----------------------------------------------------------------
// Compiler macros common to all targets in this application
#ifndef _DEFINES_H
#define _DEFINES_H

//#define TIMING              // Not currently used
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Susy stuff
//#define DEBUG_CHECK       // Print lambdas, offsets, etc.

// Whether or not to project determinant out of given observable
#define NODET 0
#define YESDET 1
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
