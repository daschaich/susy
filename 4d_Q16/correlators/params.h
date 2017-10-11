// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;           // 1 if it is time to stop

  // Initialization parameters
  int nx, ny, nz, nt;     // Lattice dimensions

  // Use Nblock blocks each containing Nmeas measurements
  int Nblock, Nmeas;
  char cfg[MAX_CFG][MAXFILENAME];   // Lattice file names

#ifdef SMEAR
  // Parameters for APE or stout smearing
  int smearflag;                    // NONE, STOUT, APE
  int Nsmear;
  Real alpha;
#endif
} params;
#endif
// -----------------------------------------------------------------
