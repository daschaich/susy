// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;           // 1 if it is time to stop

  // Initialization parameters
  int nx, ny, nz, nt;     // Lattice dimensions

  int startflag;          // What to do for beginning lattice
  int fixflag;            // Whether to gauge fix to Coulomb gauge
  char startfile[MAXFILENAME];

  double vevK[N_K], vevS[N_K];  // Konishi and SUGRA vacuum subtractions

#ifdef SMEAR
  // Parameters for APE or stout smearing
  int smearflag;                // NONE, STOUT, APE
  int Nsmear;
  Real alpha;
#endif
} params;
#endif
// -----------------------------------------------------------------
