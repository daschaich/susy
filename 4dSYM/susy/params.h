// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;       // 1 if it is time to stop

  // Initialization parameters
  int nx, ny, nz, nt; // Lattice dimensions
  int iseed;          // For random numbers

  int warms;                // The number of warmup trajectories
  int trajecs;              // The number of real trajectories
  Real traj_length;         // The length of each trajectory
  int nsteps[2];            // Fermion and gauge steps
  int propinterval;         // Number of trajectories between measurements
  int startflag;            // What to do for beginning lattice
  int saveflag;             // What to do with lattice at end
  Real lambda, kappa;       // 't Hooft coupling and Nc/(2lambda)
  Real bmass, fmass;        // Bosonic and fermion masses
  Real kappa_u1;            // Determinant coupling
  Real Ckonishi;            // Konishi coupling
  int fixflag;              // Whether to gauge fix to Coulomb gauge

  // Inversion parameters
  int niter;                    // Maximum number of CG iterations
  Real rsqmin;                  // For deciding on convergence
  int nsrc;                     // Number of stochastic sources
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  // ILDG LFN if applicable

#ifdef EIG
  // Eigenvalue parameters
  int Nvec, maxIter;
  Real eig_tol;
#endif

#ifdef STOUT
  // Stout-smearing parameters
  int Nstout;
  Real rho;
#endif

#ifdef PHASE
  // Pfaffian parameters
  int ckpt_load, ckpt_save;
#endif
} params;
#endif
// -----------------------------------------------------------------
