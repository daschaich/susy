// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;           // 1 if it is time to stop

  // Initialization parameters
  int nx, ny, nz, nt;     // Lattice dimensions
  int PBC;                // Temporal fermion boundary condition
  int iseed;              // For random numbers

  // RHMC and multi-mass CG parameters
  // Number of Nth roots and polynomial order
  int Nroot, Norder;

  int warms;              // The number of warmup trajectories
  int trajecs;            // The number of real trajectories
  Real traj_length;       // The length of each trajectory
  int nsteps[2];          // Fermion and gauge steps
  int propinterval;       // Number of trajectories between measurements
  int startflag;          // What to do for beginning lattice
  int fixflag;            // Whether to gauge fix to Coulomb gauge
  int saveflag;           // What to do with lattice at end
  Real lambda, kappa;     // 't Hooft coupling and Nc/(2lambda)
  Real bmass, fmass;      // Bosonic and fermion masses
  Real kappa_u1;          // Plaquette determinant coupling
  Real G;                 // Q-invariant plaquette determinant coupling

  // Inversion parameters
  int niter;                    // Maximum number of CG iterations
  Real rsqmin;                  // For deciding on convergence
  char startfile[MAXFILENAME], savefile[MAXFILENAME];

#ifdef BILIN
  int nsrc;                     // Number of stochastic sources
#endif

#ifdef SMEAR
  // Parameters for APE or stout smearing
  int smearflag;                // NONE, STOUT, APE
  int Nsmear;
  Real alpha;
#endif

#ifdef EIG
  // Eigenvalue parameters
  int Nvec, maxIter;
  Real eig_tol;
#endif

#ifdef PHASE
  // Pfaffian parameters
  int ckpt_load, ckpt_save;
#endif
} params;
#endif
// -----------------------------------------------------------------
