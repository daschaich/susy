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
  Real friction;          // Parameter used for SMD algorithm
  int nsteps[2];          // Fermion and gauge steps
  int propinterval;       // Number of trajectories between measurements
  int startflag;          // What to do for beginning lattice
  int fixflag;            // Whether to gauge fix to Coulomb gauge
  int saveflag;           // What to do with lattice at end
  Real lambda, kappa;     // 't Hooft coupling and Nc/(4lambda)
  Real bmass, fmass;      // Bosonic and fermion masses
  Real kappa_u1;          // Plaquette determinant coupling
  Real G;                 // Q-invariant plaquette determinant coupling
#ifdef DIMREDUCE
  Real cWline;            // Coefficient of center-breaking term protecting
                          // single-link 'Wilson line' in reduced dir(s)
#endif

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

#ifdef CHEB
  // Chebyshev spectral density stuff
  int Nstoch;
  int cheb_order;               // Number of coefficients to compute
  Real lambda_min, lambda_max;  // Bounds on spectral range
#endif

#ifdef MODE
  // Giusti--Luescher stochastic mode number and step function stuff
#define MAX_OMEGA 100
  int Nstoch;
  int step_order;           // Selects between options in mode_coeffs.c
  int numOmega;             // Number of Omega at which to evaluate nu
  Real Omega[MAX_OMEGA];    // List of Omega at which to evaluate nu
#endif

#ifdef PHASE
  // Pfaffian parameters
  int ckpt_load, ckpt_save;
#endif
} params;
#endif
// -----------------------------------------------------------------
