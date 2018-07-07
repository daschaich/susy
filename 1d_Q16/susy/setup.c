// -----------------------------------------------------------------
// BFSS/BMN setup
#include "susy_includes.h"

#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size and seed, and send to others
int initial_set() {
  int prompt = 0, status = 0;
  if (mynode() == 0) {
    // Print banner
#ifdef BMN
    printf("BMN, Nc = %d, DIMF = %d, fermion rep = adjoint\n",
           NCOL, DIMF);
#else
    printf("BFSS, Nc = %d, DIMF = %d, fermion rep = adjoint\n",
           NCOL, DIMF);
#endif
    printf("Microcanonical simulation with refreshing\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
#ifdef HMC_ALGORITHM
    printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef PHI_ALGORITHM
    printf("Phi algorithm\n");
#else   // Quit!
    printf("Only works for phi algorithm\n");
    exit(1);
#endif
    time_stamp("start");
    status = get_prompt(stdin, &prompt);

    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);
    IF_OK status += get_i(stdin, prompt, "PBC", &par_buf.PBC);
    IF_OK status += get_i(stdin, prompt, "iseed", &par_buf.iseed);

    // Number of Nth roots to take (in addition to 1/4 power)
    IF_OK status += get_i(stdin, prompt, "Nroot", &par_buf.Nroot);

    // RHMC degree
    IF_OK status += get_i(stdin, prompt, "Norder", &par_buf.Norder);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  nt = par_buf.nt;
  PBC = par_buf.PBC;
  iseed = par_buf.iseed;

  // Lattice volume sanity checks, including dimensional reduction
  if (mynode() == 0) {
    if (nt < 1) {
      printf("nt must be positive\n");
      exit(1);
    }
  }

  // Set up stuff for RHMC and multi-mass CG
  Nroot = par_buf.Nroot;
  fnorm = malloc(sizeof fnorm * Nroot);
  max_ff = malloc(sizeof max_ff * Nroot);

  Norder = par_buf.Norder;
  amp = malloc(sizeof amp * Norder);
  amp4 = malloc(sizeof amp4 * Norder);
  amp8 = malloc(sizeof amp8 * Norder);
  shift = malloc(sizeof shift * Norder);
  shift4 = malloc(sizeof shift4 * Norder);
  shift8 = malloc(sizeof shift8 * Norder);

  this_node = mynode();
  number_of_nodes = numnodes();
  one_ov_N = 1.0 / (Real)NCOL;
  total_iters = 0;
  minus1 = cmplx(-1.0, 0.0);
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate space for fields
void make_fields() {
  Real size = (Real)(2.0 * sizeof(complex));
  FIELD_ALLOC(tr_eta, complex);
  FIELD_ALLOC(tr_dest, complex);

  size += (Real)(2.0 * NFERMION * sizeof(matrix));
  FIELD_ALLOC_VEC(src, matrix, NFERMION);
  FIELD_ALLOC_VEC(dest, matrix, NFERMION);
  FIELD_ALLOC_MAT(Gamma_X, matrix, NCHIRAL_FERMION, NCHIRAL_FERMION);

  // CG Fermions
  size += (Real)(3.0 * NFERMION * sizeof(matrix));
  FIELD_ALLOC_VEC(mpm, matrix, NFERMION);
  FIELD_ALLOC_VEC(pm0, matrix, NFERMION);
  FIELD_ALLOC_VEC(rm, matrix, NFERMION);

  // Temporary matrices and Fermions
  size += (Real)((2.0 + NFERMION + NSCALAR) * sizeof(matrix));
  FIELD_ALLOC(tempmat, matrix);
  FIELD_ALLOC(tempmat2, matrix);
  FIELD_ALLOC_VEC(temp_ferm, matrix, NFERMION);
  FIELD_ALLOC_VEC(temp_X, matrix, NSCALAR);

#if defined(EIG) || defined(PHASE)
  size += (Real)(2.0 * NFERMION * sizeof(matrix));
  FIELD_ALLOC_VEC(src, matrix, NFERMION);
  FIELD_ALLOC_VEC(res, matrix, NFERMION);
#endif

  size *= sites_on_node;
  node0_printf("Mallocing %.1f MBytes per core for fields\n", size / 1e6);
#ifdef PHASE
  // Total number of matvecs is (nt * NFERMION * DIMF)^2 / 4
  Nmatvecs = nt * NFERMION * DIMF * nt * DIMF;

  // Total size of matrix is (nt * NFERMION * DIMF) x (sites_on_node * NFERMION * DIMF)
  size = (Real)(nt * NFERMION * DIMF * NFERMION * DIMF * sizeof(complex));
  size *= sites_on_node;
  node0_printf("Q has %d columns --> %li matvecs and %.1f MBytes per core...",
               nt * NFERMION * DIMF, Nmatvecs, size / 1e6);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_bc() {
  register int i;
  register site *s;

  // Single-offset terms only
  FORALLSITES(i, s) {
    s->bc[0] = 1.0;
    s->bc[1] = 1.0;
    if (s->t + 1 > nt - 1)
      s->bc[0] = PBC;
    if (s->t - 1 < 0)
      s->bc[1] = PBC;
  }

  // BC test
  //  FORALLSITES(i, s)
  //      printf("%d : %4.2g %4.2g\n", s->t, s->bc[0], s->bc[1]);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume and seed
  prompt = initial_set();
  // Initialize the node random number generator
  initialize_prn(&node_prn, iseed, nt + mynode());
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();

  // Set up boundary conditions
  if (PBC >= 0)
    node0_printf("Periodic temporal boundary conditions\n");
  if (PBC < 0)
    node0_printf("Antiperiodic temporal boundary conditions\n");
  setup_bc();

  // Allocate space for fields
  make_fields();

  return prompt;
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Read in parameters for Monte Carlo
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;
#ifdef EIG
  int i, j;
#endif
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Warms, trajecs
    IF_OK status += get_i(stdin, prompt, "warms", &par_buf.warms);
    IF_OK status += get_i(stdin, prompt, "trajecs", &par_buf.trajecs);
    IF_OK status += get_f(stdin, prompt, "traj_length", &par_buf.traj_length);

    // Number of fermion and gauge steps
    IF_OK status += get_i(stdin, prompt, "nstep", &par_buf.nsteps[0]);
    IF_OK status += get_i(stdin, prompt, "nstep_gauge", &par_buf.nsteps[1]);

    // Trajectories between propagator measurements
    IF_OK status += get_i(stdin, prompt, "traj_between_meas",
                          &par_buf.propinterval);

    // lambda, mu
    IF_OK status += get_f(stdin, prompt, "lambda", &par_buf.lambda);
    IF_OK status += get_f(stdin, prompt, "mu", &par_buf.mu);

    // Maximum conjugate gradient iterations
    IF_OK status += get_i(stdin, prompt, "max_cg_iterations", &par_buf.niter);

    // Error per site for conjugate gradient
    IF_OK {
      status += get_f(stdin, prompt, "error_per_site", &x);
      par_buf.rsqmin = x * x;
    }

#ifdef EIG
    // Number of eigenvalues to calculate
    IF_OK status += get_i(stdin, prompt, "Nvec", &par_buf.Nvec);
    IF_OK status += get_f(stdin, prompt, "eig_tol", &par_buf.eig_tol);
    IF_OK status += get_i(stdin, prompt, "maxIter", &par_buf.maxIter);
#endif

#ifdef PHASE
    // Optional checkpointing for pfaffian computation
    IF_OK status += get_i(stdin, prompt, "ckpt_load", &par_buf.ckpt_load);
    IF_OK status += get_i(stdin, prompt, "ckpt_save", &par_buf.ckpt_save);
#endif

    // Find out what kind of starting lattice to use
    IF_OK status += ask_starting_lattice(stdin, prompt, &par_buf.startflag,
                                         par_buf.startfile);

    // Find out what to do with lattice at end
    IF_OK status += ask_ending_lattice(stdin, prompt, &(par_buf.saveflag),
                                       par_buf.savefile);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  traj_length = par_buf.traj_length;
  nsteps[0] = par_buf.nsteps[0];
  nsteps[1] = par_buf.nsteps[1];

  propinterval = par_buf.propinterval;
  niter = par_buf.niter;
  rsqmin = par_buf.rsqmin;

  lambda = par_buf.lambda;
  mu = par_buf.mu;

#ifdef BMN
  mass_so3 = mu * mu / 9.0;
  mass_so6 = 0.25 * mass_so3;
  mass_Myers = 1.41421356 * mu / 3.0; 
  mass_fermion = 0.25 * mu;
#else   // BFSS case: mass_so3 = mass_so6 = mu^2; other two shouldn't be used
  mass_so3 = mu * mu;
  mass_so6 = mass_so3;
  mass_Myers = -999.0;
  mass_fermion = -999.0;
#endif

  kappa = (Real)NCOL * 0.25 / lambda;
  node0_printf("lambda=%.4g --> kappa=Nc/(4lambda)=%.4g\n", lambda, kappa);

#ifdef EIG
  // Include some mallocs here (make_fields has already been called)
  // (This memory usage will not be reported to the user)
  Nvec = par_buf.Nvec;
  eigVal = malloc(sizeof *eigVal * Nvec);
  eigVec = malloc(sizeof *eigVec * Nvec);
  for (i = 0; i < Nvec; i++) {
    eigVec[i] = malloc(sizeof(matrix*) * NFERMION);
    for (j = 0; j < NFERMION; j++)
      FIELD_ALLOC(eigVec[i][j], matrix);
  }

  eig_tol = par_buf.eig_tol;
  maxIter = par_buf.maxIter;
#endif

#ifdef PHASE
  ckpt_load = par_buf.ckpt_load;
  ckpt_save = par_buf.ckpt_save;
#endif

  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile);

  // Allocate some more arrays to be used by LAPACK in scalar_eig.c
  // Needs to be above reunitarization
  work = malloc(sizeof *work * 4 * NCOL);
  store = malloc(sizeof *store * 2 * NCOL * NCOL);
  Rwork = malloc(sizeof *Rwork * (3 * NCOL - 2));
  eigs = malloc(sizeof *eigs * NCOL);
  ipiv = malloc(sizeof *ipiv * NCOL);   // For generic/reunitarize.c

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);

  // Compute initial Gamma_X
  build_Gamma_X();

  return 0;
}
// -----------------------------------------------------------------
