// -----------------------------------------------------------------
// Supersymmetric setup
#include "susy_includes.h"

#define IF_OK if(status==0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size and seed, and send to others
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    // stringification kludge from GNU preprocessor manual
    // http://gcc.gnu.org/onlinedocs/cpp/Stringification.html
#define XSTR(s) STR(s)
#define STR(s) #s
    // end kludge
    printf("N=4 SYM, Nc = %d, DIMF = %d, fermion rep = " XSTR(FREP) "\n",
           NCOL, DIMF);
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
    status = get_prompt(stdin,  &prompt);

    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt, "nz", &par_buf.nz);
    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);
    IF_OK status += get_i(stdin, prompt, "iseed", &par_buf.iseed);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  nx = par_buf.nx;
  ny = par_buf.ny;
  nz = par_buf.nz;
  nt = par_buf.nt;
  iseed = par_buf.iseed;

  this_node = mynode();
  number_of_nodes = numnodes();
  volume = nx * ny * nz * nt;
  total_iters = 0;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate space for fields
void make_fields() {
  double size = (double)(2 + 4 * NUMLINK + 3 * NUMLINK * NUMLINK)
                * sites_on_node * sizeof(su3_vector);

  FIELD_ALLOC_VEC(tsite, su3_vector, NUMLINK);

  FIELD_ALLOC(site_src, su3_vector);
  FIELD_ALLOC(site_dest, su3_vector);
  FIELD_ALLOC_VEC(link_src, su3_vector, NUMLINK);
  FIELD_ALLOC_VEC(link_dest, su3_vector, NUMLINK);
  FIELD_ALLOC_VEC(link_dest2, su3_vector, NUMLINK);
  FIELD_ALLOC_MAT(plaq_src, su3_vector, NUMLINK, NUMLINK);
  FIELD_ALLOC_MAT(plaq_dest, su3_vector, NUMLINK, NUMLINK);
  FIELD_ALLOC_MAT(plaq_dest2, su3_vector, NUMLINK, NUMLINK);

  // Stout smearing stuff needed for `hot-start' random configurations
  size += (double)(1 + 3 * NUMLINK) * sites_on_node * sizeof(su3_matrix_f);
  size += (double)(NUMLINK) * sites_on_node * sizeof(anti_hermitmat);
  FIELD_ALLOC_VEC(thin_link, su3_matrix_f, NUMLINK);
  FIELD_ALLOC_VEC(smeared_link, su3_matrix_f, NUMLINK);
  FIELD_ALLOC_VEC(stp, su3_matrix_f, NUMLINK);    // Staples
  FIELD_ALLOC_VEC(Q, anti_hermitmat, NUMLINK);    // To be exponentiated
  FIELD_ALLOC(tempmat, su3_matrix_f);             // Staple storage

  node0_printf("Mallocing %.1f MBytes per core for fields\n", size / 1e6);

#ifdef PHASE
  FIELD_ALLOC(src, Twist_Fermion);
  FIELD_ALLOC(res, Twist_Fermion);

  int Ndat = 16 * DIMF;
  double tr = (double)volume * Ndat * sites_on_node * Ndat;

  // Total number of matvecs is (volume * Ndat)^2 / 4
  Nmatvecs = volume * Ndat * volume * 4 * DIMF;
  node0_printf("Q has %d columns --> %d matvecs and %.1f MBytes per core...",
               volume * Ndat, Nmatvecs, tr * sizeof(complex) / 1e6);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume and seed
  prompt = initial_set();
  // Initialize the node random number generator
  initialize_prn(&node_prn, iseed, volume + mynode());
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();
  // Set up offset tables for the five paths
  setup_offset();
  setup_qclosed_offset();
  // Allocate space for fields
  make_fields();

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for Monte Carlo
int readin(int prompt) {
  // prompt=1 indicates prompts are to be given for input
  int status;
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

    // lambda, kappa_u1, bmass and fmass
    IF_OK status += get_f(stdin, prompt, "lambda", &par_buf.lambda);
    IF_OK status += get_f(stdin, prompt, "kappa_u1", &par_buf.kappa_u1);
    IF_OK status += get_f(stdin, prompt, "bmass", &par_buf.bmass);
    IF_OK status += get_f(stdin, prompt, "fmass", &par_buf.fmass);

#ifdef STOUT
    // Stout smearing stuff
    IF_OK status += get_i(stdin, prompt, "Nstout", &par_buf.Nstout);
    IF_OK status += get_f(stdin, prompt, "rho", &par_buf.rho);
#endif

    // Maximum conjugate gradient iterations
    IF_OK status += get_i(stdin, prompt, "max_cg_iterations", &par_buf.niter);

    // Error per site for conjugate gradient
    IF_OK {
      status += get_f(stdin, prompt, "error_per_site", &x);
      par_buf.rsqmin = x;
    }

    // Number of stochastic sources for fermion bilinear and susy trans
    IF_OK status += get_i(stdin, prompt, "nsrc", &par_buf.nsrc);

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
    IF_OK status += ask_starting_lattice(stdin,  prompt, &par_buf.startflag,
                                         par_buf.startfile);

    // Find out whether or not to gauge fix to Coulomb gauge
    IF_OK status += ask_gauge_fix(stdin, prompt, &par_buf.fixflag);

    // Find out what to do with lattice at end
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
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
  fixflag = par_buf.fixflag;
  niter = par_buf.niter;
  rsqmin = par_buf.rsqmin;

  lambda = par_buf.lambda;
  kappa_u1 = par_buf.kappa_u1;
  bmass = par_buf.bmass;
  fmass = par_buf.fmass;
  kappa = (Real)NCOL * 0.5 / lambda;
  node0_printf("lambda=%.4g --> kappa=Nc/(2lambda)=%.4g\n",
               lambda, kappa);
  node0_printf("C2=%.4g\n", C2);

  nsrc = par_buf.nsrc;
#ifdef EIG
  Nvec = par_buf.Nvec;
  eig_tol = par_buf.eig_tol;
  maxIter = par_buf.maxIter;
#endif
#ifdef STOUT
  Nstout = par_buf.Nstout;
  rho = par_buf.rho;
#endif
#ifdef PHASE
  ckpt_load = par_buf.ckpt_load;
  ckpt_save = par_buf.ckpt_save;
#endif

  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile);
  strcpy(stringLFN, par_buf.savefile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);
  // Generate the irrep lattice
  fermion_rep();
  return 0;
}
// -----------------------------------------------------------------
