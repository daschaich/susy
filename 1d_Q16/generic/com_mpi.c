// -----------------------------------------------------------------
// Communications routines for MPI
// Exported functions:
// initialize_machine()   Do machine dependent setup at the very beginning
// normal_exit()          Close communications and exit
// terminate()            Halt program abruptly and exit
// machine_type()         Return string describing communications architecture
// mynode()               Return node number of this node
// numnodes()             Return number of nodes
// g_sync()               Provide a synchronization point for all nodes
// g_floatsum()           Sum a Real over all nodes
// g_vecfloatsum()        Sum a vector of Reals over all nodes
// g_doublesum()          Sum a double over all nodes
// g_vecdoublesum()       Sum a vector of doubles over all nodes
// g_complexsum()         Sum a generic precision complex number over all nodes
// g_veccomplexsum()      Sum a vector of generic precision complex numbers
//                          over all nodes
// g_dcomplexsum()        Sum a double precision complex number over all nodes
// g_vecdcomplexsum()     Sum a vector of double_complex over all nodes
// g_xor32()              Find global exclusive or of 32-bit word
// g_floatmax()           Find maximum Real over all nodes
// g_doublemax()          Find maximum double over all nodes
// broadcast_float()      Broadcast a generic precision number from
//                          node 0 to all nodes
// broadcast_double()     Broadcast a double precision number
// broadcast_complex()    Broadcast a generic precision complex number
// broadcast_dcomplex()   Broadcast a double precision complex number
// broadcast_bytes()      Broadcast a number of bytes
// send_integer()         Send an integer to one other node
// receive_integer()      Receive an integer
// send_field()           Send a field to one other node
// get_field()            Receive a field from some other node
// dclock()               Return a double precision time, with arbitrary zero
// time_stamp()           Print wall clock time with message
// make_nn_gathers()      Make all necessary lists for communications with
//                          nodes containing neighbor sites
// make_gather()          Calculate and store necessary communications lists
//                          for a given gather mapping
// declare_gather_site()  Create a message tag that defines specific
//                          details of a site gather to be used later
// declare_gather_field() Create a message tag that defines specific
//                          details of a field gather to be used later
// prepare_gather()       Allocate buffers for a previously declared gather
//                          Will automatically be called from do_gather()
//                          if not done before
// do_gather()            Execute a previously declared gather
// wait_gather()          Wait for receives to finish,
//                          ensuring that the data have actually arrived,
//                          and set pointers to received data
// cleanup_gather()       Free all the buffers that were allocated
//                          NB: The gathered data may soon disappear
// accumulate_gather()    Combine gathers into single message tag
// declare_accumulate_gather_site()
//                        Do declare_gather_site() and accumulate_gather()
//                          in single step
// declare_accumulate_gather_field()
//                        Do declare_gather_field() and accumulate_gather()
//                          in single step
// start_gather_site()    Declare/prepare/do site gather in a single step
// start_gather_field()   Declare/prepare/do field gather in a single step
// start_general_gather_site()
//                        Start asynchronous sends and receives required
//                          to gather site data at arbitrary displacement
// start_general_gather_field()
//                        Start asynchronous sends and receives required to
//                          gather field data at arbitrary displacement
// wait_general_gather()  Wait for receives to finish, insuring that the
//                          data have actually arrived,
//                          and set pointers to received data
// cleanup_general_gather()
//                        Free all the buffers that were allocated
//                          NB: the gathered data may soon disappear
#include <time.h>
#include "generic_includes.h"
#include <mpi.h>
#if PRECISION == 1
#define OUR_MPI_REAL MPI_FLOAT
#else
#define OUR_MPI_REAL MPI_DOUBLE
#endif

#define NOWHERE -1  // Not an index in array of fields

/* message types used here */
#define SEND_INTEGER_ID    1  /* send an integer to one other node */
#define SEND_FIELD_ID      2  /* id of field sent from one node to another */
#define GENERAL_GATHER_ID  3  /* id used by general_gather routines */
#define GATHER_BASE_ID     4  /* ids greater than or equal to this are used
                                 by the gather routines */

// Macro to compute the message id
#define GATHER_ID(x) (GATHER_BASE_ID+(x))
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Internal data types
/* "comlink" is the basic structure used in gathering neighboring sites.
   Each node will maintain one such structure for each direction for each
   (other) node that contains sites that are neighbors of the sites on
   this node.  For example, if the neighbors of sites on this node
   are found on two other nodes, then this node will maintain a linked
   list of two comlink structures for gathering
*/
typedef struct comlink {
  struct comlink *nextcomlink;  /* pointer to next in list, NULL if last */
  int othernode;                /* number of the node to which we connect */
  int n_subl_connected[3];
  /* Number of sites on this node that have neighbors on other node connected
     by this "comlink" of certain parity of the receiver.
     The indices 0 and 1 refer to specific parities and the
     index 2 refers to all parities */
  int *sitelist[3];
  // Address of list of indices of a certain receiver parity whose
  // neighbors are found through this comlink.  The index is the same as above.
  /* Different comlink structures may point to the same list.
     For example, the receive list for one gather may be a send list for
     the opposite gather. */
} comlink;

/* Linked list type to store id offsets for the sender.
   Needed to match the id that receiver is expecting */
typedef struct id_list_t {
  int id_offset;           /* id offset */
  struct id_list_t *next;  /* linked list */
} id_list_t;

// Structure to hold all necessary info for a gather
typedef struct gather_t {
  int *neighbor;    /* keeps track if gather neighbor is on our node or not */
  comlink *neighborlist;         /* comlink for receiving messages */
  comlink *neighborlist_send;    /* comlink for sending messages */
  id_list_t *id_list;            /* list of id offsets for sending */
  int n_recv_msgs, n_send_msgs;  /* number of messages to receive and send */
  int offset_increment; /* total number of message ids used for this gather */
} gather_t;

/* structure to keep track of details of a declared gather */
typedef struct gmem_t {
  char *mem;            /* source (destination) address for send (receive) */
  int size;             /* size of sent field */
  int stride;           /* stride of source/destination field */
  int num;              /* number of sites in sitelist */
  int *sitelist;        /* sites gathered to/from */
  struct gmem_t *next;  /* linked list */
} gmem_t;

/* Structure to keep track of outstanding sends and receives */
typedef struct {
  int msg_node;         /* node sending or receiving message */
  int id_offset;        /* id offset for this message */
  int msg_size;         /* size of message in bytes */
  char *msg_buf;        /* address of buffer malloc'd for message */
  gmem_t *gmem;         /* linked list explaining detailed usage for buffer */
  MPI_Request msg_req;  /* message handle returned by system call */
} msg_sr_t;

/* structure to store declared gathers
   this is the actual structure used internally
   it has the same name as the typedef which contains this structure which
   the user sees */
struct msg_tag {
  int *ids;          /* array of message ids used in gather */
  int nids;          /* number of message ids used in gather */
  int nrecvs;        /* number of messages to receive in gather */
  int nsends;        /* number of messages to send in gather */
  msg_sr_t *recv_msgs;  /* array of messages to receive */
  msg_sr_t *send_msgs;  /* array of messages to send */
};

// Global variables for communications stuff
/* message ids for gather encode a sequence number for the gather
   so that if several gathers are going at once, you can read
   the message corresponding to the right one. */
/* for computing message id in gather */
/* not needed anymore, but may be used for a check later */
static int id_offset;           // Label gathers by round-robin
static int num_gather_ids;      // Number of id offsets allowed
static int *id_array;           // Keep track of used ids
static gather_t *gather_array;  // Array storing gather setup info

// Number of gathers (mappings) that have been set up
static int n_gathers, gather_array_len;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Basic communications functions
void err_func(MPI_Comm *comm, int *stat, ...) {
  int len;
  char err_string[MPI_MAX_ERROR_STRING];

  printf("MPI error number: %i\n", *stat);
  MPI_Error_string(*stat, err_string, &len);
  printf("%s\n", err_string);
  terminate(*stat);
}

// Machine initialization
void initialize_machine(int *argc, char ***argv) {
  int i, flag, *tag_ub;
  MPI_Comm comm;
  MPI_Errhandler errhandler;

  flag = MPI_Init(argc, argv);
  comm = MPI_COMM_WORLD;
  if (flag)
    err_func(&comm, &flag);
  flag = MPI_Comm_create_errhandler(err_func, &errhandler);
  if (flag)
    err_func(&comm, &flag);
  flag = MPI_Comm_set_errhandler(MPI_COMM_WORLD, errhandler);
  if (flag)
    err_func(&comm, &flag);

  // Check if 32 bit int is set correctly
#ifdef SHORT_IS_32BIT
  if (sizeof(unsigned short) != 4) {
    printf("node%i: SHORT_IS_32BIT is set but sizeof(unsigned short) = %i\n",
           mynode(), sizeof(unsigned short));
    terminate(1);
  }
#else
  if (sizeof(unsigned int) != 4) {
    printf("node%i: SHORT_IS_32BIT is not set but sizeof(unsigned int) = %i\n",
           mynode(), (int)sizeof(unsigned int));
    terminate(1);
  }
#endif

  // Get the number of message types
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag);
  num_gather_ids = *tag_ub + 1 - GATHER_BASE_ID;
  if (num_gather_ids > 1024)
    num_gather_ids = 1024;

  id_offset = 0;
  id_array = malloc(sizeof *id_array * num_gather_ids);
  for (i = 0; i < num_gather_ids; ++i)
    id_array[i] = 0;

  n_gathers = 0;
  gather_array_len = 0;
  gather_array = NULL;
}

// Normal exit for multinode processes
void normal_exit(int status) {
  time_stamp("exit");
  g_sync();
  MPI_Finalize();
  exit(status);
}

// Terminate for multinode processes -- kill all nodes
void terminate(int status) {
  time_stamp("termination");
  printf("Termination: node%d, status = %d\n", this_node, status);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, 0);
  exit(status);
}

// Tell what kind of machine we are on
static char name[]="MPI (portable)";
char* machine_type() {
  return name;
}

// Return this node number
int mynode() {
  int node;
  MPI_Comm_rank(MPI_COMM_WORLD, &node);
  return node;
}

// Return number of nodes
int numnodes() {
  int nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &nodes);
  return nodes;
}

// Synchronize all nodes
void g_sync() {
  MPI_Barrier(MPI_COMM_WORLD);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Global sums
// Sum signed integer over all nodes
void g_intsum(int *ipt) {
  int work;
  MPI_Allreduce(ipt, &work, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  *ipt = work;
}

// Sum unsigned 32-bit integer type
void g_uint32sum(u_int32type *pt) {
  u_int32type work;
#ifdef SHORT_IS_32BIT
  MPI_Allreduce(pt, &work, 1, MPI_UNSIGNED_SHORT,
     MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Allreduce(pt, &work, 1, MPI_UNSIGNED,
     MPI_SUM, MPI_COMM_WORLD);
#endif
  *pt = work;
}

// Sum Real over all nodes
void g_floatsum(Real *fpt) {
  Real work;
  MPI_Allreduce(fpt, &work, 1, OUR_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
  *fpt = work;
}

// Sum a vector of Reals over all nodes
void g_vecfloatsum(Real *fpt, int length) {
  int i;
  Real *work = malloc(sizeof *work * length);
  MPI_Allreduce(fpt, work, length, OUR_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (i = 0; i<length; i++)
    fpt[i] = work[i];
  free(work);
}

// Sum double over all nodes
void g_doublesum(double *dpt) {
  double work;
  MPI_Allreduce(dpt, &work, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *dpt = work;
}

// Sum a vector of doubles over all nodes
void g_vecdoublesum(double *dpt, int length) {
  int i;
  double *work = malloc(sizeof *work * length);
  MPI_Allreduce(dpt, work, length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (i = 0; i < length; i++)
    dpt[i] = work[i];
  free(work);
}

// Sum complex over all nodes
void g_complexsum(complex *cpt) {
  complex work;
  MPI_Allreduce(cpt, &work, 2, OUR_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
  *cpt = work;
}

// Sum a vector of complex over all nodes
void g_veccomplexsum(complex *cpt, int length) {
  int i;
  complex *work = malloc(sizeof *work * length);
  MPI_Allreduce(cpt, work, 2 * length, OUR_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (i = 0; i < length; i++)
    cpt[i] = work[i];
  free(work);
}

// Sum double_complex over all nodes
void g_dcomplexsum(double_complex *cpt) {
  double_complex work;
  MPI_Allreduce(cpt, &work, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *cpt = work;
}

// Sum a vector of double_complex over all nodes
void g_vecdcomplexsum(double_complex *cpt, int length) {
  int i;
  double_complex *work = malloc(sizeof *work * length);
  MPI_Allreduce(cpt, work, 2 * length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (i = 0; i < length; i++)
    cpt[i] = work[i];
  free(work);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Global xor and maxima
// Global exclusive or acting on u_int32type, for checksums
void g_xor32(u_int32type *pt) {
  u_int32type work;
#ifdef SHORT_IS_32BIT
  MPI_Allreduce(pt, &work, 1, MPI_UNSIGNED_SHORT,
                MPI_BXOR, MPI_COMM_WORLD);
#else
  MPI_Allreduce(pt, &work, 1, MPI_UNSIGNED,
                MPI_BXOR, MPI_COMM_WORLD);
#endif
  *pt = work;
}

// Find maximum of Real over all nodes
void g_floatmax(Real *fpt) {
  Real work;
  MPI_Allreduce(fpt, &work, 1, OUR_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
  *fpt = work;
}

// Find maximum of double over all nodes
void g_doublemax(double *dpt) {
  double work;
  MPI_Allreduce(dpt, &work, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  *dpt = work;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Broadcasts
// Broadcast Real from node zero
void broadcast_float(Real *fpt) {
  MPI_Bcast(fpt, 1, OUR_MPI_REAL, 0, MPI_COMM_WORLD);
}

// Broadcast double from node zero
void broadcast_double(double *dpt) {
  MPI_Bcast(dpt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

// Broadcast generic precision complex number from node zero
void broadcast_complex(complex *cpt) {
  MPI_Bcast(cpt, 2, OUR_MPI_REAL, 0, MPI_COMM_WORLD);
}

// Broadcast double precision complex number from node zero
void broadcast_dcomplex(double_complex *cpt) {
  MPI_Bcast(cpt, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

// Broadcast bytes from node 0 to all others
void broadcast_bytes(char *buf, int size) {
  MPI_Bcast(buf, size, MPI_BYTE, 0, MPI_COMM_WORLD);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Send and receive integer
// Send an integer to one other node
// To be called only by the node doing the sending
void send_integer(int tonode, int *address) {
  MPI_Send(address, 1, MPI_INT, tonode, SEND_INTEGER_ID, MPI_COMM_WORLD);
}

// Receive an integer from another node
void receive_integer(int fromnode, int *address) {
  MPI_Status status;
  MPI_Recv(address, 1, MPI_INT, fromnode, SEND_INTEGER_ID,
           MPI_COMM_WORLD, &status);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Send and receive field
// To be called only by the node doing the sending
void send_field(char *buf, int size, int tonode) {
  MPI_Send(buf, size, MPI_BYTE, tonode, SEND_FIELD_ID, MPI_COMM_WORLD);
}

// To be called only by the node to which the field was sent
void get_field(char *buf, int size, int fromnode) {
  MPI_Status status;
  MPI_Recv(buf, size, MPI_BYTE, fromnode, SEND_FIELD_ID, MPI_COMM_WORLD,
           &status);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Timing routines
// Double precision CPU time in seconds
double dclock_cpu() {
  long fine;
  fine = clock();
  return (((double)fine) / CLOCKS_PER_SEC);
}

// Double precision wall clock time in seconds
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
double dclock() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}
#else
double dclock() {
  return dclock_cpu();
}
#endif

// Print time stamp
void time_stamp(char *msg) {
  time_t time_stamp;

  if (mynode() == 0) {
    time(&time_stamp);
    printf("%s: %s\n", msg, ctime(&time_stamp));
    fflush(stdout);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Functions used for gathers
// Find coordinates of neighbor
// Used by make_gather for nearest neighbor gathers
static void neighbor_coords_special(
  int t,                     // Coordinates of site
  int *dirpt,                       // Direction (e.g., TUP)
  int fb,                           // Forwards/backwards
  int *t2p)
                                    // Pointers to coordinates of neighbor
{
  int dir;

  dir = (fb==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
  *t2p = t;
  switch(dir) {
    case TUP   : *t2p = (t + 1) % nt;       break;
    case TDOWN : *t2p = (t + nt - 1) % nt;  break;
    default: printf("BOTCH: bad direction\n"); terminate(1);
  }
}

// Set up comlink structures needed by nearest neighbor gather routines
// make_lattice() must be called first
void make_nn_gathers() {
  int i, gather_parity;

  if (n_gathers != 0) {
    printf("error: make_nn_gathers must come before any make_gather\n");
    terminate(1);
  }

  gather_array_len = 8;
  gather_array = malloc(sizeof *gather_array * gather_array_len);
  if (gather_array == NULL) {
    printf("make_nn_gathers: node%d can't malloc gather_array\n", this_node);
    terminate(1);
  }

  if (nt&1)
    gather_parity = SCRAMBLE_PARITY;
  else
    gather_parity = SWITCH_PARITY;

  i = TUP;
  make_gather(neighbor_coords_special, &i, WANT_INVERSE,
              ALLOW_EVEN_ODD, gather_parity);
  // Already in desired order, no need to sort
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Functions used to make gathers
#define RECEIVE 0
#define SEND    1

static int parity_function(int t) {
  return t&1;
}

// Copy a linked list of comlinks, switching send and receive parity
static comlink* copy_list_switch(comlink *old_compt, int *send_subl) {
  if (old_compt == NULL)
    return NULL;

  int r_subl, s_subl;
  comlink *compt, *firstpt = malloc(sizeof *firstpt);
  compt = firstpt;

  do {
    compt->othernode = old_compt->othernode;
    for (r_subl=0; r_subl<2; r_subl++) {
      s_subl = send_subl[r_subl];
      compt->n_subl_connected[s_subl] = old_compt->n_subl_connected[r_subl];
      compt->sitelist[s_subl] = old_compt->sitelist[r_subl];
    }
    compt->n_subl_connected[2] = old_compt->n_subl_connected[2];
    compt->sitelist[2] = old_compt->sitelist[2];
    if (old_compt->nextcomlink != NULL)
      compt->nextcomlink = malloc(sizeof(comlink));
    else
      compt->nextcomlink = NULL;
    old_compt = old_compt->nextcomlink;
    compt = compt->nextcomlink;
  } while (old_compt != NULL);
  return firstpt;
}

// Sort a list of sites according to the order of the sites on the
// node with which they communicate
static void sort_site_list(
  int n,    /* number of elements in list */
  int *list,          /* pointer to list */
  void (*func)(int, int*, int, int*),
                        /* function which defines mapping */
  int *args,          /* arguments to pass to function */
  int forw_back)  /* look forwards or backwards in map */
{
  register int j, k, in1, in2, flag;
  register site *s;
  int t;
  int *key;

  if (n == 0)
    return;
  key = malloc(sizeof *key * n);
  if (key == NULL) {
    printf("sort_site_list: node%d can't malloc key\n", mynode());
    terminate(1);
  }

  // Construct sort key
  for (j = 0; j < n; j++) {
    s = &(lattice[list[j]]);
    func(s->t, args, forw_back, &t);
    key[j] = node_index(t);
  }

  // Bubble sort can be improved if it is too slow
  for (j = n - 1; j > 0; j--) {
    flag = 0;
    for (k = 0; k < j; k++) {
      in1 = key[k];
      in2 = key[k+1];
      if (in1 > in2) {
        flag = 1;
        key[k] = in2;
        key[k + 1] = in1;
        in1 = list[k];
        list[k] = list[k + 1];
        list[k + 1] = in1;
      }
    }
    if (flag == 0)
      break;
  }
  free(key);
}

// Make comlink for send or receive
static comlink* make_send_receive_list(
  void (*func)(int, int*, int, int*),
            /* function which defines sites to gather from */
  int *args,    /* list of arguments, to be passed to function */
  int want_even_odd,  /* ALLOW_EVEN_ODD or NO_EVEN_ODD */
  int forw_back,  /* FORWARDS or BACKWARDS */
  int send_recv,        /* SEND or RECEIVE list */
  int *n_msgs)          /* returns number of messages in list */
{
  int i, j, subl;
  site *s;
  int t;
  int *sbuf[2];  /* to be malloc'd */
  int *tbuf = malloc(sizeof *tbuf * numnodes());
  // Remember where comlinks are
  comlink **combuf = malloc(sizeof(comlink*) * numnodes());
  comlink *compt, **comptpt;
  comlink *firstpt;

  /* make temporary buffers of numnodes() integers to count numbers of
     neighbors in each sublattice on each node */
  for (subl = 0; subl < 2; subl++) {
    sbuf[subl] = malloc(sizeof(int) * numnodes());
    // Clear neighbor_numbers
    for (i = 0; i < numnodes(); i++)
      sbuf[subl][i] = 0;
  }
  for (i = 0; i < numnodes(); i++)
    tbuf[i] = 0;

  /* scan sites in lattice */
  FORALLSITES(i, s) {
    /* find coordinates, node, and sublattice of receiving site */
    if (send_recv == RECEIVE) {
      func(s->t, args, forw_back, &t);
      subl = parity_function(s->t);
    }
    else {  /* SEND */
      func(s->t, args, -forw_back, &t);
      subl = parity_function(t);
    }
    j = node_number(t);

    /* if site is off node, increment neighbor_counter */
    if (j != mynode()) {
      ++tbuf[j];
      if (want_even_odd == NO_EVEN_ODD)
        subl = 0;
      ++sbuf[subl][j];
    }
  }

  *n_msgs = 0;
  firstpt = NULL;
  comptpt = &firstpt;
  /* for each neighbor_counter that is nonzero, create a comlink */
  for (j = 0; j < numnodes(); j++) {
    if (j == mynode())
      continue;  /* not for local node */
    if (tbuf[j] == 0)
      continue;   /* no neighbors on this node */

    compt = malloc(sizeof *compt);
    *comptpt = compt;
    combuf[j] = compt;  /* to make it easy to find again */
    compt->nextcomlink = NULL;  /* currently terminates list */
    compt->othernode = j;
    compt->n_subl_connected[2] = tbuf[j];
    for (subl = 0; subl < 2; subl++) {
      compt->n_subl_connected[subl] = sbuf[subl][j];
    }
    compt->sitelist[2] = malloc(sizeof(int) * tbuf[j]);
    compt->sitelist[0] = compt->sitelist[2];
    for (subl = 1; subl < 2; subl++)
      compt->sitelist[subl] = (compt->sitelist[subl - 1]) + sbuf[subl - 1][j];
    /* sitelist[...] must be filled in later */
    comptpt = &(compt->nextcomlink);  /* linked list, if we
          extend it this will get address of next comlink. */
    ++(*n_msgs);
  }

  /* clear neighbor_numbers, to be used as counters now */
  for (subl = 0; subl < 2; subl++) {
    for (i = 0; i < numnodes(); i++)
      sbuf[subl][i] = 0;
  }

  /* scan sites in node again */
  FORALLSITES(i, s) {
    /* find coordinates, node, and sublattice of receiving site */
    if (send_recv == RECEIVE) {
      func(s->t, args, forw_back, &t);
      subl = parity_function(s->t);
    }
    else {  /* SEND */
      func(s->t, args, -forw_back, &t);
      subl = parity_function(t);
    }
    j = node_number(t);

    /* if neighbor is offnode, add to list in appropriate comlink */
    if (j != mynode()) {
      if (want_even_odd == NO_EVEN_ODD) subl = 0;
      combuf[j]->sitelist[subl][sbuf[subl][j]] = i;
      ++sbuf[subl][j];
    }
  }
  /* sort the lists of links according to the ordering of their
     even neighbors in the lower numbered node.  The list of sites
     on the lower numbered node is already in order. */
  for (compt = firstpt; compt != NULL; compt = compt->nextcomlink) {
    if (compt->othernode > this_node)
      continue;
    /* this is lower numbered node, so don't sort */
    if (send_recv == RECEIVE)
      i = forw_back;
    else i = -forw_back;
    for (subl = 0; subl < 2; subl++)
      sort_site_list(compt->n_subl_connected[subl],
          compt->sitelist[subl], func, args, i);
  }

  // Free temporary storage
  free(combuf);
  free(tbuf);
  for (subl = 0; subl < 2; subl++)
    free(sbuf[subl]);

  return firstpt;
}

// Determine tag offsets needed by sender
static id_list_t* make_id_list(
  comlink *recv,       /* neighborlist */
  int n_recv,          /* number of receives */
  comlink *send)       /* neighborlist_send */
{
  int i, *buf = malloc(sizeof *buf * n_recv);
  id_list_t *tol_top, *tol, **tol_next;
  MPI_Request sreq, *req = malloc(sizeof *req * n_recv);
  MPI_Status stat;

  for (i = 0; recv != NULL; ++i, recv=recv->nextcomlink) {
    buf[i] = i;
    MPI_Isend(&buf[i], 1, MPI_INT, recv->othernode, 0, MPI_COMM_WORLD,
         &req[i]);
  }
  if (i!=n_recv) {printf("error i!=n_recv\n"); terminate(1);}

  tol_next = &tol_top;
  while (send != NULL) {
    *tol_next = malloc(sizeof(id_list_t));
    tol = *tol_next;
    MPI_Irecv(&i, 1, MPI_INT, send->othernode, 0, MPI_COMM_WORLD, &sreq);
    MPI_Wait(&sreq, &stat);
    tol->id_offset = i;
    tol_next = &(tol->next);
    send = send->nextcomlink;
  }
  *tol_next = NULL;

  for (i = 0; i < n_recv; ++i)
    MPI_Wait(&req[i], &stat);

  free(req);
  free(buf);

  return tol_top;
}

// Determine max number of ids needed for gather
static int get_max_receives(int n_recv) {
  int work;
  MPI_Allreduce(&n_recv, &work, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return work;
}

// Add another gather to the list of tables
int make_gather(
  void (*func)(int, int*, int, int*),
                        /* function which defines sites to gather from */
  int *args,    /* list of arguments, to be passed to function */
  int inverse,    /* OWN_INVERSE, WANT_INVERSE, or NO_INVERSE */
  int want_even_odd,  /* ALLOW_EVEN_ODD or NO_EVEN_ODD */
  int parity_conserve)  /* {SAME,SWITCH,SCRAMBLE}_PARITY */
{
  int i, j, subl;
  site *s;
  int dir, t;
  int *send_subl;       /* sublist of sender for a given receiver */

  // We will have one or two more gathers
  if (inverse == WANT_INVERSE)
    n_gathers += 2;
  else
    n_gathers += 1;

  // If necessary, lengthen gather array to add more gathers
  if (n_gathers > gather_array_len) {
    gather_array_len = n_gathers;
    gather_array = realloc(gather_array,
                           gather_array_len * sizeof(*gather_array));
  }

  dir = n_gathers - 1;  // Index of gather we are working on
  gather_array[dir].neighbor = malloc(sizeof(int) * sites_on_node);
  if (gather_array[dir].neighbor == NULL) {
    printf("make_gather: node%d: no room for neighbor vector\n", this_node);
    terminate(1);
  }
  if (inverse == WANT_INVERSE) {
    dir = n_gathers - 2;  // Index of gather we are working on
    gather_array[dir].neighbor = malloc(sizeof(int) * sites_on_node);
    if (gather_array[dir].neighbor == NULL) {
      printf("make_gather: node%d no room for neighbor vector\n", this_node);
      terminate(1);
    }
  }

  if (want_even_odd == ALLOW_EVEN_ODD && parity_conserve != SCRAMBLE_PARITY) {
    send_subl = malloc(sizeof *send_subl * 2);
    if (send_subl == NULL) {
      printf("node%d: no room for send_subl\n", this_node);
      terminate(1);
    }
    for (subl = 0; subl < 2; subl++)
      send_subl[subl] = NOWHERE;
  }
  else
    send_subl = NULL;

  // Check to see if mapping has advertised parity and inverse properties
  // Also check to see if it returns legal values for coordinates
  FORALLSITES(i, s) {
    // Find coordinates of neighbor who sends us data
    func(s->t, args, FORWARDS, &t);

    if (t < 0 || t >= nt) {
      printf("Gather mapping does not stay in lattice\n");
      printf("It mapped %d to %d\n", s->t, t);
      terminate(1);
    }

    if (parity_conserve != SCRAMBLE_PARITY) {
      int r_subl, s_subl;
      r_subl = parity_function(s->t);
      s_subl = parity_function(t);

      if (want_even_odd == ALLOW_EVEN_ODD) {
        if (send_subl[r_subl] == NOWHERE)
          send_subl[r_subl] = s_subl;
        else if (send_subl[r_subl] != s_subl) {
          printf("Gather mixes up sublattices: %d vs %d\n",
                 send_subl[r_subl], s_subl);
          printf("on mapping %d -> %d\n", s->t, t);
          terminate(1);
        }
      }

      if (parity_conserve == SAME_PARITY && s_subl != r_subl) {
        printf("Gather mapping does not obey claimed SAME_PARITY\n");
        printf("It mapped %d with %d to %d with %d\n",
               s->t, r_subl, t, s_subl);
        terminate(1);
      }
      if (parity_conserve == SWITCH_PARITY && s_subl == r_subl) {
        printf("Gather mapping does not obey claimed SWITCH_PARITY\n");
        printf("It mapped %d with %d to %d with %d\n",
               s->t, r_subl, t, s_subl);
        terminate(1);
      }

      if (inverse == OWN_INVERSE) {
        int t2;
        func(t, args, FORWARDS, &t2);
        if (s->t != t2) {
          printf("Gather mapping is not its own inverse\n");
          printf("Its square mapped %d to %d\n", s->t, t2);
          terminate(1);
        }
      }
    }
  }

  // Receive lists: fill in pointers to sites which are on this node
  // NOWHERE if they are off-node
  FORALLSITES(i, s) {
    // Find coordinates of neighbor who sends us data
    func(s->t, args, FORWARDS, &t);
    j = node_number(t); /* node for neighbor site */
    /* if neighbor is on node, set up pointer */
    if (j == mynode())
      gather_array[dir].neighbor[i] = node_index(t);
    else
      gather_array[dir].neighbor[i] = NOWHERE;
  }

  // Make lists of sites which get data from other nodes
  gather_array[dir].neighborlist =
    make_send_receive_list(func, args, want_even_odd, FORWARDS, RECEIVE,
          &gather_array[dir].n_recv_msgs);

  /* SEND LISTS: */
  /* Now make lists of sites to which we send */
  /* Under some conditions, if mapping is its own inverse we can use
     the lists we have already made */
  if (inverse==OWN_INVERSE &&
      (want_even_odd!=ALLOW_EVEN_ODD || parity_conserve!=SCRAMBLE_PARITY)) {
    if (want_even_odd==NO_EVEN_ODD || parity_conserve==SAME_PARITY) {
      gather_array[dir].neighborlist_send = gather_array[dir].neighborlist;
      gather_array[dir].n_send_msgs = gather_array[dir].n_recv_msgs;
    }
    else {
      gather_array[dir].neighborlist_send =
  copy_list_switch(gather_array[dir].neighborlist, send_subl);
      gather_array[dir].n_send_msgs = gather_array[dir].n_recv_msgs;
    }
  }
  else {
    /* Make new linked list of comlinks for send lists */
    gather_array[dir].neighborlist_send =
      make_send_receive_list(func, args, want_even_odd, FORWARDS, SEND,
            &gather_array[dir].n_send_msgs);
  } /* End general case for send lists */

  gather_array[dir].id_list = make_id_list(gather_array[dir].neighborlist,
                                           gather_array[dir].n_recv_msgs,
                                           gather_array[dir].neighborlist_send);

  gather_array[dir].offset_increment = get_max_receives(gather_array[dir].n_recv_msgs);

  if (inverse != WANT_INVERSE) {
    free(send_subl);
    return dir;
  }

  // Now, if necessary, make inverse gather
  /* In most cases, we can use the same lists as the gather, in one
     form or another.  Of course, by the time you get to here
     you know that inverse = WANT_INVERSE */
  dir++;  /* inverse gather has direction one more than original */

  /* Always set up pointers to sites on this node */
  /* scan sites in lattice */
  FORALLSITES(i, s) {
    // Find coordinates of neighbor who sends us data
    func(s->t, args, BACKWARDS, &t);
    j = node_number(t); /* node for neighbor site */

    /* if neighbor is on node, set up pointer */
    if (j == mynode())
      gather_array[dir].neighbor[i] = node_index(t);
    else
      gather_array[dir].neighbor[i] = NOWHERE;
  }

  if (parity_conserve == SAME_PARITY || want_even_odd == NO_EVEN_ODD) {
    /* Use same comlinks as inverse gather, switching send and receive.
       Nearest neighbor gathers are an example of this case. */
    gather_array[dir].neighborlist = gather_array[dir-1].neighborlist_send;
    gather_array[dir].neighborlist_send = gather_array[dir-1].neighborlist;
    gather_array[dir].n_recv_msgs = gather_array[dir-1].n_send_msgs;
    gather_array[dir].n_send_msgs = gather_array[dir-1].n_recv_msgs;
  } else if (parity_conserve == SWITCH_PARITY) {
    /* make new comlinks, but use same lists as inverse gather, switching
       send and receive, switching even and odd. */
    gather_array[dir].neighborlist =
      copy_list_switch(gather_array[dir-1].neighborlist_send, send_subl);
    gather_array[dir].neighborlist_send =
      copy_list_switch(gather_array[dir-1].neighborlist, send_subl);
    gather_array[dir].n_recv_msgs = gather_array[dir-1].n_send_msgs;
    gather_array[dir].n_send_msgs = gather_array[dir-1].n_recv_msgs;
  } else {  /* general case.  Really only get here if ALLOW_EVEN_ODD
         and SCRAMBLE_PARITY */
    /* RECEIVE LISTS */
    gather_array[dir].neighborlist =
      make_send_receive_list(func, args, want_even_odd, BACKWARDS, RECEIVE,
            &gather_array[dir].n_recv_msgs);
    /* SEND LISTS */
    gather_array[dir].neighborlist_send =
      make_send_receive_list(func, args, want_even_odd, BACKWARDS, SEND,
            &gather_array[dir].n_send_msgs);
  } /* End making new lists for inverse gather */

  gather_array[dir].id_list =
    make_id_list(gather_array[dir].neighborlist,
      gather_array[dir].n_recv_msgs,
      gather_array[dir].neighborlist_send);

  gather_array[dir].offset_increment =
    get_max_receives(gather_array[dir].n_recv_msgs);

  free(send_subl);
  return (dir - 1);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gather routines
/* declare_strided_gather() returns a pointer to msg_tag which will
   be used as input to subsequent prepare_gather() (optional), do_gather(),
   wait_gather() and cleanup_gather() calls.

   This handles gathers from both the site structure and an array of
   fields and is not called directly by the user.  Instead they should
   call declare_gather_site() or declare_gather_field().

 prepare_gather() allocates buffers needed for the gather.  This call is
   optional since it will automatically be called from do_gather() if
   not explicitly called before.

 do_gather() starts the actual gather.  This may be repeated after a
    wait_gather() to repeat the exact same gather.

 wait_gather() waits for the gather to finish.

 cleanup_gather() frees memory allocated for the gather including the msg_tag.

   example:
  msg_tag *tag;
  tag = declare_gather_site(F_OFFSET(src), sizeof(vector), TUP,
                            EVEN, gen_pt[0]);
  prepare_gather(tag);  ** this step is optional **
  do_gather(tag);
    ** do other stuff, but don't modify tag or gen_pt[0] **
  wait_gather(tag);
    ** gen_pt[0][i] now contains the address of the src
     vector (or a copy thereof) on the neighbor of site i in the
     TUP direction for all even sites i.
     Do whatever you want with it here, but don't modify tag or
     gen_pt[0].
     Do modify the source field src **
  do_gather(tag);
    ** do other stuff **
  wait_gather(tag);
    ** gen_pt[0][i] now contains the address of the modified src.
     The restart-wait may be repeated as often as desired
  cleanup_gather(tag);
    ** subsequent calls will overwrite the gathered fields. but if you
     don't clean up, you will eventually run out of space
*/

// Return msg_tag containing details for specific gather
// Handle gathers from both the site structure and an array of fields
msg_tag* declare_strided_gather(
  void *field,          /* source buffer aligned to desired field */
  int stride,           /* bytes between fields in source buffer */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int index,    /* direction to gather from. eg TUP - index into
         neighbor tables */
  int subl,   /* subl of sites whose neighbors we gather.
         It is EVENANDODD, if all sublattices are done. */
  char **dest)   /* one of the vectors of pointers */
{
  int i;          /* scratch */
  site *s;          /* scratch pointer to site */
  msg_tag *mtag;  /* message tag structure we will return a pointer to */
  msg_sr_t *mrecv, *msend; /* arrays for send and receive lists */
  gmem_t *gmem;
  comlink *compt; /* pointer to current comlink */
  gather_t *gt;         /* pointer to current gather */
  id_list_t *idl;

  gt = &gather_array[index];

  /* set pointers in sites whose neighbors are on this node.  (If all
     neighbors are on this node, this is the only thing done.) */
  if (subl == EVENANDODD) {
    FORALLSITES(i, s) {
      if (gt->neighbor[i] != NOWHERE)
        dest[i] = (char *)field + gt->neighbor[i] * stride;
    }
  }
  else {
    FORSOMEPARITY(i, s, subl) {
      if (gt->neighbor[i] != NOWHERE)
        dest[i] = (char *)field + gt->neighbor[i] * stride;
    }
  }

  switch(subl) {
    case EVEN:        subl = 0; break;
    case ODD:         subl = 1; break;
    case EVENANDODD:  subl = 2; break;
    default:  printf("ERROR: bad sublattice\n"); terminate(subl);
  }

  // Allocate the message tag
  mtag = malloc(sizeof *mtag);
  mtag->nids = gt->offset_increment;
  mtag->ids = NULL;

  /* allocate a buffer for the msg_sr_t's.  This is dynamically allocated
     because there may be an arbitrary number of gathers in progress
     in any direction. */
  i = 0;
  for (compt = gt->neighborlist; compt != NULL;
       compt = compt->nextcomlink) {
    if (compt->n_subl_connected[subl]!=0)
      ++i;
  }
  mtag->nrecvs = i;
  if (gt->n_recv_msgs == 0)
    mrecv = NULL;
  else {
    mrecv = malloc(sizeof *mrecv * gt->n_recv_msgs);
    if (mrecv == NULL) {
      printf("declare_strided_gather: node%d can't malloc mrecv\n", this_node);
      terminate(1);
    }
  }
  mtag->recv_msgs = mrecv;

  for (i = 0, compt = gt->neighborlist_send; compt != NULL;
       compt = compt->nextcomlink) {
    if (compt->n_subl_connected[subl]!=0) ++i;
  }
  mtag->nsends = i;
  if (gt->n_send_msgs == 0)
    msend = NULL;
  else {
    msend = malloc(sizeof(msg_sr_t) * gt->n_send_msgs);
    if (msend == NULL) {
      printf("NO ROOM for msend, node%d\n", mynode());
      terminate(1);
    }
  }
  mtag->send_msgs = msend;

  /* for each node which has neighbors of my sites */
  for (i = 0, compt = gt->neighborlist; compt != NULL;
       i++, compt = compt->nextcomlink) {
    if (compt->n_subl_connected[subl] == 0) continue;
    mrecv[i].msg_node = compt->othernode;
    mrecv[i].id_offset = i;
    mrecv[i].msg_size = size*compt->n_subl_connected[subl];
    mrecv[i].msg_buf = NULL;
    gmem = malloc(sizeof(gmem_t));
    mrecv[i].gmem = gmem;
    gmem->num = compt->n_subl_connected[subl];
    gmem->sitelist = compt->sitelist[subl];
    gmem->mem = (char *)dest;
    gmem->stride = sizeof(char *);
    gmem->size = size;
    gmem->next = NULL;
  }

  /* for each node whose neighbors I have */
  idl = gt->id_list;
  for (i=0, compt = gt->neighborlist_send; compt != NULL;
       i++, compt = compt->nextcomlink, idl = idl->next) {
    if (compt->n_subl_connected[subl] == 0) continue;
    msend[i].msg_node = compt->othernode;
    msend[i].id_offset = idl->id_offset;
    msend[i].msg_size = size*compt->n_subl_connected[subl];
    msend[i].msg_buf = NULL;
    gmem = malloc(sizeof(gmem_t));
    msend[i].gmem = gmem;
    gmem->num = compt->n_subl_connected[subl];
    gmem->sitelist = compt->sitelist[subl];
    gmem->mem = field;
    gmem->stride = stride;
    gmem->size = size;
    gmem->next = NULL;
  }

  return mtag;
}

// Allocate buffers for gather
void prepare_gather(msg_tag *mtag) {
  int i, j, nids;
  int *ids;
  msg_sr_t *mrecv,*msend;
  gmem_t *gmem;
  char *tpt;

  if (mtag->ids != NULL) {
    printf("error: already prepared\n");
    terminate(1);
  }

  nids = mtag->nids;
  if (nids != 0) {
    ids = malloc(sizeof *ids * nids);
    mtag->ids = ids;
    for (i = 0, j = id_offset; i < nids; i++, j = (j + 1) % num_gather_ids) {
      // Find next available type
      while (id_array[j] != 0) {
        j = (j + 1) % num_gather_ids;
        if (j==id_offset) {
          printf("error: not enough message ids\n");
          terminate(1);
        }
      }
      ids[i] = j;
      id_array[j] = 1;
    }
    id_offset = j;
  }

  mrecv = mtag->recv_msgs;
  /* for each node which has neighbors of my sites */
  for (i = 0; i < mtag->nrecvs; ++i) {
    if (mrecv[i].msg_size == 0) {
      node0_printf("error: unexpected zero msg_size\n");
      terminate(1);
    }
    tpt = malloc(mrecv[i].msg_size);
    mrecv[i].msg_buf = tpt;
    if (tpt == NULL) {
      printf("NO ROOM for msg_buf, node%d\n", mynode());
      terminate(1);
    }
    // Set pointers in sites to correct location
    gmem = mrecv[i].gmem;
    do {
      for (j = 0; j < gmem->num; ++j,tpt += gmem->size) {
        ((char **)gmem->mem)[gmem->sitelist[j]] = tpt;
      }
    } while ((gmem=gmem->next) != NULL);
  }

  msend = mtag->send_msgs;
  // For each node whose neighbors I have
  for (i = 0; i < mtag->nsends; ++i) {
    msend[i].msg_buf = malloc(msend[i].msg_size);
    if (msend[i].msg_buf == NULL) {
      printf("NO ROOM for msg_buf, node%d\n", mynode());
      terminate(1);
    }
  }
}

// Actually execute the gather using mtag returned by start_gather_site
void do_gather(msg_tag *mtag) {
  register int i, j;
  register char *tpt; /* scratch pointer in buffers */
  msg_sr_t *mbuf;
  gmem_t *gmem;

  if ((mtag->ids == NULL) && (mtag->nids != 0))
    prepare_gather(mtag);

  mbuf = mtag->recv_msgs;
  // For each node which has neighbors of my sites
  for (i = 0; i < mtag->nrecvs; i++) {
    // Post receive
    MPI_Irecv(mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE, MPI_ANY_SOURCE,
              GATHER_ID(mtag->ids[mbuf[i].id_offset]), MPI_COMM_WORLD,
              &mbuf[i].msg_req);
  }

  mbuf = mtag->send_msgs;
  // For each node whose neighbors I have
  for (i = 0; i < mtag->nsends; ++i) {
    // Gather data into the buffer
    tpt = mbuf[i].msg_buf;
    gmem = mbuf[i].gmem;
    do {
      for (j = 0; j < gmem->num; ++j, tpt += gmem->size) {
        memcpy(tpt, gmem->mem + gmem->sitelist[j] * gmem->stride, gmem->size);
      }
    } while ((gmem = gmem->next) != NULL);
    // Start the send
    MPI_Isend(mbuf[i].msg_buf, mbuf[i].msg_size, MPI_BYTE, mbuf[i].msg_node,
              GATHER_ID(mtag->ids[mbuf[i].id_offset]), MPI_COMM_WORLD,
              &mbuf[i].msg_req);
  }
}

// Wait for gather to finish
void wait_gather(msg_tag *mtag) {
  MPI_Status status;
  int i;

  /* wait for all receive messages */
  for (i = 0; i < mtag->nrecvs; i++) {
    MPI_Wait(&mtag->recv_msgs[i].msg_req, &status);
  }

  /* wait for all send messages */
  for (i = 0; i < mtag->nsends; i++) {
    MPI_Wait(&mtag->send_msgs[i].msg_req, &status);
  }
}

// Free buffers associated with message tag
void cleanup_gather(msg_tag *mtag) {
  int i;
  gmem_t *gmem, *next;

  if (mtag->ids != NULL)
    for (i = 0; i < mtag->nids; ++i)
      id_array[mtag->ids[i]] = 0;

  // Free all receive buffers
  for (i = 0; i < mtag->nrecvs; i++) {
    free(mtag->recv_msgs[i].msg_buf);
    gmem = mtag->recv_msgs[i].gmem;
    do {
      next = gmem->next;
      free(gmem);
      gmem = next;
    } while (gmem != NULL);
  }

  // Free all send buffers
  for (i = 0; i < mtag->nsends; i++) {
    free(mtag->send_msgs[i].msg_buf);
    gmem = mtag->send_msgs[i].gmem;
    do {
      next = gmem->next;
      free(gmem);
      gmem = next;
    } while (gmem != NULL);
  }

  // Free the msg_tag buffer
  free(mtag->recv_msgs);
  free(mtag->send_msgs);
  free(mtag->ids);
  free(mtag);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Convenience routines for gathers
// Declare gather with a field offset
msg_tag* declare_gather_site(
  field_offset field, /* which field? Some member of structure "site" */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int index,    /* direction to gather from. eg TUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  return declare_strided_gather((char *)lattice + field, sizeof(site), size,
         index, parity, dest);
}

// Old style gather routine: declare and start in one call
msg_tag* start_gather_site(
  field_offset field, /* which field? Some member of structure "site" */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int index,    /* direction to gather from. eg TUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  msg_tag *mt;

  mt = declare_strided_gather((char *)lattice + field, sizeof(site), size,
                              index, parity, dest);
  prepare_gather(mt);
  do_gather(mt);

  return mt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gather routines from an array of fields
// Declare a gather from an array of fields
msg_tag* declare_gather_field(
  void *field,   /* which field? Pointer returned by malloc() */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int index,    /* direction to gather from. eg TUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  return declare_strided_gather(field, size, size, index, parity, dest);
}

// Old style gather routine: declare and start in one call
msg_tag* start_gather_field(
  void *field,   /* which field? Pointer returned by malloc() */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int index,    /* direction to gather from. eg TUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  msg_tag *mt;

  mt = declare_strided_gather(field, size, size, index, parity, dest);
  prepare_gather(mt);
  do_gather(mt);

  return mt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Multi-gather routines
/* accumulate_gather(msg_tag **mtag, msg_tag *tag)
   Joins declared gathers together under a single msg_tag.
   The second argument (tag) would be merged with the first (mtag).
   If mtag is NULL then this just copies tag into mtag.

 declare_accumulate_gather_site() declares and joins gathers.

 example:

   msg_tag *tag1, *tag2, *mtag;

   tag1 = declare_gather_site(F_OFFSET(phi), sizeof(vector), TUP,
                              EVEN, gen_pt1);
   tag2 = declare_gather_site(F_OFFSET(phi), sizeof(vector), TDOWN,
                              EVEN, gen_pt2);
   mtag = NULL;
   accumulate_gather(&mtag, tag1);
   accumulate_gather(&mtag, tag2);
   prepare_gather(mtag);  ** optional **
   do_gather(mtag);
   wait_gather(mtag);
   ** stuff **
   do_gather(tag1);     ** this is valid as long as the combined gather
   wait_gather(tag1);      (mtag) has been waited on **
   ** stuff **
   do_gather(mtag);
   wait_gather(mtag);
   cleanup_gather(mtag);
   cleanup_gather(tag1);
   cleanup_gather(tag2);

 Note that mtag must be set to NULL first in this case.
 If there is no need to use the single gathers alone one could do:

   msg_tag *mtag;

   mtag = NULL;
   declare_accumulate_gather_site(&mtag, F_OFFSET(phi), sizeof(vector),
                                  TUP, EVEN, gen_pt1);
   declare_accumulate_gather_site(&mtag, F_OFFSET(phi), sizeof(vector),
                                  TDOWN, EVEN, gen_pt2);
   prepare_gather(mtag);  ** optional **
   do_gather(mtag);
   wait_gather(mtag);
   ** stuff **
   do_gather(mtag);
   wait_gather(mtag);
   cleanup_gather(mtag);

 one could also replace
   mtag = NULL;
   declare_accumulate_gather_site(&mtag, F_OFFSET(phi), sizeof(vector),
                                  TUP, EVEN, gen_pt1);
 with
   mtag = declare_gather_site(F_OFFSET(phi), sizeof(vector), TUP,
                    EVEN, gen_pt1);
 since they do the same thing, however the first form is a bit more uniform
 in the given example.
*/

// Copy the gmem_t structure
static void copy_gmem(gmem_t **dest, gmem_t *src) {
  while (*dest != NULL)
    dest = &((*dest)->next);
  do {
    *dest = malloc(sizeof(gmem_t));
    if (*dest == NULL) {
      printf("copy_gmem: node%d can't malloc dest\n", this_node);
      terminate(1);
    }
    memcpy(*dest, src, sizeof(gmem_t));
    dest = &((*dest)->next);
    src = src->next;
  } while (src != NULL);
  *dest = NULL;
}

// Merge a source msg_sr_t structure into the dest
static void add_msgt(msg_sr_t **dest, int *ndest,
                     msg_sr_t *src, int nsrc, int nids) {

  int i, j, n = 0;

  for (i = 0; i < nsrc; ++i) {
    for (j = 0; j < *ndest; ++j) {
      if ((*dest)[j].msg_node == src[i].msg_node) {
        ++n;
        break;
      }
    }
  }
  n = *ndest + nsrc - n;

  if (n != 0) {
    *dest = realloc(*dest, n * sizeof(msg_sr_t));
    if (*dest == NULL) {
      printf("add_msgt: node%d can't realloc dest\n", this_node);
      terminate(1);
    }
    for (i = 0; i < nsrc; ++i) {
      for (j = 0; j < *ndest; ++j) {
        if ((*dest)[j].msg_node==src[i].msg_node) break;
      }
      if (j < *ndest) {
        (*dest)[j].msg_size += src[i].msg_size;
        copy_gmem(&((*dest)[j].gmem), src[i].gmem);
      }
      else {
        (*dest)[*ndest+i].msg_node = src[i].msg_node;
        (*dest)[*ndest+i].id_offset = nids + src[i].id_offset;
        (*dest)[*ndest+i].msg_size = src[i].msg_size;
        (*dest)[*ndest+i].msg_buf = NULL;
        (*dest)[*ndest+i].gmem = NULL;
        copy_gmem(&((*dest)[*ndest+i].gmem), src[i].gmem);
      }
    }
  }
  *ndest = n;
}

// Merge an already declared gather
void accumulate_gather(msg_tag **mmtag, msg_tag *mtag) {
  msg_tag *amtag;

  if (*mmtag == NULL) {
    amtag = malloc(sizeof *amtag);
    if (amtag == NULL) {
      printf("accumulate_gather: node%d can't malloc amtag\n",mynode());
      terminate(1);
    }
    amtag->nids = 0;
    amtag->ids = NULL;
    amtag->nrecvs = 0;
    amtag->recv_msgs = NULL;
    amtag->nsends = 0;
    amtag->send_msgs = NULL;
    *mmtag = amtag;
  }
  else
    amtag = *mmtag;

  add_msgt(&(amtag->recv_msgs), &(amtag->nrecvs),
      mtag->recv_msgs, mtag->nrecvs, amtag->nids);
  add_msgt(&(amtag->send_msgs), &(amtag->nsends),
      mtag->send_msgs, mtag->nsends, amtag->nids);
  amtag->nids += mtag->nids;
}

// Declare and merge gather, handling both site structure and array of fields
static void declare_accumulate_strided_gather(
  msg_tag **mmtag,      /* tag to accumulate gather into */
  void *field,          /* which field? Some member of structure "site" */
  int stride,           /* bytes between fields in source buffer */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int index,    /* direction to gather from. eg TUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  msg_tag *mtag;

  mtag = declare_strided_gather(field, stride, size, index, parity, dest);
  if (*mmtag == NULL)
    *mmtag = mtag;
  else {
    accumulate_gather(mmtag, mtag);
    cleanup_gather(mtag);
  }
}

// Declare and merge gather from field offset
void declare_accumulate_gather_site(
  msg_tag **mmtag,
  field_offset field, /* which field? Some member of structure "site" */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int index,    /* direction to gather from. eg TUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  declare_accumulate_strided_gather(mmtag, (char *)lattice + field,
             sizeof(site), size, index, parity, dest);
}

// Declare and merge gather from an array of fields
void declare_accumulate_gather_field(
  msg_tag **mmtag,
  void *field,   /* which field? Pointer returned by malloc() */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int index,    /* direction to gather from. eg TUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  declare_accumulate_strided_gather(mmtag, field, size, size, index, parity,
             dest);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// General gather routines
/* start_general_gather_site() returns a msg_tag which will
   be used as input to subsequent wait_general_gather() and
   cleanup_general_gather() calls.

   usage: tag = start_general_gather_site(source, size, displacement, parity, dest)
   example:
  msg_tag *tag;
  int disp = 1;
  tag = start_general_gather_site(F_OFFSET(phi), sizeof(vector), disp,
      EVEN, gen_pt[0]);
    ** do other stuff **
  wait_general_gather(tag);
    ** gen_pt[0][i] now contains the address of the phi
     vector (or a copy thereof) on the neighbor of site i in the
     TUP direction for all even sites i.
     Do whatever you want with it here.
    **
  cleanup_general_gather(tag);
    ** subsequent calls will overwrite the gathered fields. but if you
     don't clean up, you will eventually run out of space **
*/

struct msg_tmp { int node, count; }; /* temporary structure for keeping track
          of messages to be sent or received */
static struct msg_tmp *to_nodes, *from_nodes; /* arrays for messages */
static int g_gather_flag=0; /* flag to tell if general gather in progress */
static int tsize;     /* size of entry in messages =2*sizeof(int)+size */
static char **tdest;     /* tdest is copy of dest */
/* from_nodes, tsize and tdest are global because they are set in
   start_general_gather_site() and used in wait_general_gather().  This
   works because we allow only one general_gather in progress at a
   time. */

msg_tag* start_general_strided_gather(
  char *field,          /* source buffer aligned to desired field */
  int stride,           /* bytes between fields in source buffer */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int displacement,  /* displacement to gather from */
  int parity,   /* parity of sites to which we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  register int i, j;
  register site *s;
  register char *tpt; /* scratch pointer in buffers */
  int nsites;   /* number of sites in this receive or send */
  int disp_parity;  /* parity of displacement vector */
  int send_parity;  /* parity of sites that may be sent */
  int tt;       /* temporary coordinates */
  int othernode;    /* node sent to or received from */
  msg_sr_t *mrecv, *msend;
  msg_tag *mtag;    /* message tag, to be returned */
  int n_send_msgs, n_recv_msgs;

  /* check for gather already in progress */
  if (g_gather_flag != 0) {
    printf("ERROR: node%d, two general_gathers() at once!\n", mynode());
    terminate(1);
  }
  n_recv_msgs = 0;
  n_send_msgs = 0;
  tsize = 2 * sizeof(int) + size;   // Align pointer to double word
  tdest = dest;
  /* find parity of sites that may be sent */
  if (displacement % 2 == 0)
    disp_parity = EVEN;
  else
    disp_parity = ODD;
  switch(parity) {
    case EVEN:
      if (disp_parity == EVEN)
        send_parity = EVEN;
      else
        send_parity = ODD;
      break;
    case ODD:
      if (disp_parity == EVEN)
        send_parity = ODD;
      else
        send_parity = EVEN;
      break;
    default: // EVENANDODD
      if (parity != EVENANDODD) {
  printf("ERROR: bad parity\n");
  terminate(parity);
      }
      send_parity = EVENANDODD;
      break;
  }

  /* set pointers in sites whose neighbors are on this node.  (If all
     neighbors are on this node, this is the only thing done.) Make
     list of nodes from whom we expect messages */
  FORSOMEPARITY(i, s, parity) {
    if (displacement != 0)
      tt = (s->t + displacement + nt) % nt;
    else
      tt = s->t;
    othernode = node_number(tt);
    if (othernode == this_node)
      dest[i] = field + node_index(tt) * stride;
    else {
      for (j=0;j<n_recv_msgs;j++) if (from_nodes[j].node==othernode) break;
      if (j < n_recv_msgs)
  from_nodes[j].count++;
      else {
  if (n_recv_msgs == 0) {
    from_nodes = malloc(sizeof *from_nodes);
    from_nodes[0].node = othernode;
    from_nodes[0].count = 1;
    n_recv_msgs++;
  }
  else {
    from_nodes = realloc(from_nodes,
                         (n_recv_msgs + 1) * sizeof(*from_nodes));
    from_nodes[j].node = othernode;
    from_nodes[j].count = 1;
    n_recv_msgs++;
  }
      }
    }
  }

  /* scan sites of parity we are sending, make list of nodes to which
     we must send messages and the number of messages to each. */
  FORSOMEPARITY(i, s, send_parity) {
    if (displacement != 0)
      tt = (s->t - displacement + nt) % nt;
    else
      tt = s->t;
    othernode = node_number(tt);
    if (othernode != this_node) {
      for (j = 0; j < n_send_msgs; j++) if (to_nodes[j].node == othernode) break;
      if (j < n_send_msgs) {
  to_nodes[j].count++;
      }
      else {
  if (n_send_msgs == 0) {
    to_nodes = malloc(sizeof *to_nodes);
    to_nodes[0].node = othernode;
    to_nodes[0].count = 1;
    n_send_msgs++;
  }
  else {
    to_nodes = realloc(to_nodes, (n_send_msgs + 1) * sizeof(*to_nodes));
    to_nodes[j].node = othernode;
    to_nodes[j].count = 1;
    n_send_msgs++;
  }
      }
    }
  }

  mtag = malloc(sizeof *mtag);
  if (n_recv_msgs == 0)
    mrecv = NULL;
  else {
    mrecv = malloc(sizeof *mrecv * n_recv_msgs);
    if (mrecv == NULL) {
      printf("start_general_strided_gather: node%d can't malloc mrecv\n",
             mynode());
      terminate(1);
    }
  }
  if (n_send_msgs == 0)
    msend = NULL;
  else {
    msend = malloc(sizeof *msend * n_send_msgs);
    if (msend == NULL) {
      printf("start_general_strided_gather: node%d can't malloc msend\n",
             mynode());
      terminate(1);
    }
  }
  mtag->recv_msgs = mrecv;
  mtag->send_msgs = msend;

  mtag->nrecvs = n_recv_msgs;
  mtag->nsends = n_send_msgs;

  /* for each node which has neighbors of my sites */
  for (i = 0; i < n_recv_msgs; i++) {
    /* allocate buffer to receive neighbors */
    nsites = from_nodes[i].count;
    mrecv[i].msg_node = from_nodes[i].node;
    mrecv[i].msg_size = tsize * nsites;
    mrecv[i].msg_buf = malloc(mrecv[i].msg_size);
    if (mrecv[i].msg_buf == NULL) {
      printf("start_general_strided_gather: node%d can't malloc msg_buf\n",
             mynode());
      terminate(1);
    }
    /* post receive */
    MPI_Irecv(mrecv[i].msg_buf, nsites*tsize, MPI_BYTE,
         from_nodes[i].node, GENERAL_GATHER_ID,
         MPI_COMM_WORLD, &mrecv[i].msg_req);
  }

  /* for each node whose neighbors I have */
  for (i = 0; i < n_send_msgs; i++) {
    /* Allocate buffer to gather data. */
    tpt = malloc(tsize * to_nodes[i].count);
    if (tpt == NULL) {
      printf("start_general_strided_gather: node%d can't malloc tpt\n",
             mynode());
      terminate(1);
    }
    msend[i].msg_node = to_nodes[i].node;
    msend[i].msg_size = to_nodes[i].count*tsize;
    msend[i].msg_buf = tpt;
  }

  /* reset to_node counters */
  for (i = 0; i < n_send_msgs; i++)
    to_nodes[i].count = 0;
  /* gather data into the buffers. Each entry in the buffers consists
     of the index of the site to which the data is sent, followed by
     the actual data */
  FORSOMEPARITY(i, s, send_parity) {
    tt = (s->t - displacement + nt) % nt;
    othernode = node_number(tt);
    if (othernode != this_node) {
      for (j = 0; j < n_send_msgs; j++) {
        if (to_nodes[j].node == othernode)
          break;
      }
      tpt = msend[j].msg_buf + to_nodes[j].count*tsize;
      *(int *)tpt = node_index(tt);
      /* index of site on other node */
      memcpy(tpt + 2 * sizeof(int), field + i * stride, size);
      to_nodes[j].count++;
    }
  }

  /* start the sends */
  for (i = 0; i < n_send_msgs; i++) {
    nsites = to_nodes[i].count;
    MPI_Isend(msend[i].msg_buf, nsites*tsize, MPI_BYTE,
              to_nodes[i].node, GENERAL_GATHER_ID,
              MPI_COMM_WORLD, &msend[i].msg_req);
  }

  // Free temporary arrays
  if (n_send_msgs > 0)
    free(to_nodes);
  /* mark gather in progress and return */
  g_gather_flag = 1;

  return mtag;
}

msg_tag* start_general_gather_site(
  field_offset field, /* which field? Some member of structure "site" */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int displacement,  /* displacement to gather from */
  int parity,   /* parity of sites to which we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  return start_general_strided_gather((char *)lattice + field, sizeof(site),
                                      size, displacement, parity, dest);
}

msg_tag* start_general_gather_field(
  void *field,         /* which field? Pointer returned by malloc() */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int displacement,  /* displacement to gather from */
  int parity,   /* parity of sites to which we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  return start_general_strided_gather(field, size, size,
                                      displacement, parity, dest);
}

// Wait for a general gather to complete
void wait_general_gather(msg_tag *mtag) {
  int i, j, k;
  MPI_Status status;

  g_gather_flag = 0;
  for (i = 0; i < mtag->nrecvs; i++) {
    MPI_Wait(&mtag->recv_msgs[i].msg_req, &status);
    /* set pointers in sites to correct location */
    for (j=0; j<from_nodes[i].count; j++) {
      /* k = index of site on this node, sent in message */
      k = *(int *)(mtag->recv_msgs[i].msg_buf + j * tsize);
      tdest[k] = mtag->recv_msgs[i].msg_buf + j * tsize + 2 * sizeof(int);
    }
  }
  if (i > 0)
    free(from_nodes);
}

// Free memory associated with general gather
void cleanup_general_gather(msg_tag *mtag) {
  int i;
  MPI_Status status;

  /* free all receive buffers */
  for (i = 0; i < mtag->nrecvs; i++) {
    free(mtag->recv_msgs[i].msg_buf);
  }
  /* wait for all send messages, free all send buffers */
  for (i = 0; i < mtag->nsends; i++) {
    MPI_Wait(&mtag->send_msgs[i].msg_req, &status);
    free(mtag->send_msgs[i].msg_buf);
  }
  /* free the msg_tag buffer */
  free(mtag->recv_msgs);
  free(mtag->send_msgs);
  free(mtag);
}
// -----------------------------------------------------------------
