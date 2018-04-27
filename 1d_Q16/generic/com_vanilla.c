// -----------------------------------------------------------------
// Communications routines for single processor machines
// Trivial version: there is only one node, every site is on node
/* Exported Functions:
   initialize_machine()   Do machine dependent setup at the very beginning
   normal_exit()          Close communications and exit
   terminate()            Halt program abruptly and exit
   machine_type()         Return string describing communications architecture
   mynode()               Return node number of this node
   numnodes()             Return number of nodes
   g_sync()               Provide a synchronization point for all nodes
   g_floatsum()           Sum a Real over all nodes
   g_vecfloatsum()        Sum a vector of Reals over all nodes
   g_doublesum()          Sum a double over all nodes
   g_vecdoublesum()       Sum a vector of doubles over all nodes
   g_complexsum()         Sum a generic precision complex number over all nodes
   g_veccomplexsum()      Sum a vector of generic precision complex numbers
                            over all nodes
   g_dcomplexsum()        Sum a double precision complex number over all nodes
   g_vecdcomplexsum()     Sum a vector of double_complex over all nodes
   g_xor32()              Find global exclusive or of 32-bit word
   g_floatmax()           Find maximum Real over all nodes
   g_doublemax()          Find maximum double over all nodes
   broadcast_float()      Broadcast a generic precision number from
                            node 0 to all nodes
   broadcast_double()     Broadcast a double precision number
   broadcast_complex()    Broadcast a generic precision complex number
   broadcast_dcomplex()   Broadcast a double precision complex number
   broadcast_bytes()      Broadcast a number of bytes
   send_integer()         Send an integer to one other node
   receive_integer()      Receive an integer
   send_field()           Send a field to one other node
   get_field()            Receive a field from some other node
   dclock()               Return a double precision time, with arbitrary zero
   time_stamp()           Print wall clock time with message
   make_nn_gathers()      Make all necessary lists for communications with
                            nodes containing neighbor sites
   make_gather()          Calculate and store necessary communications lists
                            for a given gather mapping
   declare_gather_site()  Create a message tag that defines specific details
                            of a gather to be used later
   declare_gather_field() Create a message tag that defines specific
                            details of a gather from field to be used later
   prepare_gather()       Optional call that allocates buffers for a previously
                            declared gather.  Will automatically be called from
                            do_gather() if not done before
   do_gather()            Execute a previously declared gather
   wait_gather()          Wait for receives to finish, insuring that the
                            data has actually arrived
   cleanup_gather()       Free all the buffers that were allocated, WHICH
                            MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR
   accumulate_gather()    Combine gathers into single message tag
   declare_accumulate_gather_site()   Do declare_gather_site() and
                                      accumulate_gather() in single step
   declare_accumulate_gather_field()  Do declare_gather_field() and
                                      accumulate_gather() in single step
   start_gather_site()    Older function which does declare/prepare/do_gather
                           in a single step
   start_gather_field()   Older function which does
                               declare/prepare/do_gather_field
   start_general_gather_site()  starts asynchronous sends and receives required
                             to gather fields at arbitrary displacement.
   start_general_gather_field() starts asynchronous sends and receives
                             required to gather neighbors from an
                                array of fields
   wait_general_gather()   waits for receives to finish, insuring that the
                             data has actually arrived, and sets pointers to
           received data.
   cleanup_general_gather()  frees all the buffers that were allocated, WHICH
                               MEANS THAT THE GATHERED DATA MAY SOON DISAPPEAR.
*/

#include <time.h>
#include "generic_includes.h"

#define NOWHERE -1      // Not an index in array of fields
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Internal data types
// Structure to hold all necessary info for a gather
typedef struct gather_t {
  int *neighbor;    // Keep track of where the neighbor is
} gather_t;

// Structure to keep track of outstanding sends and receives
struct msg_tag {
 int dummy;         // Don't need anything
};
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Global variables for the communications stuff
// Array storing gather setup info
static gather_t *gather_array;

// Number of gathers (mappings) that have been set up
static int n_gathers, gather_array_len;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Basic communications functions -- mostly null
// Machine initialization
void initialize_machine(int *argc, char ***argv) {
  // Check if 32-bit int is set correctly
#ifdef SHORT_IS_32BIT
  if (sizeof(unsigned short) != 4) {
    printf("node%d: SHORT_IS_32BIT is set but sizeof(unsigned short) = %d\n",
           mynode(), sizeof(unsigned short));
    terminate(1);
  }
#else
  if (sizeof(unsigned int) != 4) {
    printf("node%d: SHORT_IS_32BIT is not set but sizeof(unsigned int) = %d\n",
           mynode(), (int)sizeof(unsigned int));
    terminate(1);
  }
#endif

  n_gathers = 0;
  gather_array_len = 0;
  gather_array = NULL;
}

// Normal exit for scalar processes
void normal_exit(int status) {
  time_stamp("exit");
  exit(status);
}

// Terminate for scalar processes
void terminate(int status) {
  time_stamp("termination");
  printf("Termination: node%d, status = %d\n", this_node, status);
  fflush(stdout);
  exit(status);
}

// Tell what kind of machine we are on
static char name[]="Scalar processor";
char* machine_type() {
  return name;
}

// Return this node number
int mynode() {
  return 0;
}

// Return number of nodes
int numnodes() {
  return 1;
}

// Synchronize all nodes
void g_sync() {
}

// Sum signed integer over all nodes
void g_intsum(int *ipt) {
}

// Sum unsigned 32-bit integer type
void g_uint32sum(u_int32type *pt) {
}

// Sum Real over all nodes
void g_floatsum(Real *fpt) {
}

// Sum a vector of Reals over all nodes
void g_vecfloatsum(Real *fpt, int length) {
}

// Sum double over all nodes
void g_doublesum(double *dpt) {
}

// Sum a vector of doubles over all nodes
void g_vecdoublesum(double *dpt, int ndoubles) {
}

// Sum complex over all nodes
void g_complexsum(complex *cpt) {
}

// Sum a vector of complex over all nodes
void g_veccomplexsum(complex *cpt, int ncomplex) {
}

// Sum double_complex over all nodes
void g_dcomplexsum(double_complex *cpt) {
}

// Sum a vector of double_complex over all nodes
void g_vecdcomplexsum(double_complex *cpt, int ncomplex) {
}

// Global exclusive or acting on u_int32type
void g_xor32(u_int32type *pt) {
}

// Find maximum of Real over all nodes
void g_floatmax(Real *fpt) {
}

// Find maximum of double over all nodes
void g_doublemax(double *dpt) {
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Broadcasts
// Broadcast Real from node zero
void broadcast_float(Real *fpt) {
}

// Broadcast double from node zero
void broadcast_double(double *dpt) {
}

// Broadcast generic-precision complex number from node zero
void broadcast_complex(complex *cpt) {
}

// Broadcast double-precision complex number from node zero
void broadcast_dcomplex(double_complex *cpt) {
}

// Broadcast bytes from node 0 to all others
void broadcast_bytes(char *buf, int size) {
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Send and receive integer
// Send an integer to one other node
// To be called only by the node doing the sending
void send_integer(int tonode, int *address) {
  printf("ERROR: called send_integer() in com_vanilla.c\n");
  terminate(1);
}

// Receive an integer from another node
void receive_integer(int fromnode, int *address) {
  printf("ERROR: called receive_integer() in com_vanilla.c\n");
  terminate(1);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Send and receive field
// To be called only by the node doing the sending
void send_field(char *buf, int size, int tonode) {
  printf("ERROR: called send_field() in com_vanilla.c\n");
  terminate(1);
}

// To be called only by the node to which the field was sent
void get_field(char *buf, int size, int fromnode) {
  printf("ERROR: called get_field() in com_vanilla.c\n");
  terminate(1);
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

  time(&time_stamp);
  printf("%s: %s\n", msg, ctime(&time_stamp));
  fflush(stdout);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Functions used for gathers
// Find coordinates of neighbor
// Used by make_gather for nearest neighbor gathers
static void neighbor_coords_special(
  int t,                            // Coordinates of site
  int *dirpt,                       // Direction (eg TUP)
  int fb,                           // Forwards/backwards
  int *t2p)                         // Pointers to coordinates of neighbor
{
  int dir;

  dir = (fb==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
  *t2p = t;
  switch(dir) {
    case TUP   : *t2p = (t + 1) % nt;      break;
    case TDOWN : *t2p = (t + nt - 1) % nt; break;
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

  gather_array_len = 2;
  gather_array = malloc(sizeof *gather_array * gather_array_len);
  if (gather_array == NULL) {
    printf("make_nn_gathers: node%d can't malloc gather_array\n", this_node);
    terminate(1);
  }

  if (nt&1)
    gather_parity = SCRAMBLE_PARITY;
  else
    gather_parity = SWITCH_PARITY;

  i = TUP; // Hack to reuse higher dimensional functions
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

// Add another gather to the list of tables
int make_gather(
  void (*func)(int, int*, int, int*),
                        /* function which defines sites to gather from */
  int *args,    /* list of arguments, to be passed to function */
  int inverse,    /* OWN_INVERSE, WANT_INVERSE, or NO_INVERSE */
  int want_even_odd,  /* ALLOW_EVEN_ODD or NO_EVEN_ODD */
  int parity_conserve)  /* {SAME,SWITCH,SCRAMBLE}_PARITY */
{
  int i, subl;
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
                           sizeof(*gather_array) * gather_array_len);
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

  // Receive lists: fill in pointers to sites
  FORALLSITES(i, s) {
    // Find coordinates of neighbor who sends us data
    func(s->t, args, FORWARDS, &t);
    gather_array[dir].neighbor[i] = node_index(t);
  }

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
    /* set up pointer */
    gather_array[dir].neighbor[i] = node_index(t);
  }

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
  tag = declare_gather_site(F_OFFSET(phi), sizeof(vector), TUP,
                        EVEN, gen_pt[0]);
        prepare_gather(tag);  ** this step is optional **
        do_gather(tag);
    ** do other stuff, but don't modify tag or gen_pt[0] **
  wait_gather(tag);
    ** gen_pt[0][i] now contains the address of the phi
     vector (or a copy thereof) on the neighbor of site i in the
     TUP direction for all even sites i.
     Do whatever you want with it here, but don't modify tag or
     gen_pt[0].
     Do modify the source field phi. **
  do_gather(tag);
    ** do other stuff **
  wait_gather(tag);
    ** gen_pt[0][i] now contains the address of the modified phi.
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
  gather_t *gt;         /* pointer to current gather */

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
  return NULL;
}

// Allocate buffers for gather
void prepare_gather(msg_tag *mtag) {
}

// Actually execute the gather using mtag returned by start_gather_site
void do_gather(msg_tag *mtag) {
}

// Wait for gather to finish
void wait_gather(msg_tag *mtag) {
}

// Free buffers associated with message tag
void cleanup_gather(msg_tag *mtag) {
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

// Merge already declared gather
void accumulate_gather(msg_tag **mmtag, msg_tag *mtag) {
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
  declare_accumulate_strided_gather(mmtag, (char *)field, size, size, index,
                                    parity, dest);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// General gather routines
// start_general_gather_site returns a msg_tag to be used as input
// to subsequent wait_general_gather and cleanup_general_gather calls

// Usage: tag = start_general_gather_site(src, size, disp, parity, dest)
// Example:
//  msg_tag *tag;
//  int disp = 1;    // Displacement
//  tag = start_general_gather_site(F_OFFSET(phi), sizeof(vector), disp,
//                                  EVEN, gen_pt[0]);
//
//  // Can do other stuff that doesn't depend on phi
//  wait_general_gather(tag);
//  // gen_pt[0][i] now contains the address of the phi vector
//  // (or a copy thereof) on the neighbor of site i in the TUP direction
//  // for all even sites i
//  // Do whatever you want with it here.
//  cleanup_general_gather(tag);
//  // Subsequent calls will overwrite the gathered fields
//  // If you don't clean up, you will eventually run out of space
static int g_gather_flag = 0; /* flag to tell if general gather in progress */

msg_tag* start_general_strided_gather(
  char *field,          /* source buffer aligned to desired field */
  int stride,           /* bytes between fields in source buffer */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int displacement,  /* displacement to gather from. four components */
  int subl,   /* subl of sites whose neighbors we gather.
         It is EVENANDODD, if all sublattices are done. */
  char **dest)   /* one of the vectors of pointers */
{
  site *s;      /* scratch pointer to site */
  int i;          /* scratch */
  int tt;  /* temporary coordinates */

  // Check for gather already in progress
  if (g_gather_flag != 0) {
    fprintf(stderr, "ERROR: node%d, two general_gathers() at once!\n",
            mynode());
    exit(1);
  }

  /* set pointers in sites whose neighbors are on this node.  (If all
     neighbors are on this node, this is the only thing done.) Make
     list of nodes from whom we expect messages */
  if (subl == EVENANDODD) {
    FORALLSITES(i, s) {
      if (displacement != 0)
        tt = (s->t + displacement + nt) % nt;
      else
        tt = s->t;
      dest[i] = field + stride * node_index(tt);
    }
  }
  else {
    FORSOMEPARITY(i, s, subl) {
      if (displacement != 0)
        tt = (s->t + displacement + nt) % nt;
      else
        tt = s->t;
      dest[i] = field + stride * node_index(tt);
    }
  }

  // Mark gather in progress and return
  g_gather_flag = 1;
  return NULL;
}

msg_tag* start_general_gather_site(
  field_offset field, /* which field? Some member of structure "site" */
  int size,   /* size in bytes of the field (eg sizeof(vector))*/
  int displacement,  /* displacement to gather from. four components */
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
  int displacement,  /* displacement to gather from. four components */
  int parity,   /* parity of sites to which we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest)   /* one of the vectors of pointers */
{
  return start_general_strided_gather((char *)field, size, size,
               displacement, parity, dest);
}

// Wait for a general gather to complete
void wait_general_gather(msg_tag *mtag) {
  g_gather_flag = 0;
}

// Free memory associated with general gather
void cleanup_general_gather(msg_tag *mtag) {
}
// -----------------------------------------------------------------
