// -----------------------------------------------------------------
// Macros and declarations for communications routines com_*.c
// Assume that the lattice is stored as an array of site structs
#ifndef _COMDEFS_H
#define _COMDEFS_H
#include "../include/int32type.h"
#include "../include/su3.h"
#include "../include/complex.h"

// Arguments to the make_gather() routine
#define FORWARDS         1
#define BACKWARDS      (-1)   // BACKWARDS = -FORWARDS
#define OWN_INVERSE      0
#define WANT_INVERSE     1
#define NO_INVERSE       2
#define ALLOW_EVEN_ODD   0
#define NO_EVEN_ODD      1
#define SAME_PARITY      0
#define SWITCH_PARITY    1
#define SCRAMBLE_PARITY  2

// msg_tag structure used for gathers
// Actual structure defined in individual com_*.c files
typedef struct msg_tag msg_tag;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Declarations for all exported routines in the com_*.c files
void initialize_machine(int *argc, char ***argv);
void normal_exit(int status);
void terminate(int status);
char* machine_type();
int mynode();
int numnodes();
void g_sync();
void g_intsum(int *ipt);
void g_uint32sum(u_int32type *pt);
void g_floatsum(Real *fpt);
void g_vecfloatsum(Real *fpt, int nReals);
void g_doublesum(double *dpt);
void g_vecdoublesum(double *dpt, int ndoubles);
void g_complexsum(complex *cpt);
void g_veccomplexsum(complex *cpt, int ncomplex);
void g_dcomplexsum(double_complex *cpt);
void g_vecdcomplexsum(double_complex *cpt, int ncomplex);
void g_xor32(u_int32type *pt );
void g_floatmax(Real *fpt);
void g_doublemax(double *dpt);
void broadcast_float(Real *fpt);
void broadcast_double(double *dpt);
void broadcast_complex(complex *cpt);
void broadcast_dcomplex(double_complex *cpt);
void broadcast_bytes(char *buf, int size);
void send_integer(int tonode, int *address);
void send_field(char *buf, int size, int tonode);
void receive_integer(int fromnode, int *address);
void get_field(char *buf, int size, int fromnode);

double dclock_cpu();
double dclock();
void time_stamp(char *msg);

void make_nn_gathers();
void sort_eight_gathers(int index);
#define sort_eight_neighborlists sort_eight_gathers

int make_gather(
  void (*func)(int, int, int, int, int*, int, int*, int*, int*, int*),
            /* function which defines sites to gather from */
  int *args,    /* list of arguments, to be passed to function */
  int inverse,    /* OWN_INVERSE, WANT_INVERSE, or NO_INVERSE */
  int want_even_odd,  /* ALLOW_EVEN_ODD or NO_EVEN_ODD */
  int parity_conserve); /* {SAME,SWITCH,SCRAMBLE}_PARITY */

msg_tag* declare_gather_site(
  field_offset field, /* which field? Some member of structure "site" */
  int size,   /* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,    /* direction to gather from. eg XUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest);  /* one of the vectors of pointers */

void prepare_gather(msg_tag *mbuf);
void do_gather(msg_tag *mbuf);
void wait_gather(msg_tag *mbuf);
void cleanup_gather(msg_tag *mbuf);

msg_tag* start_gather_site(
  field_offset field, /* which field? Some member of structure "site" */
  int size,   /* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,    /* direction to gather from. eg XUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest);  /* one of the vectors of pointers */

msg_tag* declare_gather_field(
  void *field,   /* which field? pointer returned by malloc() */
  int size,   /* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,    /* direction to gather from. eg XUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest);  /* one of the vectors of pointers */

msg_tag* declare_strided_gather(
  void *field,          /* source buffer aligned to desired field */
  int stride,           /* bytes between fields in source buffer */
  int size,   /* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,    /* direction to gather from. eg XUP - index into
         neighbor tables */
  int subl,   /* subl of sites whose neighbors we gather.
         It is EVENANDODD, if all sublattices are done. */
  char **dest);  /* one of the vectors of pointers */

msg_tag* start_gather_field(
  void *field,   /* which field? pointer returned by malloc() */
  int size,   /* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,    /* direction to gather from. eg XUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest);  /* one of the vectors of pointers */

void accumulate_gather(
  msg_tag **mmtag,      /* msg_tag to accumulate into */
  msg_tag *mtag);       /* msg_tag to add to the gather */

void declare_accumulate_gather_site(
  msg_tag **mmtag,      /* msg_tag to accumulate into */
  field_offset field, /* which field? Some member of structure "site" */
  int size,   /* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,    /* direction to gather from. eg XUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest);  /* one of the vectors of pointers */

void declare_accumulate_gather_field(
  msg_tag **mmtag,      /* msg_tag to accumulate into */
  void *field,   /* which field? pointer returned by malloc() */
  int size,   /* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,    /* direction to gather from. eg XUP - index into
         neighbor tables */
  int parity,   /* parity of sites whose neighbors we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest);  /* one of the vectors of pointers */

msg_tag* start_general_gather_site(
  field_offset field, /* which field? Some member of structure "site" */
  int size,   /* size in bytes of the field (eg sizeof(su3_vector))*/
  int *displacement,  /* displacement to gather from. four components */
  int parity,   /* parity of sites to which we gather.
         one of EVEN, ODD or EVENANDODD. */
  char **dest);  /* one of the vectors of pointers */

msg_tag* start_general_gather_field(
 void *field,          /* which field? Pointer returned by malloc() */
 int size,    /* size in bytes of the field (eg sizeof(su3_vector))*/
 int *displacement, /* displacement to gather from. four components */
 int parity,    /* parity of sites to which we gather.
         one of EVEN, ODD or EVENANDODD. */
 char **dest);   /* one of the vectors of pointers */

void wait_general_gather(msg_tag *mbuf);
void cleanup_general_gather(msg_tag *mbuf);
#endif
// -----------------------------------------------------------------
