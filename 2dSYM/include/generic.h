// -----------------------------------------------------------------
// Macros and declarations for miscellaneous generic routines
#ifndef _GENERIC_H
#define _GENERIC_H

// Other generic directory declarations are elsewhere:
//   See comdefs.h for communications
//   See io_lat.h for I/O
#include <stdio.h>
#include "../include/int32type.h"
#include "../include/complex.h"
#include "../include/susy.h"
#include "../include/macros.h"
#include "../include/random.h"
#include "../include/file_types.h"
#include "../include/io_lat.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// gaugefix.c
void gaugefix(int gauge_dir,Real relax_boost,int max_gauge_iter,
              Real gfix_tol, field_offset diffmat, field_offset sumvec);

// io_helpers.c
gauge_file *save_lattice(int flag, char *filename);
gauge_file *reload_lattice(int flag, char *filename);
int ask_starting_lattice(FILE *fp, int prompt, int *flag, char *filename);
int ask_ending_lattice(FILE *fp, int prompt, int *flag, char *filename);
int ask_gauge_fix(FILE *fp, int prompt, int *flag);
void coldlat();
void randomlat();
void funnylat();
int get_f(FILE *fp, int prompt, char *variable_name_string, Real *value);
int get_i(FILE *fp, int prompt, char *variable_name_string, int *value);
int get_vi(FILE *fp, int prompt, char *variable_name_string,
           int *value, int nvalues);
int get_vf(FILE *fp, int prompt, char *variable_name_string,
           Real *value, int nvalues);
int get_s(FILE *fp, int prompt, char *variable_name_string, char *value);
int get_prompt(FILE *fp, int *value);

// layout_hyper_prime.c
void setup_layout();
int node_number(int x, int t);
int node_index(int x, int t);
size_t num_sites(int node);
const int *get_logical_dimensions();
const int *get_logical_coordinate();

// make_lattice.c
void make_lattice();
void free_lattice();

// nersc_cksum.c
void d_linktrsum(double_complex *linktrsum);
u_int32type nersc_cksum();

// remap_stdio_from_args.c
int remap_stdio_from_args(int argc, char *argv[]);

// ranstuff.c
void initialize_prn(double_prn *prn_pt, int seed, int index);
Real myrand(double_prn *prn_pt);
#endif
// -----------------------------------------------------------------
