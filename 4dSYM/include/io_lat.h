// -----------------------------------------------------------------
// Macros, structures and prototypes for gauge configuration input and output
#ifndef _IO_LAT_H
#define _IO_LAT_H

// Definitions of restore and save lattice flags used in io_helpers.c
#define CONTINUE      10
#define FRESH         11
#define RANDOM        12
#define RELOAD_SERIAL 13
#define FORGET        40
#define SAVE_SERIAL   42

#ifdef HAVE_UNISTD_H
#include <unistd.h>     // For write, close and off_t
#endif
#include <sys/types.h>  // For off_t
#include <stdio.h>
#include "../include/int32type.h"
#include "../include/macros.h"    // For MAXFILENAME
#include "../include/file_types.h"
#include "../include/su3.h"
#include "../include/dirs.h"      // For NDIMS
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Binary lattice format
#define MAX_TIME_STAMP 64

// 1. Header
typedef struct {
  int32type magic_number;               // Identifies file format
  char time_stamp[MAX_TIME_STAMP]; /* Date and time stamp - used to
          check consistency between the
          ASCII header file and the lattice file */
  int32type dims[NDIMS];                // Full lattice dimensions
  int32type header_bytes;               /* NOT WRITTEN TO THE FILE but
           helpful for finding the data */
  int32type order;                      /* 0 means no coordinate list is
                attached and the values are in
                coordinate serial order.
                Non-zero means that a
                coordinate list is attached,
                specifying the order of values */
} gauge_header;

// 2. Checksum
typedef struct {
  u_int32type sum31;
  u_int32type sum29;
} gauge_check;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Info file format
// List of admissible keywords for info file
// More can be added as desired
#ifdef CONTROL
char *gauge_info_keyword[] = {
      "magic_number",
      "time_stamp",
      "checksums",
      "nx",
      "ny",
      "nz",
      "nt",
      "gauge.previous.filename",
      "gauge.previous.time_stamp",
      "gauge.previous.checksums",
      ""       // Last entry MUST be a zero-length keyword
};
#else
extern char *gauge_info_keyword[];
#endif

// Used to create info file name
#define ASCII_GAUGE_INFO_EXT ".info"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Struct for gauge file
typedef struct {
  FILE *fp;             // File pointer
  gauge_header *header; // Pointer to header for file
  char *filename;       // Pointer to file name string
  int byterevflag;      // Byte reverse flag, used only for reading
  int32type *rank2rcv;  // File site list, used only for serial reading
  gauge_check check;    // Checksum
} gauge_file;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for generic/io_lat4.c
void read_lat_dim_gf(char *filename, int *ndim, int dims[]);
gauge_file *restore_serial(char *filename);
gauge_file *save_serial(char *filename);
int write_gauge_info_item( FILE *fpout, /* ascii file pointer */
           char *keyword,   /* keyword */
           char *fmt,       /* output format -
                must use s, d, f, or e */
           char *src,       /* address of starting data */
           int count,       /* number of data items if > 1 */
           int stride);     /* byte stride of data if
             count > 1 */
int sprint_gauge_info_item(
  char *string,    /* character string */
  size_t nstring,     /* string length */
  char *keyword,   /* keyword */
  char *fmt,       /* output format -
          must use s, d, e, f, or g */
  char *src,       /* address of starting data
          floating point data must be
          of type (Real) */
  int count,       /* number of data items if > 1 */
  int stride);     /* byte stride of data if
          count > 1 */
gauge_file *setup_output_gauge_file();
gauge_file *setup_input_gauge_file(char *filename);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// In <application>/gauge_info.c
void write_appl_gauge_info(FILE *fp);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for I/O routine interface in generic/io_ansi.c
FILE *g_open(const char *filename, const char *mode);
int g_seek(FILE *stream, off_t offset, int whence);
size_t g_write(const void *ptr, size_t size, size_t nmemb,FILE *stream);
size_t g_read(void *ptr, size_t size, size_t nmemb, FILE *stream);
int g_close(FILE *stream);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for generic/io_lat_util.c routines
void swrite_data(FILE* fp, void *src, size_t size, char *myname, char *descrip);
int sread_data(FILE* fp, void *src, size_t size, char *myname, char *descrip);
int sread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, char *myname, char *descrip);
void swrite_gauge_hdr(FILE *fp, gauge_header *gh);
int write_gauge_info_item( FILE *fpout,    /* ascii file pointer */
           char *keyword,   /* keyword */
           char *fmt,       /* output format -
                must use s, d, e, f, or g */
           char *src,       /* address of starting data
             floating point data must be
             of type (Real) */
           int count,       /* number of data items if > 1 */
           int stride);      /* byte stride of data if
                                           count > 1 */
int sprint_gauge_info_item(
  char *string,    /* character string */
  size_t nstring,     /* string length */
  char *keyword,   /* keyword */
  char *fmt,       /* output format -
          must use s, d, e, f, or g */
  char *src,       /* address of starting data
          floating point data must be
          of type (Real) */
  int count,       /* number of data items if > 1 */
  int stride);      /* byte stride of data if
          count > 1 */
void write_gauge_info_file(gauge_file *gf);
gauge_file *setup_input_gauge_file(char *filename);
gauge_file *setup_output_gauge_file();
void read_checksum(gauge_file *gf, gauge_check *test_gc);
void write_checksum(gauge_file *gf);
void read_site_list(gauge_file *gf);
int read_gauge_hdr(gauge_file *gf);
gauge_file *r_serial_i(char *filename);
void w_serial_f(gauge_file *gf);
void r_serial_f(gauge_file *gf);
void byterevn(int32type w[], int n);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// IO for columns and diagonal elements of matrix Q
// in pfaffian phase calculation
#ifdef PHASE
void loadQ(complex **Q, int ckpt_load);
void saveQ(complex **Q, int ckpt_save);
void load_diag(complex *diag, int ckpt_load);
void save_diag(complex *diag, int ckpt_save);
#endif
#endif
// -----------------------------------------------------------------
