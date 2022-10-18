// -----------------------------------------------------------------
// Routines for gauge configuration I/O
// Works for most machines; I/O wrappers are in io_ansi.c
#include "generic_includes.h"
#include "../include/io_lat.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>

// This is dangerous -- it assumes off_t = long for this compilation
#ifndef HAVE_FSEEKO
#define fseeko fseek
#endif

#define EPS 1e-6

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

/* Checksums
   The dataset from which each checksum is computed is the full gauge
   configuration for lattice files and for propagator files, the
   propagator for a single source spin-color combination.  Data in these
   files appear as a series of 32-bit floating point numbers.  We treat
   the 32-bit values as unsigned integers v(i) where i = 0,...,N-1 ranges
   over the values in the order of appearance on the file and N is the
   total count of values in the data set.  The checksum is obtained by a
   combination of bit-rotations and exclusive or operations.  It is
   designed to be commutative and associative, unlike the BSD sum
   operation, so that the data set can be read in parallel with checksum
   contributions computed for each portion read, and then combined
   afterwards.  The sum29 checksum does a left bit rotation through i mod
   29 bits and forms an exclusive or with the accumulated checksum.  The
   sum31 checksum does the same thing, but with i mod 31 bits.

   In writing the file the bit rotation is done on the number as
   represented on the architecture and consequently as written on the
   file.  In reading and checking file integrity on an architecure with a
   relatively byte-reversed representation, byte reversal of the data
   must be done before doing the bit rotation and the resulting checksum
   must be compared with the checksum recorded on the file after
   byte-reversal.
*/
// For checksums we want a 32 bit unsigned int, for which
// we have u_int32type defined in include/int32type.h which is
// included in include/io_lat.h

#define SUCCESS  0
#define FAILURE -1
#define MAX_LINE_LENGTH 1024
#define MAX_TOKENS 512

#define TOL 0.0000001   // Tolerance for floating point checks
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void swrite_data(FILE* fp, void *src, size_t size,
                 char *myname, char *descrip) {

  if (fwrite(src, size, 1, fp) != 1) {
    printf("%s: Node %d %s write error %d\n",
           myname, this_node, descrip, errno);
    fflush(stdout);
    terminate(1);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int sread_data(FILE* fp, void *src, size_t size, char *myname, char *descrip) {
  if (fread(src, size, 1, fp) != 1) {
    printf("%s: Node %d %s read error %d\n",
           myname, this_node, descrip, errno);
    fflush(stdout);
    return 1;
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int sread_byteorder(int byterevflag, FILE* fp, void *src, size_t size,
                    char *myname, char *descrip) {

  int status = sread_data(fp,src,size,myname,descrip);
  if (byterevflag == 1)
    byterevn((int32type *)src, size / sizeof(int32type));

  return status;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Write the gauge configuration header structure in serial
void swrite_gauge_hdr(FILE *fp, gauge_header *gh) {
  char myname[] = "swrite_gauge_hdr";

  swrite_data(fp, (void *)&gh->magic_number, sizeof(gh->magic_number),
              myname, "magic_number");
  swrite_data(fp, (void *)&gh->nt, sizeof(gh->nt),
              myname, "nt");
  swrite_data(fp, (void *)gh->time_stamp, sizeof(gh->time_stamp),
              myname, "time_stamp");
  swrite_data(fp, &gh->order, sizeof(gh->order), myname, "order");

  // Header byte length
  gh->header_bytes = sizeof(gh->magic_number) + sizeof(gh->nt)
                   + sizeof(gh->time_stamp) + sizeof(gh->order);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Write a data item to the gauge info file
int write_gauge_info_item(FILE *fpout,    /* ascii file pointer */
           char *keyword,   /* keyword */
           char *fmt,       /* output format -
                must use s, d, e, f, lu, or g */
           char *src,       /* address of starting data
             floating point data must be
             of type (Real) */
           int count,       /* number of data items if > 1 */
           int stride)      /* byte stride of data if
                                           count > 1 */
{

  int i, k, n;
  char *data;
  float tt;

  // Check for valid keyword
  for (i = 0; strlen(gauge_info_keyword[i]) > 0 &&
      strcmp(gauge_info_keyword[i], keyword) != 0; i++);
  if (strlen(gauge_info_keyword[i]) == 0)
    printf("write_gauge_info_item: WARNING: keyword %s not in table\n",
           keyword);

  // Write keyword
  fprintf(fpout, "%s =", keyword);

  // Write count if more than one item
  if (count > 1)
    fprintf(fpout, "[%d]", count);

  n = count;
  if (n == 0)
    n = 1;

  // Write data
  for (k = 0, data = (char *)src; k < n; k++, data += stride) {
    fprintf(fpout," ");
    if (strstr(fmt,"s") != NULL)
      fprintf(fpout,fmt,data);
    else if (strstr(fmt,"d") != NULL)
      fprintf(fpout,fmt,*(int *)data);
    else if (strstr(fmt,"lu") != NULL)
      fprintf(fpout,fmt,*(unsigned long *)data);
    else if (strstr(fmt,"e") != NULL || strstr(fmt,"f") != NULL ||
                                        strstr(fmt,"g") != NULL) {
      tt = *(Real *)data;
      fprintf(fpout,fmt,tt);
    }
    else {
      printf("write_gauge_info_item: Unrecognized data type %s\n",fmt);
      return 1;
    }
  }
  fprintf(fpout,"\n");
  return 0;
}

/*------------------------------------------------------------------------*/

/* Write a data item to a character string */
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
  int stride)      /* byte stride of data if
          count > 1 */
{

  int i,k,n;
  size_t bytes;
  char *data;
  float tt;

  /* Check for valid keyword */
  for (i = 0; strlen(gauge_info_keyword[i]) > 0 &&
      strcmp(gauge_info_keyword[i], keyword) != 0; i++);
  if (strlen(gauge_info_keyword[i]) == 0)
    printf("write_gauge_info_item: WARNING: keyword %s not in table\n",
           keyword);

  /* Write keyword */
  bytes = 0;

  snprintf(string,nstring-bytes,"%s =",keyword);
  bytes = strlen(string);
  if (bytes >= nstring)return 1;

  /* Write count if more than one item */
  if (count > 1) {
    snprintf(string+bytes, nstring-bytes, "[%d]",count);
    bytes = strlen(string);
    if (bytes >= nstring)return 1;
  }

  n = count; if (n == 0)n = 1;

  /* Write data */
  for (k = 0, data = (char *)src; k < n; k++, data += stride) {
    snprintf(string+bytes, nstring-bytes," ");
    bytes = strlen(string);
    if (bytes >= nstring)return 1;

    if (strstr(fmt,"s") != NULL) {
      snprintf(string+bytes,nstring-bytes, fmt,data);
      bytes = strlen(string);
      if (bytes >= nstring)return 1;
    }
    else if (strstr(fmt,"d") != NULL) {
      snprintf(string+bytes,nstring-bytes,fmt,*(int *)data);
      bytes = strlen(string);
      if (bytes >= nstring)return 1;
    }
    else if (strstr(fmt,"lu") != NULL) {
      snprintf(string + bytes, nstring - bytes, fmt, *(unsigned long *)data);
      bytes = strlen(string);
      if (bytes >= nstring)return 1;
    }
    else if (strstr(fmt, "e") != NULL || strstr(fmt, "f") != NULL ||
                                         strstr(fmt, "g") != NULL) {
      tt = *(Real *)data;
      snprintf(string+bytes,nstring-bytes,fmt,tt);
      bytes = strlen(string);
      if (bytes >= nstring)return 1;
    }
    else {
      printf("write_gauge_info_item: Unrecognized data type %s\n",fmt);
      return 1;
    }
  }
  snprintf(string + bytes, nstring - bytes, "\n");
  bytes = strlen(string);
  if (bytes >= nstring)return 1;

  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Open, write, and close the ASCII info file
void write_gauge_info_file(gauge_file *gf) {
  FILE *info_fp;
  gauge_header *gh;
  char info_filename[256];
  char sums[20];

  gh = gf->header;

  /* Construct header file name from lattice file name
   by adding filename extension to lattice file name */
  strcpy(info_filename, gf->filename);
  strcat(info_filename, ASCII_GAUGE_INFO_EXT);

  // Open header file
  if ((info_fp = fopen(info_filename,"w")) == NULL) {
    printf("write_gauge_info_file: Can't open ascii info file %s\n",
           info_filename);
    return;
  }

  // Write required information
  write_gauge_info_item(info_fp, "magic_number", "%d", (char *)&gh->magic_number, 0, 0);
  write_gauge_info_item(info_fp, "time_stamp", "\"%s\"", gh->time_stamp, 0, 0);
  sprintf(sums, "%x %x", gf->check.sum29, gf->check.sum31);
  write_gauge_info_item(info_fp, "checksums", "\"%s\"", sums, 0, 0);
  write_gauge_info_item(info_fp, "nt", "%d", (char *)&nt, 0, 0);
  write_appl_gauge_info(info_fp);
  fclose(info_fp);

  printf("Wrote info file %s\n", info_filename);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up the input gauge file and gauge header structures
gauge_file* setup_input_gauge_file(char *filename) {
  char myname[] = "setup_input_gauge_file";
  gauge_file *gf = malloc(sizeof(*gf));
  gauge_header *gh = malloc(sizeof(*gh));

  // Check that file structure and header structure were set up successfully
  if (gf == NULL) {
    printf("%s: Can't malloc gf\n", myname);
    terminate(1);
  }
  if (gh == NULL) {
    printf("%s: Can't malloc gh\n", myname);
    terminate(1);
  }

  gf->filename = filename;

  // Make sure compilation gave us a 32 bit integer type
  assert(sizeof(int32type) == 4);

  gf->header = gh;
  gf->rank2rcv = NULL;
  gf->check.sum29 = 0;
  gf->check.sum31 = 0;

  return gf;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up the output gauge file and gauge header structure
gauge_file* setup_output_gauge_file() {
  int i;
  gauge_file *gf = malloc(sizeof *gf);
  gauge_header *gh = malloc(sizeof *gh);
  time_t time_stamp;

  // Check that file structure and header structure were set up successfully
  if (gf == NULL) {
    printf("setup_output_gauge_file: Can't malloc gf\n");
    terminate(1);
  }
  if (gh == NULL) {
    printf("setup_output_gauge_file: Can't malloc gh\n");
    terminate(1);
  }

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  /* Load header pointer and file name */
  gf->header = gh;

  /* Initialize */
  gf->check.sum29 = 0;
  gf->check.sum31 = 0;

  /* Load header values */
  gh->magic_number = GAUGE_VERSION_NUMBER;

  gh->nt = nt;

  // Get date and time stamp using local time on node0
  if (this_node == 0) {
    time(&time_stamp);
    strcpy(gh->time_stamp, ctime(&time_stamp));
    /* For aesthetic reasons, don't leave trailing junk bytes here to be
       written to the file */
    for (i = strlen(gh->time_stamp) + 1; i < (int)sizeof(gh->time_stamp); i++)
      gh->time_stamp[i] = '\0';

    /* Remove trailing end-of-line character */
    if (gh->time_stamp[strlen(gh->time_stamp) - 1] == '\n')
      gh->time_stamp[strlen(gh->time_stamp) - 1] = '\0';
  }

  // Broadcast to all nodes
  broadcast_bytes(gh->time_stamp, sizeof(gh->time_stamp));

  return gf;
}
/*---------------------------------------------------------------------------*/
/* Read checksum and compare.  It is assumed that the file is already
   correctly positioned.
   Should be called only by one node */
void read_checksum(gauge_file *gf, gauge_check *test_gc) {
  char myname[] = "read_checksum";
  int stat;

  // Read checksums with byte reversal
  stat = sread_byteorder(gf->byterevflag, gf->fp, &gf->check.sum29,
                         sizeof(gf->check.sum29), myname, "checksum");
  if (stat != 0)
    terminate(1);

  stat = sread_byteorder(gf->byterevflag, gf->fp, &gf->check.sum31,
                         sizeof(gf->check.sum31), myname, "checksum");
  if (stat != 0)
    terminate(1);

  if (gf->check.sum29 != test_gc->sum29 ||
      gf->check.sum31 != test_gc->sum31) {
    printf("read_checksum: Checksum violation. Computed %x %x.  Read %x %x.\n",
           test_gc->sum29, test_gc->sum31,
           gf->check.sum29, gf->check.sum31);
    terminate(1);
  }
  else {
    printf("Checksums %x %x OK\n", gf->check.sum29, gf->check.sum31);
    fflush(stdout);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Write checksum to lattice file
// Assumed that the file is correctly positioned
// Call only from one node
void write_checksum(gauge_file *gf) {
  char myname[] = "write_checksum";

  swrite_data(gf->fp, &gf->check.sum29, sizeof(gf->check.sum29),
              myname, "checksum");
  swrite_data(gf->fp, &gf->check.sum31, sizeof(gf->check.sum31),
              myname, "checksum");
  printf("Checksums %x %x\n", gf->check.sum29, gf->check.sum31);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Subroutine for reading site list from gauge configuration file
// Only node0 reads this list
void read_site_list(gauge_file *gf) {
  /* All nodes allocate space for site list table, if file is not in
     natural order */

  if (gf->header->order != NATURAL_ORDER) {
    gf->rank2rcv = malloc(sizeof(int32type) * nt);
    if (gf->rank2rcv == NULL) {
      printf("read_site_list: Can't malloc rank2rcv table\n");
      terminate(1);
    }

    // Only node0 reads the site list
    if (this_node == 0) {
      /* Reads receiving site coordinate if file is not in natural order */
      if ((int)fread(gf->rank2rcv,sizeof(int32type), nt,gf->fp) != nt) {
        printf("read_site_list: Node %d site list read error %d\n",
               this_node, errno);
        terminate(1);
      }

      if (gf->byterevflag == 1)
        byterevn(gf->rank2rcv, nt);
    }

    // Broadcast result to all nodes
    broadcast_bytes((char *)gf->rank2rcv, nt * sizeof(int32type));
  }
  else
    gf->rank2rcv = NULL;  // If no site list
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Access from node0 only
int read_gauge_hdr(gauge_file *gf) {
  FILE *fp = gf->fp;
  gauge_header *gh = gf->header;
  int32type tmp, btmp;
  int stat, byterevflag = 0;
  char myname[] = "read_gauge_hdr";

  // Read header, do byte reversal, if necessary, and check consistency
  // Read and verify magic number
  stat = sread_data(fp, &gh->magic_number, sizeof(gh->magic_number),
                    myname, "magic number");
  if (stat != 0)
    terminate(1);

  tmp = gh->magic_number;
  btmp = gh->magic_number;
  byterevn((int32type *)&btmp, 1);

  /** See if header chunk is BEGI = 1111836489 for big endian
    or the byte reverse 1229407554 for little endian **/
  if (tmp == GAUGE_VERSION_NUMBER)
    byterevflag = 0;
  else if (btmp == GAUGE_VERSION_NUMBER) {
    byterevflag = 1;
    gh->magic_number = btmp;
//    printf("Reading with byte reversal\n");
    if (sizeof(float) != sizeof(int32type)) {
      printf("read_gauge_hdr: Can't byte reverse\n");
      printf("requires size of int32type(%d) = size of float(%d)\n",
             (int)sizeof(int32type), (int)sizeof(float));
      terminate(1);
    }
  }
  else if (tmp == LIME_MAGIC_NO || btmp == LIME_MAGIC_NO) {
    // LIME format suggests a SciDAC file
    // Print error, set flag and return
    printf("%s: Looks like a SciDAC-formatted file\n", myname);
    gh->magic_number = LIME_MAGIC_NO;
    return 0;
  }
  else {
    // End of the road
    printf("read_gauge_hdr: Unrecognized magic number in gauge header\n");
    printf("Expected %x but read %x\n", GAUGE_VERSION_NUMBER, tmp);
    printf("Expected %s but read %d\n", (char *)GAUGE_VERSION_NUMBER, tmp);
    terminate(1);
    return byterevflag;
  }

  // Read and process header information
  // Get lattice dimensions
  stat = sread_byteorder(byterevflag, fp, &gh->nt, sizeof(gh->nt),
                         myname, "dimensions");
  if (stat != 0)
    terminate(1);

  // Check lattice dimensions for consistency
  if (gh->nt != nt) {
    /* So we can use this routine to discover nt
     we provide that if nt = -1 initially
       we don't die */
    if (nt != -1) {
      printf("%s: Incorrect lattice dimensions %d\n",myname, gh->nt);
      fflush(stdout);
      terminate(1);
    }
    else
      nt = gh->nt;
  }
  // Read date and time stamp
  stat = sread_data(fp, gh->time_stamp, sizeof(gh->time_stamp),
                    myname, "time stamp");
  if (stat != 0)
    terminate(1);

  // Read header byte length
  gh->header_bytes = sizeof(gh->magic_number) + sizeof(gh->nt)
                   + sizeof(gh->time_stamp) + sizeof(gh->order);

  // Read data order
  stat = sread_byteorder(byterevflag, fp, &gh->order, sizeof(gh->order),
                         myname, "order parameter");
  if (stat != 0)
    terminate(1);

  return byterevflag;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Close the file and write the info file
void w_serial_f(gauge_file *gf) {
  g_sync();
  if (this_node == 0) {
    fclose(gf->fp);
    write_gauge_info_file(gf);
  }
  // Do not free gf and gf->header so calling program can use them
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return file descriptor for opened file
gauge_file* r_serial_i(char *filename) {
  gauge_header *gh;
  gauge_file *gf;
  FILE *fp;
  int byterevflag;
  char editfilename[513];

  /* All nodes set up a gauge file and gauge header structure for reading */
  gf = setup_input_gauge_file(filename);
  gh = gf->header;

  /* Node 0 alone opens the file and reads the header */
  g_sync();
  if (this_node == 0) {
    fp = fopen(filename, "rb");
    if (fp == NULL) {
      /* If this is a partition format SciDAC file the node0 name
         has an extension ".vol0000".  So try again. */
      printf("r_serial_i: Node %d can't open file %s, error %d\n",
          this_node,filename,errno);fflush(stdout);
      strncpy(editfilename,filename,504);
      editfilename[504] = '\0';  /* Just in case of truncation */
      strcat(editfilename,".vol0000");
      printf("r_serial_i: Trying SciDAC partition volume %s\n",editfilename);
      fp = fopen(editfilename, "rb");
      if (fp == NULL) {
        printf("r_serial_i: Node %d can't open file %s, error %d\n",
               this_node,editfilename,errno);fflush(stdout);terminate(1);
      }
      printf("r_serial_i: Open succeeded\n");
    }

    gf->fp = fp;

    byterevflag = read_gauge_hdr(gf);
  }

  else gf->fp = NULL;  /* The other nodes don't know about this file */

  // Broadcast the byterevflag and header structure from node0 to all nodes
  broadcast_bytes((char *)&byterevflag, sizeof(byterevflag));
  gf->byterevflag = byterevflag;
  broadcast_bytes((char *)gh,sizeof(gauge_header));

  // No further processing here if this is a SciDAC file
  if (gh->magic_number == LIME_MAGIC_NO)
    return gf;

  // Read site list and broadcast to all nodes
  read_site_list(gf);

  return gf;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Close the file on node0
void r_serial_f(gauge_file *gf) {
  g_sync();
  if (this_node == 0)
    fclose(gf->fp);

  if (gf->rank2rcv != NULL)
    free(gf->rank2rcv);

  // Do not free gf and gf->header so calling program can use them
}
// -----------------------------------------------------------------
