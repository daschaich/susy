// -----------------------------------------------------------------
// Application-dependent routine for writing gauge info file

// This file is an ASCII companion to the gauge configuration file
// and contains information about the action used to generate it.
// This information is consistently written in the pattern
//     keyword value
// or
//     keyword[n] value1 value2 ... valuen
// where n is an integer.
//
// Possible keywords are listed in io_lat.h
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Write the ASCII info file
// Call this from one of the lattice output routines in io_lat4.c
void write_appl_gauge_info(FILE *fp) {
  char sums[20];

  // The file has already been opened
  // The required magic number, time stamp, and lattice dimensions
  // have already been written
  // The rest are optional
  if (startlat_p != NULL) {
    // To retain some info about the original (or previous) configuration
    write_gauge_info_item(fp, "gauge.previous.filename","\"%s\"",
                          startlat_p->filename, 0, 0);
    write_gauge_info_item(fp, "gauge.previous.time_stamp","\"%s\"",
                          startlat_p->header->time_stamp, 0, 0);
    sprintf(sums, "%x %x", startlat_p->check.sum29, startlat_p->check.sum31);
    write_gauge_info_item(fp, "gauge.previous.checksums", "\"%s\"", sums, 0, 0);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#define INFOSTRING_MAX 2048
// Follow USQCD style for record XML
char *create_QCDML() {
  size_t bytes = 0;
  char *info = malloc(INFOSTRING_MAX);
  size_t max = INFOSTRING_MAX;
  char begin[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><usqcdInfo><version>1.0</version>";
  char begin_info[] = "<info>";
  char end_info[] = "</info>";
  char end[] = "</usqcdInfo>";
  Real myplaq = g_plaq;                             // Precision conversion
  Real nersc_linktr = linktrsum.real * one_ov_N;    // Convention and precision
  char sums[20];

  snprintf(info + bytes, max - bytes, "%s", begin);
  bytes = strlen(info);

  snprintf(info + bytes, max - bytes, "<plaq>%e</plaq>", myplaq);
  bytes = strlen(info);

  bytes = strlen(info);
  snprintf(info + bytes, max - bytes, "<linktr>%e</linktr>", nersc_linktr);

  bytes = strlen(info);
  snprintf(info + bytes, max - bytes, "%s", begin_info);

  // The rest are optional
  if (startlat_p != NULL) {
    // To retain some info about the original (or previous) configuration
    bytes = strlen(info);
    sprint_gauge_info_item(info + bytes, max - bytes, "gauge.previous.filename",
                           "%s", startlat_p->filename, 0, 0);

    bytes = strlen(info);
    sprint_gauge_info_item(info + bytes, max - bytes, "gauge.previous.time_stamp",
                           "%s", startlat_p->header->time_stamp, 0, 0);
    sprintf(sums, "%x %x", startlat_p->check.sum29, startlat_p->check.sum31);

    bytes = strlen(info);
    sprint_gauge_info_item(info + bytes, max - bytes, "gauge.previous.checksums",
                           "%s", sums, 0, 0);
  }

  bytes = strlen(info);
  snprintf(info + bytes, max - bytes, "%s", end_info);

  bytes = strlen(info);
  snprintf(info + bytes, max - bytes, "%s", end);
  return info;
}

void free_QCDML(char *info) {
  if (info != NULL)
    free(info);
}
// -----------------------------------------------------------------
