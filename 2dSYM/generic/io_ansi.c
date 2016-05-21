// -----------------------------------------------------------------
// Wrappers for parallel file access using ANSI standard I/O
// These are patterned after stdio fopen, fseek, fwrite, fread, fclose
// Needed in case system doesn't use ANSI standard calls
#include "generic_includes.h"
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>

FILE *g_open(const char *filename, const char *mode) {
  return fopen(filename,mode);
}

int g_seek(FILE *stream, off_t offset, int whence) {
  return fseeko(stream,offset,whence);
}

size_t g_write(const void *ptr, size_t size, size_t nmemb,FILE *stream) {
  return fwrite(ptr,size,nmemb,stream);
}

size_t g_read(void *ptr, size_t size, size_t nmemb, FILE *stream) {
  return fread(ptr,size,nmemb,stream);
}

int g_close(FILE *stream) {
  return fclose(stream);
}
// -----------------------------------------------------------------
