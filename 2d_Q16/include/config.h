// -----------------------------------------------------------------
// Collect macros for preprocessor tweaks
// Accommodate differences in compilers, architecture and OS
// NOT generated automatically by configure
#ifndef _CONFIG_H
#define _CONFIG_H
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compiler/processor-dependent macros
/* Specify the unsigned 32 bit integer base type for this compiler */
/* Run the script "getint.sh" to find out what to use */
/* One and only one of these should be defined */
#define INT_IS_32BIT 1  /* Most present systems */
#undef SHORT_IS_32BIT   /* Needed on T3E UNICOS, for example */
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compiler/OS-dependent macros
// Define if ieeefp.h exists
// Most systems don't have this
//#define HAVE_IEEEFP_H 1

// Define if unistd.h exists
// Most systems have this (exception: NT)
#define HAVE_UNISTD_H 1

// Define if <sys/time.h> exists
// Most systems have this
#define HAVE_SYS_TIME_H 1

// Define if ANSI "fseeko" present
// Most systems have this (exceptions: T3E UNICOS)
#define HAVE_FSEEKO 1

#endif
// -----------------------------------------------------------------
