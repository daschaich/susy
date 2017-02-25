// -----------------------------------------------------------------
// Define a type to ensure 32-bit integers for binary file formats
#ifndef _TYPE32_H
#define _TYPE32_H

#include "../include/config.h"

// One and only one of SHORT_IS_32BIT and INT_IS_32BIT should be defined
#if defined(SHORT_IS_32BIT) && defined(INT_IS_32BIT)
  #error "Can't define both SHORT_IS_32BIT and INT_IS_32BIT in config.h!"
#endif

#if !defined(SHORT_IS_32BIT) && !defined(INT_IS_32BIT)
  #error "Must define one of SHORT_IS_32BIT or INT_IS_32BIT in config.h!"
#endif

#ifdef SHORT_IS_32BIT
typedef short int32type;
typedef unsigned short u_int32type;
#else
typedef int int32type;
typedef unsigned int u_int32type;
#endif

#endif
// -----------------------------------------------------------------
