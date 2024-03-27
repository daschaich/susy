// -----------------------------------------------------------------
// Directions, links, and macros to give their opposites
// MPI communications assume directions from 0 to 7
// !!! TDOWN is the same as DIR_5 -- use goffset to be safe
#ifndef _DIRS_H
#define _DIRS_H

#define NDIMS 3       // Number of dimensions
#define NODIR -1      // Not a direction
#define XUP 0
#define YUP 1
#define TUP 2
#define TDOWN 3
#define YDOWN 4
#define XDOWN 5
#define OPP_DIR(dir) (5 - (dir))  // Opposite spacetime direction

#define NUMLINK 3
#define NPLAQ 3                   // NUMLINK * (NUMLINK - 1 ) / 2
#define OPP_LDIR(dir) (5 - (dir)) // Opposite link direction

#endif
// -----------------------------------------------------------------
