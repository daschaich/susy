// -----------------------------------------------------------------
// Directions, links, and macros to give their opposites
// MPI communications assume directions from 0 to 3
#ifndef _DIRS_H
#define _DIRS_H

#define NDIMS 2       // Number of dimensions
#define NODIR -1      // Not a direction
#define XUP 0
#define TUP 1
#define TDOWN 2
#define XDOWN 3
#define OPP_DIR(dir) (3 - (dir))  // Opposite spacetime direction

#define NUMLINK 2
#define OPP_LDIR(dir) (3 - (dir)) // Opposite link direction

#endif
// -----------------------------------------------------------------
