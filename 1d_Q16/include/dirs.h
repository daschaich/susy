// -----------------------------------------------------------------
// Directions, links, and macros to give their opposites
#ifndef _DIRS_H
#define _DIRS_H

#define NODIR -1      // Not a direction
#define TUP 0
#define TDOWN 1
#define OPP_DIR(dir) (1 - (dir))  // Opposite spacetime direction

#define NSCALAR 9
#define NFERMION 16
#define NCHIRAL_FERMION 8
#define OPP_LDIR(dir) (1 - (dir)) // Opposite link direction

#endif
// -----------------------------------------------------------------
