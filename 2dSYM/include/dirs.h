// -----------------------------------------------------------------
// Directions, and a macro to give the opposite direction
// These must go from 0 to 3 because they will be used to index an array
// Also define NDIRS = number of directions
#ifndef _DIRS_H
#define _DIRS_H

#define XUP 0
#define TUP 1
#define TDOWN 2
#define XDOWN 3

#define NODIR -1      // Not a direction

#define OPP_DIR(dir) (3 - (dir))  // Opposite direction
#define NDIRS 4       // Number of directions
#define NDIMS 2       // Number of dimensions

#endif
// -----------------------------------------------------------------
