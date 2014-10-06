// -----------------------------------------------------------------
#ifndef _PRECISION_H
#define _PRECISION_H
// Generic floating point type, which defaults to single precision
#ifndef PRECISION
#define PRECISION 1
#endif

#if (PRECISION == 2)
typedef double Real;
#else
typedef float Real;
#endif

#endif
// -----------------------------------------------------------------
