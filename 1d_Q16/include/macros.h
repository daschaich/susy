// -----------------------------------------------------------------
// Macros common to all applications
#ifndef _MACROS_H
#define _MACROS_H
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Constants
#define    PI 3.14159265358979323846
#define TWOPI 6.283185307179586

// Conventions for defining checkerboard parity
#define EVEN 0x02
#define ODD 0x01
#define EVENANDODD 0x03

// ASCII string length for all file names
#define MAXFILENAME 256
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// field offset and field pointer
// Used when fields in the site are arguments to subroutines
// Usage: fo = F_OFFSET(field)
//   where "field" is the name of a field in the site struct
// Usage: address = F_PT(&site, fo)
//   where &site is the address of the site and fo is a field_offset
//   Usually, the result will have to be cast to an appropriate pointer
//   It is naturally a char*
typedef int field_offset;
#define F_OFFSET(a) \
  ((field_offset)(((char *)&(lattice[0]. a ))-((char *)&(lattice[0])) ))
#define F_PT( site , fo )  ((char *)( site ) + (fo))
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Macros for looping over directions
#define FORALLDIR(dir) for (dir = TUP; dir <= TUP; dir++)

#define FORALLUPDIR(dir) for (dir = TUP; dir <= TUP; dir++)

#define FORALLUPDIRBUT(direction, dir) \
   for (dir = TUP; dir <= TUP; dir++) if (dir != direction)

// Switches EVEN and ODD, nulls EVENANDODD
#define OPP_PAR(parity) (0x03 ^ parity)
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// printf on node zero only
#define node0_printf if(this_node==0)printf
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Macros for looping over sites on a node
// Usage:
//  int i;      // Index of the site on the node
//  site *s;    // Pointer to the current site
//  FORALLSITES(i, s) {
//    ...
//  }
//
//  int subl;   // Sublattice to loop over
//  FORSOMESUBLATTICE(i, s, subl) {
//    ...
//  }

// See loopend.h for FORSOMEPARITY definition

// Standard red-black checkerboard
#define FOREVENSITES(i,s) \
    for(i=0,s=lattice;i<even_sites_on_node;i++,s++)
#define FORODDSITES(i,s) \
    for(i=even_sites_on_node,s= &(lattice[i]);i<sites_on_node;i++,s++)
#define FORSOMEPARITY(i,s,choice) \
    for( i=((choice)==ODD ? even_sites_on_node : 0 ),  \
    s= &(lattice[i]); \
    i< ( (choice)==EVEN ? even_sites_on_node : sites_on_node); \
    i++,s++)
#define FORALLSITES(i,s) \
    for(i=0,s=lattice;i<sites_on_node;i++,s++)
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Timing switches
#ifdef TIMING
#define TIC(n) tmptime[n] = -dclock();
#define TOC(n,mytimer) tmptime[n] += dclock(); mytimer+=tmptime[n];
#else
#define TIC(n)
#define TOC(n,mytimer)
#endif

#endif
// -----------------------------------------------------------------
