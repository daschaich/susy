// -----------------------------------------------------------------
// Write out all values of NBPAVRG correlators
// The label is the member of the site structure
// The data is printed for each three-dimensional spatial separation
// Write only for t=0
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void print_var3(char *label) {
#ifdef PL_CORR
  int currentnode = 0, newnode;
  int l, x, y, z, t, node0 = 0;
  complex lbuf;

  g_sync();
  t = 0;
  for (z = 0; z < nx; z++) {
    for (y = 0; y < ny; y++) {
      for (x = 0; x < nx; x++) {
        newnode = node_number(x, y, z, t);
        if (newnode != currentnode) {  // Switch to another node
          g_sync();
          currentnode = newnode;
        }
        if (this_node == 0) {
          if (currentnode == 0) {
            l = node_index(x, y, z, t);
            lbuf = lattice[l].print_var;
          }
          else
            get_field((char *)&lbuf, sizeof(complex), currentnode);

          // See defines for avg[] above
          // for the order of the output correlation values
          if ((printf("%s %d %d %d  %.6g %.6g\n",label, x, y, z,
                      (double)lbuf.real, (double)lbuf.imag) == EOF)) {
            printf("print_var: Write error\n"); 
            terminate(1);
          }
        }
        else {  // For nodes other than 0
          if (this_node == currentnode) {
            node0_printf("This node is zero\n");
            l = node_index(x, y, z, t);
            lbuf = lattice[l].print_var;
            send_field((char *)&lbuf, sizeof(complex), node0);
          }
        }
      }
    }
  }
  g_sync();
  if (this_node == 0)
    fflush(stdout);
#endif
}
// -----------------------------------------------------------------
