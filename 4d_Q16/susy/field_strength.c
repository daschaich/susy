// -----------------------------------------------------------------
// Compute SU(3) field strength tensor
// Need to define six field strength indices
// Use tempmat and tempmat2 for temporary storage

// This cartoon shows how the plaquettes are calculated
// The path begins and ends at the 'O' at the corner
// F_munu is the sum of these plaquettes minus their adjoints

//  ^         --------<--------       --------<--------
//  |dir2     |               |       |               |
//  |         |               |       |               |
//  |         |               |       |               |
//  |         |               ^       |               ^
//  ------>   |               |       |               |
//    dir     |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            -------->-------O       O------->--------
//
//            --------<-------O       O-------<--------
//            |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            |               ^       |               ^
//            |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            -------->--------       -------->--------

// Convention: try to use gen_pt[0] and mtag for links in direction dir
// These are gathered from +/- dir2

#include "susy_includes.h"

// link_src is offset for matrix link[4] in site struct
// field_dest is offset for matrix fieldstrength[6] in site struct
void make_field_strength() {
  register int i, component, dir, dir2;
  register site *s;
  int j;
  complex cc;
  matrix tmat, tmat2;
  msg_tag *mtag, *mtag2;
  
  FORALLDIR(dir) {
    for(dir2 = dir+1; dir2<NUMLINK; dir2++)
    {
      component = plaq_index[dir][dir2];
      
      // +dir +dir2 plaquette
      mtag = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                               goffset[dir2], EVENANDODD, gen_pt[0]);
      mtag2 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                               goffset[dir], EVENANDODD, gen_pt[1]);
      
      wait_gather(mtag);
      wait_gather(mtag2);
      FORALLSITES(i, s) {
        mult_nn(&s->link[dir], (matrix *)(gen_pt[1][i]), &tmat);
        mult_na(&tmat, (matrix *)(gen_pt[0][i]), &tmat2);
        mult_na(&tmat2, &s->link[dir2], &tmat);
        adjoint(&tmat, &tmat2);
        sub_matrix(&tmat, &tmat2, &FS[component][i]);
      }
      cleanup_gather(mtag2);
      
      // -dir +dir2 plaquette
      // Reuse link[dir] gather from dir2 corresponding to mtag
      FORALLSITES(i, s) {
        mult_an(&s->link[dir2], &s->link[dir], &tmat);
        mult_an((matrix *)(gen_pt[0][i]), &tmat, &(tempmat[i]));
      }
      mtag2 = start_gather_field(tempmat, sizeof(matrix),
                                 goffset[dir]+1, EVENANDODD, gen_pt[1]);
      
      wait_gather(mtag2);
      FORALLSITES(i, s) {
        mult_nn(&s->link[dir2], (matrix *)(gen_pt[1][i]), &tmat);
        adjoint(&tmat, &tmat2);
        sum_matrix(&tmat, &FS[component][i]);
        dif_matrix(&tmat2, &FS[component][i]);
      }
      cleanup_gather(mtag);
      cleanup_gather(mtag2);
      
      // -dir -dir2 plaquette
      mtag = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                               goffset[dir]+1, EVENANDODD, gen_pt[0]);
      mtag2 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                goffset[dir2]+1, EVENANDODD, gen_pt[1]);
      
      wait_gather(mtag);
      wait_gather(mtag2);
      FORALLSITES(i,s){
        mult_nn((matrix *)(gen_pt[0][i]), &s->link[dir2], &(tempmat[i]));
        mult_nn((matrix *)(gen_pt[1][i]), &s->link[dir], &(tempmat2[i]));
      }
      cleanup_gather(mtag);
      cleanup_gather(mtag2);
      
      mtag = start_gather_field(tempmat, sizeof(matrix),
                                goffset[dir2]+1, EVENANDODD, gen_pt[0]);
      mtag2 = start_gather_field(tempmat2, sizeof(matrix),
                                 goffset[dir]+1, EVENANDODD, gen_pt[1]);
      
      wait_gather(mtag);
      wait_gather(mtag2);
      FORALLSITES(i,s){
        mult_an((matrix *)(gen_pt[1][i]),(matrix *)(gen_pt[0][i]), &tmat);
        adjoint(&tmat, &tmat2);
        sum_matrix(&tmat, &FS[component][i]);
        dif_matrix(&tmat2, &FS[component][i]);
      }
      cleanup_gather(mtag);
      cleanup_gather(mtag2);
      
      // +dir -dir2 plaquette
      mtag2 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                goffset[dir], EVENANDODD, gen_pt[1]);
      
      wait_gather(mtag2);
      FORALLSITES(i, s) {
        mult_an(&s->link[dir2], &s->link[dir], &tmat);
        mult_nn(&tmat, (matrix *)(gen_pt[1][i]), &tempmat[i]);
      }
      cleanup_gather(mtag2);
      
      mtag = start_gather_field(tempmat, sizeof(matrix),
                                goffset[dir2]+1, EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i,s){
        mult_na((matrix *)(gen_pt[0][i]), &s->link[dir], &tmat);
        adjoint(&tmat, &tmat2);
        sum_matrix(&tmat, &FS[component][i]);
        dif_matrix(&tmat2, &FS[component][i]);
      }
      cleanup_gather(mtag);
      
      // Make traceless
      FORALLSITES(i, s) {
        cc = trace(&FS[component][i]);
        CMULREAL(cc, one_ov_N, cc);
        for(j = 0; j < NCOL; j++)
          CDIF(FS[component][i].e[j][j], cc);
      }
    }
  }
}
// -----------------------------------------------------------------
