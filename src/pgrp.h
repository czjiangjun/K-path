#ifndef PGRP_H
#define PGRP_H

extern int  chk_cubic(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3]);
extern int  chk_hexagonal(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3]);
extern int  chk_tetragonal(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3]);
extern int  chk_orthorombic(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3]);
extern int  chk_monoclinic(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3]);
extern int  chk_triclinic(int lat, int nop, int *ind, double Ttrcl[3][3],int *lat_new, double Tnew[3][3]);
extern int  new_transl(int lat, int nop, int *ind, double Ttrcl[3][3], int *lat_new, double Tnew[3][3]);
extern int  det_pgrp(int nop, double sym_op[][4][3]);

#endif
