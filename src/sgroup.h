#ifndef SGROUP_H
#define SGROUP_H

extern void transform_coor(t_cell *cl_in, double T[3][3]);
extern int get_Telpr(int lat, double Telpr[3][3]);
extern int get_centering_translations(int lat, int *nvec, double vec[][3]);
extern void transform_group_opr(int nop, double op[][4][3], double A[3][3]);
extern void rotate_group_transl(int nop, double op[][4][3], double A[3][3]);
extern void find_atom_pgroup(
  int lat, int nop_s, int npgrp_s, double sym_op[][4][3], double Rsh[3],
  int *lat_new, int *npgrp_p, double op[][4][3], double T[3][3]
);
		      
#endif
