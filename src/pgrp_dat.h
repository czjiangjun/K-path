#ifndef   PGRP_DAT_H
#define   PGRP_DAT_H

#define  NumPgrp   37  /* Total number of point groups */

/* define Bravais lattices */
enum {
  TRICLINIC    ,
  MONOCLINIC_P ,
  MONOCLINIC_A ,
  ORTHOROMBIC_P,
  ORTHOROMBIC_I,
  ORTHOROMBIC_C,
  ORTHOROMBIC_A,
  ORTHOROMBIC_F,
  TETRAGONAL_P ,
  TETRAGONAL_I ,
  RHOMBOHEDRAL ,
  HEXAGONAL    ,
  CUBIC_P      ,
  CUBIC_I      ,
  CUBIC_F      ,
  MONOCLINIC_B ,
  ORTHOROMBIC_B=100 , /*  lattice for input only  */
  MONOCLINIC_C ,      /*  lattices for reduction only  */
  MONOCLINIC_I 

};

/*   enum for point groups  */
enum { 
  NC1,  NCi,

  NC2,  NCs,  NC2h,

  ND2,  NC2v,  ND2h,
     
  NC4,  NS4,  NC4h,  ND4,  NC4v,  ND2d,  ND2d2,  ND4h,


  NC3,  NC3i,  ND3,  ND32,  NC3v,  NC3v2,  ND3d,  ND3d2,	

  NC6,  NC3h,  NC6h,  ND6,  NC6v,  ND3h,  ND3h2,  ND6h,

  NT,  NTh,  NO,  NTd,  NOh
};

/*  symmetry operations  */
enum { 
  ONE,
  C2001,  C2010,  C2100,
  C3p111,  C3p1_1_1,  C3p_11_1,  C3p_1_11,
  C3m111,  C3m1_1_1,  C3m_11_1,  C3m_1_11,
  C2110,  C2101,  C2011,  C21_10,  C2_101,  C201_1,
  C4p001,  C4p010,  C4p100,
  C4m001,  C4m010,  C4m100,
  INV,      
  Mxy0,  Mx0z,  M0yz,
  S3p111,  S3p1_1_1,  S3p_11_1,  S3p_1_11,
  S3m111,  S3m1_1_1,  S3m_11_1,  S3m_1_11,
  Mx_xz,  M_xyx,  Mxy_y,  Mxxz,  Mxyx,  Mxyy,
  S4p001,  S4p010,  S4p100,
  S4m001,  S4m010,  S4m100,

/*** for hexagonal coordinate system ***/
  ONEh,  C3p001h,  C3m001h,  C2001h,  C6p001h,  C6m001h,  
  C2110h,  C2100h,  C2010h,  C21_10h,  C2120h,  C2210h,
  INVh,  S3p001h,  S3m001h,  Mxy0h,  S6p001h,  S6m001h,
  Mx_xzh,  Mx2xzh,  M2xxzh,  Mxxzh,  Mx0zh,  M0yzh
};

/*  numbers of space groups for each point group  */
#define   NC1_sgrp    1
#define   NCi_sgrp    1

#define   NC2_sgrp    3
#define   NCs_sgrp    4
#define   NC2h_sgrp   6

#define   ND2_sgrp    9
#define   NC2v_sgrp   22
#define   ND2h_sgrp   28

#define   NC4_sgrp    6
#define   NS4_sgrp    2
#define   NC4h_sgrp   6
#define   ND4_sgrp    10
#define   NC4v_sgrp   12
#define   ND2d_sgrp   6
#define   ND2d2_sgrp  6
#define   ND4h_sgrp   20

#define   NC3_sgrp    4
#define   NC3i_sgrp   2
#define   ND3_sgrp    3
#define   ND32_sgrp   4
#define   NC3v_sgrp   4
#define   NC3v2_sgrp  2
#define   ND3d_sgrp   2
#define   ND3d2_sgrp  4

#define   NC6_sgrp    6
#define   NC3h_sgrp   1
#define   NC6h_sgrp   2
#define   ND6_sgrp    6
#define   NC6v_sgrp   4
#define   ND3h_sgrp   2
#define   ND3h2_sgrp  2
#define   ND6h_sgrp   4

#define   NT_sgrp     5
#define   NTh_sgrp    7
#define   NO_sgrp     8
#define   NTd_sgrp    6
#define   NOh_sgrp    10

#endif
