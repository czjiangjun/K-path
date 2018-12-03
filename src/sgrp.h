#include "pgrp_dat.h"

#ifndef  SGRP_H
#define  SGRP_H

/*** TRICLINIC  ***/
extern int  lat_C1_sgrp[NC1_sgrp];
extern int  lat_Ci_sgrp[NCi_sgrp];

extern char  *comnt_C1_sgrp[NC1_sgrp];
extern char  *comnt_Ci_sgrp[NCi_sgrp];

extern double rC1_sgrp[NC1_sgrp][1][3];
extern double rCi_sgrp[NCi_sgrp][2][3];
   
/*** MONOCLINIC  ***/
extern int  lat_C2_sgrp[NC2_sgrp];
extern int  lat_Cs_sgrp[NCs_sgrp];
extern int  lat_C2h_sgrp[NC2h_sgrp];
   
extern char  *comnt_C2_sgrp[NC2_sgrp];
extern char  *comnt_Cs_sgrp[NCs_sgrp];
extern char  *comnt_C2h_sgrp[NC2h_sgrp];
   
extern double rC2_sgrp[NC2_sgrp][2][3];
extern double rCs_sgrp[NCs_sgrp][2][3];
extern double rC2h_sgrp[NC2h_sgrp][4][3];

/* the same as before only B-centered monoclinic */
extern int  lat_B_C2_sgrp[NC2_sgrp];
extern int  lat_B_Cs_sgrp[NCs_sgrp];
extern int  lat_B_C2h_sgrp[NC2h_sgrp];
   
extern char  *comnt_B_C2_sgrp[NC2_sgrp];
extern char  *comnt_B_Cs_sgrp[NCs_sgrp];
extern char  *comnt_B_C2h_sgrp[NC2h_sgrp];
   
extern double rB_C2_sgrp[NC2_sgrp][2][3];
extern double rB_Cs_sgrp[NCs_sgrp][2][3];
extern double rB_C2h_sgrp[NC2h_sgrp][4][3];

/*** ORTHOROMBIC ***/
extern int  lat_D2_sgrp[ND2_sgrp];
extern int  lat_C2v_sgrp[NC2v_sgrp];
extern int  lat_D2h_sgrp[ND2h_sgrp];
   
extern char  *comnt_D2_sgrp[ND2_sgrp];
extern char  *comnt_C2v_sgrp[NC2v_sgrp];
extern char  *comnt_D2h_sgrp[ND2h_sgrp];
   
extern double rD2_sgrp[ND2_sgrp][4][3];
extern double rC2v_sgrp[NC2v_sgrp][4][3];
extern double rD2h_sgrp[ND2h_sgrp][8][3];
   
/*** TETRAGONAL ***/
extern int  lat_C4_sgrp[NC4_sgrp];
extern int  lat_S4_sgrp[NS4_sgrp];
extern int  lat_C4h_sgrp[NC4h_sgrp];
extern int  lat_D4_sgrp[ND4_sgrp];
extern int  lat_C4v_sgrp[NC4v_sgrp];
extern int  lat_D2d_sgrp[ND2d_sgrp];
extern int  lat_D2d2_sgrp[ND2d2_sgrp];
extern int  lat_D4h_sgrp[ND4h_sgrp];
   
extern char  *comnt_C4_sgrp[NC4_sgrp];
extern char  *comnt_S4_sgrp[NS4_sgrp];
extern char  *comnt_C4h_sgrp[NC4h_sgrp];
extern char  *comnt_D4_sgrp[ND4_sgrp];
extern char  *comnt_C4v_sgrp[NC4v_sgrp];
extern char  *comnt_D2d_sgrp[ND2d_sgrp];
extern char  *comnt_D2d2_sgrp[ND2d2_sgrp];
extern char  *comnt_D4h_sgrp[ND4h_sgrp];
   
extern double rC4_sgrp[NC4_sgrp][4][3];
extern double rS4_sgrp[NS4_sgrp][4][3];
extern double rC4h_sgrp[NC4h_sgrp][8][3];
extern double rD4_sgrp[ND4_sgrp][8][3];
extern double rC4v_sgrp[NC4v_sgrp][8][3];
extern double rD2d_sgrp[ND2d_sgrp][8][3];
extern double rD2d2_sgrp[ND2d2_sgrp][8][3];
extern double rD4h_sgrp[ND4h_sgrp][16][3];

/*** TRIGONAL ***/
extern int  lat_C3_sgrp[NC3_sgrp];
extern int  lat_C3i_sgrp[NC3i_sgrp];
extern int  lat_D3_sgrp[ND3_sgrp];
extern int  lat_D32_sgrp[ND32_sgrp];
extern int  lat_C3v_sgrp[NC3v_sgrp];
extern int  lat_C3v2_sgrp[NC3v2_sgrp];
extern int  lat_D3d_sgrp[ND3d_sgrp];
extern int  lat_D3d2_sgrp[ND3d2_sgrp];
   
extern char  *comnt_C3_sgrp[NC3_sgrp];
extern char  *comnt_C3i_sgrp[NC3i_sgrp];
extern char  *comnt_D3_sgrp[ND3_sgrp];
extern char  *comnt_D32_sgrp[ND32_sgrp];
extern char  *comnt_C3v_sgrp[NC3v_sgrp];
extern char  *comnt_C3v2_sgrp[NC3v2_sgrp];
extern char  *comnt_D3d_sgrp[ND3d_sgrp];
extern char  *comnt_D3d2_sgrp[ND3d2_sgrp];
   
extern double rC3_sgrp[NC3_sgrp][3][3];
extern double rC3i_sgrp[NC3i_sgrp][6][3];
extern double rD3_sgrp[ND3_sgrp][6][3];
extern double rD32_sgrp[ND32_sgrp][6][3];
extern double rC3v_sgrp[NC3v_sgrp][6][3];
extern double rC3v2_sgrp[NC3v2_sgrp][6][3];
extern double rD3d_sgrp[ND3d_sgrp][12][3];
extern double rD3d2_sgrp[ND3d2_sgrp][12][3];
   
/*** HEXAGONAL ***/
extern int  lat_C6_sgrp[NC6_sgrp];
extern int  lat_C3h_sgrp[NC3h_sgrp];
extern int  lat_C6h_sgrp[NC6h_sgrp];
extern int  lat_D6_sgrp[ND6_sgrp];
extern int  lat_C6v_sgrp[NC6v_sgrp];
extern int  lat_D3h_sgrp[ND3h_sgrp];
extern int  lat_D3h2_sgrp[ND3h2_sgrp];
extern int  lat_D6h_sgrp[ND6h_sgrp];
   
extern char  *comnt_C6_sgrp[NC6_sgrp];
extern char  *comnt_C3h_sgrp[NC3h_sgrp];
extern char  *comnt_C6h_sgrp[NC6h_sgrp];
extern char  *comnt_D6_sgrp[ND6_sgrp];
extern char  *comnt_C6v_sgrp[NC6v_sgrp];
extern char  *comnt_D3h_sgrp[ND3h_sgrp];
extern char  *comnt_D3h2_sgrp[ND3h2_sgrp];
extern char  *comnt_D6h_sgrp[ND6h_sgrp];
   
extern double rC6_sgrp[NC6_sgrp][6][3];
extern double rC3h_sgrp[NC3h_sgrp][6][3];
extern double rC6h_sgrp[NC6h_sgrp][12][3];
extern double rD6_sgrp[ND6_sgrp][12][3];
extern double rC6v_sgrp[NC6v_sgrp][12][3];
extern double rD3h_sgrp[ND3h_sgrp][12][3];
extern double rD3h2_sgrp[ND3h2_sgrp][12][3];
extern double rD6h_sgrp[ND6h_sgrp][24][3];

/*** CUBIC  ***/
extern int  lat_T_sgrp[NT_sgrp];
extern int  lat_Th_sgrp[NTh_sgrp];
extern int  lat_O_sgrp[NO_sgrp];
extern int  lat_Td_sgrp[NTd_sgrp];
extern int  lat_Oh_sgrp[NOh_sgrp];

extern char  *comnt_T_sgrp[NT_sgrp];
extern char  *comnt_Th_sgrp[NTh_sgrp];
extern char  *comnt_O_sgrp[NO_sgrp];
extern char  *comnt_Td_sgrp[NTd_sgrp];
extern char  *comnt_Oh_sgrp[NOh_sgrp];

extern double rT_sgrp[NT_sgrp][12][3];
extern double rTh_sgrp[NTh_sgrp][24][3];
extern double rO_sgrp[NO_sgrp][24][3];
extern double rTd_sgrp[NTd_sgrp][24][3];
extern double rOh_sgrp[NOh_sgrp][48][3];

#endif
