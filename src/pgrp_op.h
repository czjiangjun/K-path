#ifndef   PGRP_OP_H
#define   PGRP_OP_H  1

#define  NumPgrp   37

extern int nop_pgrp[NumPgrp];
extern char name_pgrp[NumPgrp][32];
extern int ind_pgrp[NumPgrp][48];
     
extern double C1_pgrp[1][3][3];
extern double Ci_pgrp[2][3][3];
extern double C2_pgrp[2][3][3];
extern double Cs_pgrp[2][3][3];
extern double C2h_pgrp[4][3][3];
extern double D2_pgrp[4][3][3];
extern double C2v_pgrp[4][3][3];
extern double D2h_pgrp[8][3][3];

extern double C4_pgrp[4][3][3];
extern double S4_pgrp[4][3][3];
extern double C4h_pgrp[8][3][3];
extern double D4_pgrp[8][3][3];
extern double C4v_pgrp[8][3][3];
extern double D2d_pgrp[8][3][3];
extern double D2d2_pgrp[8][3][3];
extern double D4h_pgrp[16][3][3];

extern double C3_pgrp[3][3][3];
extern double C3i_pgrp[6][3][3];
extern double D3_pgrp[6][3][3];
extern double D32_pgrp[6][3][3];
extern double C3v_pgrp[6][3][3];
extern double C3v2_pgrp[6][3][3];
extern double D3d_pgrp[12][3][3];
extern double D3d2_pgrp[12][3][3];
   
extern double C6_pgrp[6][3][3];
extern double C3h_pgrp[6][3][3];
extern double C6h_pgrp[12][3][3];
extern double D6_pgrp[12][3][3];
extern double C6v_pgrp[12][3][3];
extern double D3h_pgrp[12][3][3];
extern double D3h2_pgrp[12][3][3];
extern double D6h_pgrp[24][3][3];

extern double T_pgrp[12][3][3];
extern double Th_pgrp[24][3][3];
extern double O_pgrp[24][3][3];
extern double Td_pgrp[24][3][3];
extern double Oh_pgrp[48][3][3];
     
#endif
