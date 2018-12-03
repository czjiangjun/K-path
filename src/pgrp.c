#include <stdio.h>
#include <stdlib.h>
#include "type_sg.h"
#include "pgrp_dat.h"
#include "pgrp_op.h"
#include "math_sg.h"

/* For point group determine Bravais lattice and redefine ELEMENTARY
 * translation  vectors.
 * int chk_cubic()      - check a point group to belong to cubic family
 * int chk_tetragonal() - check a point group to belong to tetragonal family
 * ..........
 *  column number (i) of Tnew is decomposition of (i) new vector */

int  chk_cubic(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3])
{
   int i,j=0;
   if( nop!= nop_pgrp[NT] && nop!= nop_pgrp[NTh] &&
       nop!= nop_pgrp[NO] && nop!= nop_pgrp[NTd] &&
       nop!= nop_pgrp[NOh]
     ) return(0);
   /*  find 3 axis 3  */
   for(i=0; i<nop; i++)
     if( ind[i]==C3p111 || ind[i]==C3p_11_1 || ind[i]==C3p_1_11 ||
         ind[i]==S3p111 || ind[i]==S3p_11_1 || ind[i]==S3p_1_11
        ) j++;
   if( j==3 || j==6 ) {
      *lat_new=lat;
       Tnew[0][0]=1.; Tnew[0][1]=0.; Tnew[0][2]=0.;
       Tnew[1][0]=0.; Tnew[1][1]=1.; Tnew[1][2]=0.;
       Tnew[2][0]=0.; Tnew[2][1]=0.; Tnew[2][2]=1.;
       return(1);
   }
   return(0);
}

int  chk_hexagonal(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3])
{
   int i;
   if( nop!= nop_pgrp[NC6] && nop!= nop_pgrp[NC6h] &&
       nop!= nop_pgrp[ND6h]
     ) return(0);

   if( lat != HEXAGONAL ) return(0);
   for(i=0; i<nop; i++)
     if( ind[i]==C6p001h || ind[i]==S6p001h ) {
      *lat_new=lat;
       Tnew[0][0]=1.; Tnew[0][1]=0.; Tnew[0][2]=0.;
       Tnew[1][0]=0.; Tnew[1][1]=1.; Tnew[1][2]=0.;
       Tnew[2][0]=0.; Tnew[2][1]=0.; Tnew[2][2]=1.;
       return(1);
     }
   return(0);
}

int  chk_trigonal(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3])
{
   int i;
   if( nop!= nop_pgrp[NC3] && nop!= nop_pgrp[NC3i] &&
       nop!= nop_pgrp[ND3] && nop!= nop_pgrp[NC3v] &&
       nop!= nop_pgrp[ND3d]
     ) return(0);

   *lat_new=RHOMBOHEDRAL;
   if( lat == CUBIC_P ) {
       for(i=0; i<nop; i++)
       if( ind[i]==C3p111 ) {
           /*   hexagonal setting  */
           Tnew[0][0]= 1.; Tnew[0][1]= 0.; Tnew[0][2]= 1.;
           Tnew[1][0]=-1.; Tnew[1][1]= 1.; Tnew[1][2]= 1.;
           Tnew[2][0]= 0.; Tnew[2][1]=-1.; Tnew[2][2]= 1.;
           return(1);
       } else
       if( ind[i]==C3p1_1_1 ) { /*  invertion of 2 and 3 strings  */
           Tnew[0][0]= 1.; Tnew[0][1]= 0.; Tnew[0][2]= 1.;
           Tnew[1][0]= 1.; Tnew[1][1]=-1.; Tnew[1][2]=-1.;
           Tnew[2][0]= 0.; Tnew[2][1]= 1.; Tnew[2][2]=-1.;
           return(1);
       } else
       if( ind[i]==C3p_11_1 ) {
           Tnew[0][0]=-1.; Tnew[0][1]= 0.; Tnew[0][2]=-1.;
           Tnew[1][0]=-1.; Tnew[1][1]= 1.; Tnew[1][2]= 1.;
           Tnew[2][0]= 0.; Tnew[2][1]= 1.; Tnew[2][2]=-1.;
           return(1);
       } else
       if( ind[i]==C3p_1_11 ) {
           Tnew[0][0]=-1.; Tnew[0][1]= 0.; Tnew[0][2]=-1.;
           Tnew[1][0]= 1.; Tnew[1][1]=-1.; Tnew[1][2]=-1.;
           Tnew[2][0]= 0.; Tnew[2][1]=-1.; Tnew[2][2]= 1.;
           return(1);
       }
   }
   if( lat == CUBIC_I ) {
       for(i=0; i<nop; i++)
       if( ind[i]==C3p111 ) {
           Tnew[0][0]=-1.; Tnew[0][1]= 0.; Tnew[0][2]= 0.5;
           Tnew[1][0]= 1.; Tnew[1][1]=-1.; Tnew[1][2]= 0.5;
           Tnew[2][0]= 0.; Tnew[2][1]= 1.; Tnew[2][2]= 0.5;
           return(1);
       } else
       if( ind[i]==C3p1_1_1 ) { /*  invertion of 2 and 3 strings  */
           Tnew[0][0]=-1.; Tnew[0][1]= 0.; Tnew[0][2]= 0.5;
           Tnew[1][0]=-1.; Tnew[1][1]= 1.; Tnew[1][2]=-0.5;
           Tnew[2][0]= 0.; Tnew[2][1]=-1.; Tnew[2][2]=-0.5;
           return(1);
       } else
       if( ind[i]==C3p_11_1 ) {
           Tnew[0][0]= 1.; Tnew[0][1]= 0.; Tnew[0][2]=-0.5;
           Tnew[1][0]= 1.; Tnew[1][1]=-1.; Tnew[1][2]= 0.5;
           Tnew[2][0]= 0.; Tnew[2][1]=-1.; Tnew[2][2]=-0.5;
           return(1);
       } else
       if( ind[i]==C3p_1_11 ) {
           Tnew[0][0]= 1.; Tnew[0][1]= 0.; Tnew[0][2]=-0.5;
           Tnew[1][0]=-1.; Tnew[1][1]= 1.; Tnew[1][2]=-0.5;
           Tnew[2][0]= 0.; Tnew[2][1]= 1.; Tnew[2][2]= 0.5;
           return(1);
       }
   }
   if( lat == CUBIC_F ) {
       for(i=0; i<nop; i++)
       if( ind[i]==C3p111 ) {
           Tnew[0][0]=-0.5; Tnew[0][1]= 0.;  Tnew[0][2]= 1.;
           Tnew[1][0]= 0.5; Tnew[1][1]=-0.5; Tnew[1][2]= 1.;
           Tnew[2][0]= 0.;  Tnew[2][1]= 0.5; Tnew[2][2]= 1.;
           return(1);
       } else
       if( ind[i]==C3p1_1_1 ) { /*  invertion of 2 and 3 strings  */
           Tnew[0][0]=-0.5; Tnew[0][1]= 0.;  Tnew[0][2]= 1.;
           Tnew[1][0]=-0.5; Tnew[1][1]= 0.5; Tnew[1][2]=-1.;
           Tnew[2][0]= 0.;  Tnew[2][1]=-0.5; Tnew[2][2]=-1.;
           return(1);
       } else
       if( ind[i]==C3p_11_1 ) {
           Tnew[0][0]= 0.5; Tnew[0][1]= 0.;  Tnew[0][2]=-1.;
           Tnew[1][0]= 0.5; Tnew[1][1]=-0.5; Tnew[1][2]= 1.;
           Tnew[2][0]= 0.;  Tnew[2][1]=-0.5; Tnew[2][2]=-1.;
           return(1);
       } else
       if( ind[i]==C3p_1_11 ) {
           Tnew[0][0]= 0.5; Tnew[0][1]= 0.;  Tnew[0][2]=-1.;
           Tnew[1][0]=-0.5; Tnew[1][1]= 0.5; Tnew[1][2]=-1.;
           Tnew[2][0]= 0.;  Tnew[2][1]= 0.5; Tnew[2][2]= 1.;
           return(1);
       }
   }

   if( lat==HEXAGONAL || lat==RHOMBOHEDRAL ) {
       for(i=0; i<nop; i++)
       if( ind[i]==C3p001h ) {
          *lat_new=lat;
           Tnew[0][0]=1.; Tnew[0][1]=0.; Tnew[0][2]=0.;
           Tnew[1][0]=0.; Tnew[1][1]=1.; Tnew[1][2]=0.;
           Tnew[2][0]=0.; Tnew[2][1]=0.; Tnew[2][2]=1.;
           return(1);
       }
   }
   return(0);
}

int  chk_tetragonal(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3])
{
   int i;
   if( nop!= nop_pgrp[NC4]  && nop!= nop_pgrp[NS4] &&
       nop!= nop_pgrp[NC4h] && nop!= nop_pgrp[ND4] &&
       nop!= nop_pgrp[NC4v] && nop!= nop_pgrp[ND2d] &&
       nop!= nop_pgrp[ND4h]
     ) return(0);
   if( lat==RHOMBOHEDRAL ||
       lat==ORTHOROMBIC_P || lat==ORTHOROMBIC_C ||  lat==ORTHOROMBIC_A ||
       lat==ORTHOROMBIC_I || lat==ORTHOROMBIC_F ||
       lat==MONOCLINIC_P  || lat==MONOCLINIC_A ||
       lat==TRICLINIC  )   return(0);

  *lat_new=lat;
   if( lat==CUBIC_P ) *lat_new=TETRAGONAL_P;
   else
   if( lat==CUBIC_I || lat==CUBIC_F) *lat_new=TETRAGONAL_I;

   for(i=0; i<nop; i++) {
       if( ind[i]==C4p001 || ind[i]==S4p001 ) {
           if( lat==CUBIC_P || lat==CUBIC_I ||
               lat==TETRAGONAL_P ||  lat==TETRAGONAL_I) {
               Tnew[0][0]=1.; Tnew[0][1]=0.; Tnew[0][2]=0.;
               Tnew[1][0]=0.; Tnew[1][1]=1.; Tnew[1][2]=0.;
               Tnew[2][0]=0.; Tnew[2][1]=0.; Tnew[2][2]=1.;
               return(1);
           }
           else
           if( lat==CUBIC_F) {
               Tnew[0][0]= 0.5; Tnew[0][1]= 0.5; Tnew[0][2]=0.;
               Tnew[1][0]=-0.5; Tnew[1][1]= 0.5; Tnew[1][2]=0.;
               Tnew[2][0]= 0.;  Tnew[2][1]= 0.;  Tnew[2][2]=1.;
               return(1);
           }
       }

       if( ind[i]==C4p010 || ind[i]==S4p010 ) {
           if( lat==CUBIC_P || lat==CUBIC_I ) {
               Tnew[0][0]= 0.; Tnew[0][1]= 1.; Tnew[0][2]= 0.;
               Tnew[1][0]= 0.; Tnew[1][1]= 0.; Tnew[1][2]= 1.;
               Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
               return(1);
           }
           else
           if( lat==CUBIC_F) {
               Tnew[0][0]=-0.5; Tnew[0][1]= 0.5; Tnew[0][2]=0.;
               Tnew[1][0]= 0.0; Tnew[1][1]= 0.0; Tnew[1][2]=1.;
               Tnew[2][0]= 0.5; Tnew[2][1]= 0.5; Tnew[2][2]=0.;
               return(1);
           }
       }

       if( ind[i]==C4p100 ||ind[i]==S4p100 ) {
           if( lat==CUBIC_P || lat==CUBIC_I ) {
               Tnew[0][0]= 0.; Tnew[0][1]= 0.; Tnew[0][2]= 1.;
               Tnew[1][0]= 1.; Tnew[1][1]= 0.; Tnew[1][2]= 0.;
               Tnew[2][0]= 0.; Tnew[2][1]= 1.; Tnew[2][2]= 0.;
               return(1);
           }
           else
           if( lat==CUBIC_F) {
               Tnew[0][0]= 0.0; Tnew[0][1]= 0.0; Tnew[0][2]=1.;
               Tnew[1][0]= 0.5; Tnew[1][1]= 0.5; Tnew[1][2]=0.;
               Tnew[2][0]=-0.5; Tnew[2][1]= 0.5; Tnew[2][2]=0.;
               return(1);
           }
       }
   }
   return(0);
}

int  chk_orthorombic(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3])
{

   if( nop!= nop_pgrp[ND2] && nop!= nop_pgrp[NC2v] &&
       nop!= nop_pgrp[ND2h]
     ) return(0);

   if( lat==RHOMBOHEDRAL ||
       lat==MONOCLINIC_P || lat==MONOCLINIC_A ||
       lat==TRICLINIC  ) return(0);
   if( lat==HEXAGONAL ) {
      *lat_new=ORTHOROMBIC_C;
    /* check C2001h */
       if( (ind[1]==C2001h && ind[2]==C2110h && ind[3]==C21_10h ) ||
           (ind[1]==C2001h && ind[2]==Mx_xzh && ind[3]==Mxxzh ) ) {
           Tnew[0][0]= 1.; Tnew[0][1]= 1.; Tnew[0][2]= 0.;
           Tnew[1][0]=-1.; Tnew[1][1]= 1.; Tnew[1][2]= 0.;
           Tnew[2][0]= 0.; Tnew[2][1]= 0.; Tnew[2][2]= 1.;
           return(1);
       }
       if( (ind[1]==C2001h && ind[2]==C2100h && ind[3]==C2120h ) ||
           (ind[1]==C2001h && ind[2]==Mx2xzh && ind[3]==Mx0zh ) ) {
           Tnew[0][0]= 1.; Tnew[0][1]= 1.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 2.; Tnew[1][2]= 0.;
           Tnew[2][0]= 0.; Tnew[2][1]= 0.; Tnew[2][2]= 1.;
           return(1);
       }
       if( (ind[1]==C2001h && ind[2]==C2010h && ind[3]==C2210h) ||
           (ind[1]==C2001h && ind[2]==M2xxzh && ind[3]==M0yzh ) ) {
           Tnew[0][0]= 2.; Tnew[0][1]= 0.; Tnew[0][2]= 0.;
           Tnew[1][0]= 1.; Tnew[1][1]= 1.; Tnew[1][2]= 0.;
           Tnew[2][0]= 0.; Tnew[2][1]= 0.; Tnew[2][2]= 1.;
           return(1);
       }
      *lat_new=ORTHOROMBIC_A;
    /* check C21_10h */
       if( ind[1]==C21_10h && ind[2]==Mxy0h && ind[3]==Mx_xzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]=-1.; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]=-1.; Tnew[1][2]=-1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2100h */
       if( ind[1]==C2100h && ind[2]==Mxy0h && ind[3]==Mx0zh ) {
           Tnew[0][0]= 0.; Tnew[0][1]=-1.; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]=-2.; Tnew[1][2]= 0.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2210h */
       if( ind[1]==C2210h && ind[2]==Mxy0h && ind[3]==M2xxzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 0.; Tnew[0][2]= 2.;
           Tnew[1][0]= 0.; Tnew[1][1]=-1.; Tnew[1][2]= 1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2110h */
       if( ind[1]==C2110h && ind[2]==Mxy0h && ind[3]==Mxxzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 1.; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]=-1.; Tnew[1][2]= 1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2120h */
       if( ind[1]==C2120h && ind[2]==Mxy0h && ind[3]==Mx2xzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 1.; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]= 0.; Tnew[1][2]= 2.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2010h */
       if( ind[1]==C2010h && ind[2]==Mxy0h && ind[3]==M0yzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 2.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 1.; Tnew[1][2]= 1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
   } /*    if( lat==HEXAGONAL )  */

   *lat_new=0;
   /*  check C2001  */
   if( (ind[1]==C2001 && ind[2]==C2010 && ind[3]==C2100) ||
       (ind[1]==C2001 && ind[2]==Mx0z &&  ind[3]==M0yz) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_P;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_I;
       else
       if( lat==CUBIC_F )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==ORTHOROMBIC_P || lat==ORTHOROMBIC_C || lat==ORTHOROMBIC_A ||
           lat==ORTHOROMBIC_I || lat==ORTHOROMBIC_F )  *lat_new=lat;

       if( *lat_new != 0 ) {
           Tnew[0][0]= 1.; Tnew[0][1]= 0.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 1.; Tnew[1][2]= 0.;
           Tnew[2][0]= 0.; Tnew[2][1]= 0.; Tnew[2][2]= 1.;
           return(1);
       }
   }

   if( (ind[1]==C2001 && ind[2]==C2110 && ind[3]==C21_10) ||
       (ind[1]==C2001 && ind[2]==Mx_xz && ind[3]==Mxxz) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_C;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==CUBIC_F ) {
          *lat_new=ORTHOROMBIC_I;
           Tnew[0][0]= 0.5; Tnew[0][1]= 0.5; Tnew[0][2]= 0.;
           Tnew[1][0]=-0.5; Tnew[1][1]= 0.5; Tnew[1][2]= 0.;
           Tnew[2][0]= 0. ; Tnew[2][1]= 0. ; Tnew[2][2]= 1.;
           return(1);
       }

       if( *lat_new != 0 ) {
           Tnew[0][0]= 1.0; Tnew[0][1]= 1.0; Tnew[0][2]= 0.;
           Tnew[1][0]=-1.0; Tnew[1][1]= 1.0; Tnew[1][2]= 0.;
           Tnew[2][0]= 0. ; Tnew[2][1]= 0. ; Tnew[2][2]= 1.;
           return(1);
       }
   }
   /*  check C2010  */
   if( (ind[1]==C2010 && ind[2]==Mxy0 &&  ind[3]==M0yz) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_P;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_I;
       else
       if( lat==CUBIC_F )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==ORTHOROMBIC_C ) *lat_new=ORTHOROMBIC_A;
       else
       if( lat==ORTHOROMBIC_P || lat==ORTHOROMBIC_I ||
           lat==ORTHOROMBIC_F )  *lat_new=lat;

       if( *lat_new != 0 ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 1.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 0.; Tnew[1][2]= 1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
   }

   if( (ind[1]==C2010 && ind[2]==C2101 && ind[3]==C2_101) ||
       (ind[1]==C2010 && ind[2]==M_xyx && ind[3]==Mxyx) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_C;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==CUBIC_F ) {
          *lat_new=ORTHOROMBIC_I;
           Tnew[0][0]=-0.5; Tnew[0][1]= 0.5; Tnew[0][2]= 0.;
           Tnew[1][0]= 0. ; Tnew[1][1]= 0.;  Tnew[1][2]= 1.;
           Tnew[2][0]= 0.5; Tnew[2][1]= 0.5; Tnew[2][2]= 0.;
           return(1);
       }

       if( *lat_new != 0 ) {
           Tnew[0][0]=-1.; Tnew[0][1]= 1.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 0.; Tnew[1][2]= 1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 1.; Tnew[2][2]= 0.;
           return(1);
       }
   }

   /*  check C2100  */
   if( (ind[1]==C2100 && ind[2]==Mxy0 &&  ind[3]==Mx0z) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_P;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_I;
       else
       if( lat==CUBIC_F )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==ORTHOROMBIC_C ) *lat_new=ORTHOROMBIC_A;
       else
       if( lat==ORTHOROMBIC_P || lat==ORTHOROMBIC_I ||
           lat==ORTHOROMBIC_F )  *lat_new=lat;

       if( *lat_new != 0 ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 0.; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]= 1.; Tnew[1][2]= 0.;
           Tnew[2][0]=-1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
   }

   if( (ind[1]==C2100 && ind[2]==C2011 && ind[3]==C201_1) ||
       (ind[1]==C2100 && ind[2]==Mxy_y && ind[3]==Mxyy) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_C;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==CUBIC_F ) {
          *lat_new=ORTHOROMBIC_I;
           Tnew[0][0]= 0. ; Tnew[0][1]= 0. ; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.5; Tnew[1][1]= 0.5; Tnew[1][2]= 0.;
           Tnew[2][0]=-0.5; Tnew[2][1]= 0.5; Tnew[2][2]= 0.;
           return(1);
       }

       if( *lat_new != 0 ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 0.; Tnew[0][2]= 1.;
           Tnew[1][0]= 1.; Tnew[1][1]= 1.; Tnew[1][2]= 0.;
           Tnew[2][0]=-1.; Tnew[2][1]= 1.; Tnew[2][2]= 0.;
           return(1);
       }
   }

   /*  check C2110  */
   if( (ind[1]==C2110 && ind[2]==Mxy0 &&  ind[3]==Mxxz) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_A;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==CUBIC_F ) {
          *lat_new=ORTHOROMBIC_I;
           Tnew[0][0]= 0.; Tnew[0][1]= 0.5; Tnew[0][2]= 0.5;
           Tnew[1][0]= 0.; Tnew[1][1]=-0.5; Tnew[1][2]= 0.5;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.;  Tnew[2][2]= 0.;
           return(1);
       }
       if( *lat_new != 0 ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 1.; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]=-1.; Tnew[1][2]= 1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
   }

   /*  check C21_10  */
   if( (ind[1]==C21_10 && ind[2]==Mxy0 &&  ind[3]==Mx_xz) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_A;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==CUBIC_F ) {
          *lat_new=ORTHOROMBIC_I;
           Tnew[0][0]= 0.0; Tnew[0][1]=-0.5; Tnew[0][2]= 0.5;
           Tnew[1][0]= 0.0; Tnew[1][1]=-0.5; Tnew[1][2]=-0.5;
           Tnew[2][0]= 1.0; Tnew[2][1]= 0.0; Tnew[2][2]= 0.0;
           return(1);
       }
       if( *lat_new != 0 ) {
           Tnew[0][0]= 0.0; Tnew[0][1]=-1.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 0.0; Tnew[1][1]=-1.0; Tnew[1][2]=-1.0;
           Tnew[2][0]= 1.0; Tnew[2][1]= 0.0; Tnew[2][2]= 0.0;
           return(1);
       }
   }

   /*  check C2101  */
   if( (ind[1]==C2101 && ind[2]==Mx0z &&  ind[3]==Mxyx) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_A;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==CUBIC_F ) {
          *lat_new=ORTHOROMBIC_I;
           Tnew[0][0]= 0.; Tnew[0][1]=-0.5; Tnew[0][2]= 0.5;
           Tnew[1][0]= 1.; Tnew[1][1]= 0.;  Tnew[1][2]= 0.;
           Tnew[2][0]= 0.; Tnew[2][1]= 0.5; Tnew[2][2]= 0.5;
           return(1);
       }
       if( *lat_new != 0 ) {
           Tnew[0][0]= 0.; Tnew[0][1]=-1.; Tnew[0][2]= 1.;
           Tnew[1][0]= 1.; Tnew[1][1]= 0.; Tnew[1][2]= 0.;
           Tnew[2][0]= 0.; Tnew[2][1]= 1.; Tnew[2][2]= 1.;
           return(1);
       }
   }
   /*  check C2_101  */
   if( (ind[1]==C2_101 && ind[2]==Mx0z &&  ind[3]==M_xyx) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_A;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==CUBIC_F ) {
          *lat_new=ORTHOROMBIC_I;
           Tnew[0][0]= 0.0; Tnew[0][1]=-0.5; Tnew[0][2]=-0.5;
           Tnew[1][0]= 1.0; Tnew[1][1]= 0.0; Tnew[1][2]= 0.0;
           Tnew[2][0]= 0.0; Tnew[2][1]=-0.5; Tnew[2][2]= 0.5;
           return(1);
       }
       if( *lat_new != 0 ) {
           Tnew[0][0]= 0.0; Tnew[0][1]=-1.0; Tnew[0][2]=-1.0;
           Tnew[1][0]= 1.0; Tnew[1][1]= 0.0; Tnew[1][2]= 0.0;
           Tnew[2][0]= 0.0; Tnew[2][1]=-1.0; Tnew[2][2]= 1.0;
           return(1);
       }
   }

   /*  check C2011  */
   if( (ind[1]==C2011 && ind[2]==M0yz &&  ind[3]==Mxyy) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_A;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==CUBIC_F ) {
          *lat_new=ORTHOROMBIC_I;
           Tnew[0][0]= 1.; Tnew[0][1]= 0.;  Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 0.5; Tnew[1][2]= 0.5;
           Tnew[2][0]= 0.; Tnew[2][1]=-0.5; Tnew[2][2]= 0.5;
           return(1);
       }
       if( *lat_new != 0 ) {
           Tnew[0][0]= 1.; Tnew[0][1]= 0.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 1.; Tnew[1][2]= 1.;
           Tnew[2][0]= 0.; Tnew[2][1]=-1.; Tnew[2][2]= 1.;
           return(1);
       }
   }

   /*  check C201_1  */
   if( (ind[1]==C201_1 && ind[2]==M0yz &&  ind[3]==Mxy_y) ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P )
          *lat_new=ORTHOROMBIC_A;
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I )
          *lat_new=ORTHOROMBIC_F;
       else
       if( lat==CUBIC_F ) {
          *lat_new=ORTHOROMBIC_I;
           Tnew[0][0]= 1.0; Tnew[0][1]= 0.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]=-0.5; Tnew[1][2]= 0.5;
           Tnew[2][0]= 0.0; Tnew[2][1]=-0.5; Tnew[2][2]=-0.5;
           return(1);
       }
       if( *lat_new != 0 ) {
           Tnew[0][0]= 1.0; Tnew[0][1]= 0.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]=-1.0; Tnew[1][2]= 1.0;
           Tnew[2][0]= 0.0; Tnew[2][1]=-1.0; Tnew[2][2]=-1.0;
           return(1);
       }
   }
   return(0);
}

int  chk_monoclinic(int lat, int nop, int *ind, int *lat_new, double Tnew[3][3])
{
   if( nop != nop_pgrp[NC2] && nop != nop_pgrp[NCs] &&
       nop != nop_pgrp[NC2h]
     ) return(0);
   if( lat==TRICLINIC  ) return(0);

   if( lat==HEXAGONAL ) {
      *lat_new=MONOCLINIC_P;
    /* check C2001h */
       if( ind[1]==C2001h || ind[1]==Mxy0h ) {
           Tnew[0][0]= 1.; Tnew[0][1]= 0.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 1.; Tnew[1][2]= 0.;
           Tnew[2][0]= 0.; Tnew[2][1]= 0.; Tnew[2][2]= 1.;
           return(1);
       }
      *lat_new=MONOCLINIC_A; /*  axes like in ORTHOROMBIC  */
    /* check C21_10h */
       if( ind[1]==C21_10h || ind[1]==Mxxzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]=-1.; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]=-1.; Tnew[1][2]=-1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2100h */
       if( ind[1]==C2100h || ind[1]==Mx2xzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]=-1.; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]=-2.; Tnew[1][2]= 0.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2210h */
       if( ind[1]==C2210h || ind[1]==M0yzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 0.; Tnew[0][2]= 2.;
           Tnew[1][0]= 0.; Tnew[1][1]=-1.; Tnew[1][2]= 1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2110h */
       if( ind[1]==C2110h || ind[1]==Mx_xzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 1.; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]=-1.; Tnew[1][2]= 1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2120h */
       if( ind[1]==C2120h || ind[1]==Mx0zh ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 1.; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]= 0.; Tnew[1][2]= 2.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2010h */
       if( ind[1]==C2010h || ind[1]==M2xxzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 2.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 1.; Tnew[1][2]= 1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 0.; Tnew[2][2]= 0.;
           return(1);
       }
   } /*    if( lat==HEXAGONAL )  */

   if( lat==RHOMBOHEDRAL ) { /*  input basis is HEXAGONAL  <>*/
      *lat_new=MONOCLINIC_A;
    /* check C2100h */
       if( ind[1]==C2100h || ind[1]==Mx2xzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 1./3;  Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]= 2./3;  Tnew[1][2]= 0.;
           Tnew[2][0]=-1.; Tnew[2][1]= 2./3;  Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2110h */
       if( ind[1]==C2110h || ind[1]==Mx_xzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]= 1./3; Tnew[0][2]= 1.;
           Tnew[1][0]= 0.; Tnew[1][1]=-1./3; Tnew[1][2]= 1.;
           Tnew[2][0]= 1.; Tnew[2][1]= 2./3; Tnew[2][2]= 0.;
           return(1);
       }
    /* check C2010h */
       if( ind[1]==C2010h || ind[1]==M2xxzh ) {
           Tnew[0][0]= 0.; Tnew[0][1]=-2./3; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]=-1./3; Tnew[1][2]= 1.;
           Tnew[2][0]=-1.; Tnew[2][1]= 2./3; Tnew[2][2]= 0.;
           return(1);
       }
   } /*    if( lat==RHOMBOHEDRAL )  */

   /*  check C2001  */
   if( ind[1]==C2001 || ind[1]==Mxy0 ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P || lat==ORTHOROMBIC_P ||
           lat==MONOCLINIC_P ) {
          *lat_new=MONOCLINIC_P;
           Tnew[0][0]= 1.; Tnew[0][1]= 0.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 1.; Tnew[1][2]= 0.;
           Tnew[2][0]= 0.; Tnew[2][1]= 0.; Tnew[2][2]= 1.;
           return(1);
       }
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I || lat==ORTHOROMBIC_I) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 1.; Tnew[0][1]=-1.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 1.; Tnew[1][2]= 0.;
           Tnew[2][0]= 0.; Tnew[2][1]= 0.; Tnew[2][2]= 1.;
           return(1);
       }
       else
       if( lat==CUBIC_F || lat==ORTHOROMBIC_F) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.5; Tnew[0][1]=-1.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.5; Tnew[1][1]= 0.0; Tnew[1][2]= 0.0;
           Tnew[2][0]= 0.0; Tnew[2][1]= 0.0; Tnew[2][2]= 1.0;
           return(1);
       }
       else
       if( lat==ORTHOROMBIC_C ) {
          *lat_new=MONOCLINIC_P;
           Tnew[0][0]= 0.5; Tnew[0][1]=-1.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.5; Tnew[1][1]= 0.0; Tnew[1][2]= 0.0;
           Tnew[2][0]= 0.0; Tnew[2][1]= 0.0; Tnew[2][2]= 1.0;
           return(1);
       }
       if( lat==ORTHOROMBIC_A || lat==MONOCLINIC_A ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 1.; Tnew[0][1]= 0.; Tnew[0][2]= 0.;
           Tnew[1][0]= 0.; Tnew[1][1]= 1.; Tnew[1][2]= 0.;
           Tnew[2][0]= 0.; Tnew[2][1]= 0.; Tnew[2][2]= 1.;
           return(1);
       }
   }
   /*  check C2010; in this case lat != MONOCLINIC  */
   if( ind[1]==C2010 || ind[1]==Mx0z ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P || lat==ORTHOROMBIC_P ) {
          *lat_new=MONOCLINIC_P;
           Tnew[0][0]= 1.0; Tnew[0][1]= 0.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]= 0.0; Tnew[1][2]= 1.0;
           Tnew[2][0]= 0.0; Tnew[2][1]=-1.0; Tnew[2][2]= 0.0;
           return(1);
       }
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I || lat==ORTHOROMBIC_I) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 1.0; Tnew[0][1]=-1.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]= 0.0; Tnew[1][2]= 1.0;
           Tnew[2][0]= 0.0; Tnew[2][1]=-1.0; Tnew[2][2]= 0.0;
           return(1);
       }
       else
       if( lat==CUBIC_F || lat==ORTHOROMBIC_F) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.5; Tnew[0][1]=-1.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]= 0.0; Tnew[1][2]= 1.0;
           Tnew[2][0]=-0.5; Tnew[2][1]= 0.0; Tnew[2][2]= 0.0;
           return(1);
       }
       else
       if( lat==ORTHOROMBIC_C ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]= 1.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]= 0.0; Tnew[1][2]= 1.0;
           Tnew[2][0]= 1.0; Tnew[2][1]= 0.0; Tnew[2][2]= 0.0;
           return(1);
       }
       if( lat==ORTHOROMBIC_A ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 1.0; Tnew[0][1]= 0.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]= 0.0; Tnew[1][2]= 1.0;
           Tnew[2][0]= 0.0; Tnew[2][1]=-1.0; Tnew[2][2]= 0.0;
           return(1);
       }
   }
   /*  check C2100; in this case lat != MONOCLINIC  */
   if( ind[1]==C2100 || ind[1]==M0yz ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P || lat==ORTHOROMBIC_P ) {
          *lat_new=MONOCLINIC_P;
           Tnew[0][0]= 0.0; Tnew[0][1]= 0.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 0.0; Tnew[1][1]= 1.0; Tnew[1][2]= 0.0;
           Tnew[2][0]=-1.0; Tnew[2][1]= 0.0; Tnew[2][2]= 0.0;
           return(1);
       }
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I || lat==ORTHOROMBIC_I) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]= 0.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 0.0; Tnew[1][1]= 1.0; Tnew[1][2]= 0.0;
           Tnew[2][0]=-1.0; Tnew[2][1]= 1.0; Tnew[2][2]= 0.0;
           return(1);
       }
       else
       if( lat==CUBIC_F || lat==ORTHOROMBIC_F) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]= 0.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 0.5; Tnew[1][1]= 0.0; Tnew[1][2]= 0.0;
           Tnew[2][0]=-0.5; Tnew[2][1]= 1.0; Tnew[2][2]= 0.0;
           return(1);
       }
       else
       if( lat==ORTHOROMBIC_C ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]= 0.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 0.0; Tnew[1][1]= 1.0; Tnew[1][2]= 0.0;
           Tnew[2][0]=-1.0; Tnew[2][1]= 0.0; Tnew[2][2]= 0.0;
           return(1);
       }
       if( lat==ORTHOROMBIC_A ) {
          *lat_new=MONOCLINIC_P;
           Tnew[0][0]= 0.0; Tnew[0][1]= 0.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 0.5; Tnew[1][1]= 0.0; Tnew[1][2]= 0.0;
           Tnew[2][0]=-0.5; Tnew[2][1]= 1.0; Tnew[2][2]= 0.0;
           return(1);
       }
   }
   /*  check C2110; in this case lat != MONOCLINIC != ORTHOROMBIC  */
   if( ind[1]==C2110 || ind[1]==Mx_xz ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]= 1.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 0.0; Tnew[1][1]=-1.0; Tnew[1][2]= 1.0;
           Tnew[2][0]= 1.0; Tnew[2][1]= 0.0; Tnew[2][2]= 0.0;
           return(1);
       }
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]=-0.5; Tnew[0][1]= 0.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 0.5; Tnew[1][1]= 0.0; Tnew[1][2]= 1.0;
           Tnew[2][0]=-0.5; Tnew[2][1]= 1.0; Tnew[2][2]= 0.0;
           return(1);
       }
       else
       if( lat==CUBIC_F ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]= 0.5; Tnew[0][2]= 0.5;
           Tnew[1][0]= 0.0; Tnew[1][1]=-0.5; Tnew[1][2]= 0.5;
           Tnew[2][0]= 1.0; Tnew[2][1]=-1.0; Tnew[2][2]= 0.0;
           return(1);
       }
   }
   /*  check C21_10; in this case lat != MONOCLINIC != ORTHOROMBIC  */
   /*  mul. C2110 on C4_001; x'=y  y'=-x  z'=z  */
   if( ind[1]==C21_10 || ind[1]==Mxxz ) {
       if( lat==CUBIC_P || lat==TETRAGONAL_P ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]=-1.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 0.0; Tnew[1][1]=-1.0; Tnew[1][2]=-1.0;
           Tnew[2][0]= 1.0; Tnew[2][1]= 0.0; Tnew[2][2]= 0.0;
           return(1);
       }
       else
       if( lat==CUBIC_I || lat==TETRAGONAL_I ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.5; Tnew[0][1]= 0.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 0.5; Tnew[1][1]= 0.0; Tnew[1][2]=-1.0;
           Tnew[2][0]=-0.5; Tnew[2][1]= 1.0; Tnew[2][2]= 0.0;
           return(1);
       }
       else
       if( lat==CUBIC_F ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]=-0.5; Tnew[0][2]= 0.5;
           Tnew[1][0]= 0.0; Tnew[1][1]=-0.5; Tnew[1][2]=-0.5;
           Tnew[2][0]= 1.0; Tnew[2][1]=-1.0; Tnew[2][2]= 0.0;
           return(1);
       }
   }
   /*  check C2101; in this case lat != MONOCLINIC != ORTHOROMBIC != TETRAGONAL  */
   /*   from C2110 mul. C3-x,x,x;   x'=y  y'=z  z'=x  */
   if( ind[1]==C2101 || ind[1]==M_xyx ) {
       if( lat==CUBIC_P ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]=-1.0; Tnew[0][2]= 1.0;
           Tnew[1][0]= 1.0; Tnew[1][1]= 0.0; Tnew[1][2]= 0.0;
           Tnew[2][0]= 0.0; Tnew[2][1]= 1.0; Tnew[2][2]= 1.0;
           return(1);
       }
       else
       if( lat==CUBIC_I ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.5; Tnew[0][1]= 0.0; Tnew[0][2]= 1.0;
           Tnew[1][0]=-0.5; Tnew[1][1]= 1.0; Tnew[1][2]= 0.0;
           Tnew[2][0]=-0.5; Tnew[2][1]= 0.0; Tnew[2][2]= 1.0;
           return(1);
       }
       else
       if( lat==CUBIC_F ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]=-0.5; Tnew[0][2]= 0.5;
           Tnew[1][0]= 1.0; Tnew[1][1]=-1.0; Tnew[1][2]= 0.0;
           Tnew[2][0]= 0.0; Tnew[2][1]= 0.5; Tnew[2][2]= 0.5;
           return(1);
       }
   }
   /*  check C2_101; in this case lat != MONOCLINIC != ORTHOROMBIC != TETRAGONAL  */
   /*   4-0y0;  x'=-z  y'=y  z'=x  */
   if( ind[1]==C2_101 || ind[1]==Mxyx ) {
       if( lat==CUBIC_P ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]=-1.0; Tnew[0][2]=-1.0;
           Tnew[1][0]= 1.0; Tnew[1][1]= 0.0; Tnew[1][2]= 0.0;
           Tnew[2][0]= 0.0; Tnew[2][1]=-1.0; Tnew[2][2]= 1.0;
           return(1);
       }
       else
       if( lat==CUBIC_I ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.5; Tnew[0][1]= 0.0; Tnew[0][2]=-1.0;
           Tnew[1][0]=-0.5; Tnew[1][1]= 1.0; Tnew[1][2]= 0.0;
           Tnew[2][0]= 0.5; Tnew[2][1]= 0.0; Tnew[2][2]= 1.0;
           return(1);
       }
       else
       if( lat==CUBIC_F ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 0.0; Tnew[0][1]=-0.5; Tnew[0][2]=-0.5;
           Tnew[1][0]= 1.0; Tnew[1][1]=-1.0; Tnew[1][2]= 0.0;
           Tnew[2][0]= 0.0; Tnew[2][1]=-0.5; Tnew[2][2]= 0.5;
           return(1);
       }
   }
   /*  check C2011; in this case lat != MONOCLINIC != ORTHOROMBIC != TETRAGONAL  */
   /*  from C2110 mul. C3+x,x,x;   x'=z  y'=x  z'=y  */
   if( ind[1]==C2011 || ind[1]==Mxy_y ) {
       if( lat==CUBIC_P ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 1.0; Tnew[0][1]= 0.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]= 1.0; Tnew[1][2]= 1.0;
           Tnew[2][0]= 0.0; Tnew[2][1]=-1.0; Tnew[2][2]= 1.0;
           return(1);
       }
       else
       if( lat==CUBIC_I ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]=-0.5; Tnew[0][1]= 1.0; Tnew[0][2]= 0.0;
           Tnew[1][0]=-0.5; Tnew[1][1]= 0.0; Tnew[1][2]= 1.0;
           Tnew[2][0]= 0.5; Tnew[2][1]= 0.0; Tnew[2][2]= 1.0;
           return(1);
       }
       else
       if( lat==CUBIC_F ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 1.0; Tnew[0][1]=-1.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]= 0.5; Tnew[1][2]= 0.5;
           Tnew[2][0]= 0.0; Tnew[2][1]=-0.5; Tnew[2][2]= 0.5;
           return(1);
       }
   }
   /* check C201_1; in this case lat != MONOCLINIC != ORTHOROMBIC != TETRAGONAL
    * mul. C4_100; x,y,z -> x,z,-y;  x'=x  y'=z  z'=-y */
   if( ind[1]==C201_1 || ind[1]==Mxyy ) {
       if( lat==CUBIC_P ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 1.0; Tnew[0][1]= 0.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]=-1.0; Tnew[1][2]= 1.0;
           Tnew[2][0]= 0.0; Tnew[2][1]=-1.0; Tnew[2][2]=-1.0;
           return(1);
       }
       else
       if( lat==CUBIC_I ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]=-0.5; Tnew[0][1]= 1.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.5; Tnew[1][1]= 0.0; Tnew[1][2]= 1.0;
           Tnew[2][0]= 0.5; Tnew[2][1]= 0.0; Tnew[2][2]=-1.0;
           return(1);
       }
       else
       if( lat==CUBIC_F ) {
          *lat_new=MONOCLINIC_A;
           Tnew[0][0]= 1.0; Tnew[0][1]=-1.0; Tnew[0][2]= 0.0;
           Tnew[1][0]= 0.0; Tnew[1][1]=-0.5; Tnew[1][2]= 0.5;
           Tnew[2][0]= 0.0; Tnew[2][1]=-0.5; Tnew[2][2]=-0.5;
           return(1);
       }
   }
   return(0);
}

int  chk_triclinic(int lat, int nop, int *ind, double Ttrcl[3][3],
                   int *lat_new, double Tnew[3][3])
{
   if( nop!= nop_pgrp[NC1] && nop!= nop_pgrp[NCi] )  return(0);
   *lat_new=TRICLINIC;
   asgn_n(Tnew[0],Ttrcl[0],9);
   return(1);
}


int new_transl(int lat, int nop, int *ind, double Ttrcl[3][3],
               int *lat_new, double Tnew[3][3])
{
   sorti(ind,nop);
   if( lat==CUBIC_P || lat==CUBIC_I || lat==CUBIC_F ) {
       if( chk_cubic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_trigonal(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_tetragonal(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_orthorombic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_monoclinic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_triclinic(lat,nop,ind,Ttrcl,lat_new,Tnew) ) return 0;
       fprintf(stderr,"Error lat==CUBIC in det_pgrp()\n");
       return 1;
   }

   if( lat==HEXAGONAL ) {
       if( chk_hexagonal(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_trigonal(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_orthorombic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_monoclinic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_triclinic(lat,nop,ind,Ttrcl,lat_new,Tnew) ) return 0;
       fprintf(stderr,"Error lat==HEXAGONAL in det_pgrp()\n");
       return 1;
   }

   if( lat==RHOMBOHEDRAL ) {
       if( chk_trigonal(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_monoclinic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_triclinic(lat,nop,ind,Ttrcl,lat_new,Tnew) ) return 0;
       fprintf(stderr,"Error lat==RHOMBOHEDRAL in det_pgrp()\n");
       return 1;
   }

   if( lat==TETRAGONAL_P || lat==TETRAGONAL_I) {
       if( chk_tetragonal(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_orthorombic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_monoclinic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_triclinic(lat,nop,ind,Ttrcl,lat_new,Tnew) ) return 0;
       fprintf(stderr,"Error lat==TETRAGONAL in det_pgrp()\n");
       return 1;
   }

   /* ORTHOROMBIC_A only for find_atom_pgroup() */
   if( lat==ORTHOROMBIC_P || lat==ORTHOROMBIC_I || lat==ORTHOROMBIC_A ||
       lat==ORTHOROMBIC_C || lat==ORTHOROMBIC_F ) {
       if( chk_orthorombic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_monoclinic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_triclinic(lat,nop,ind,Ttrcl,lat_new,Tnew) ) return 0;
       fprintf(stderr,"Error lat==ORTHOROMBIC in det_pgrp()\n");
       return 1;
   }

   if( lat==MONOCLINIC_P || lat==MONOCLINIC_A ) {
       if( chk_monoclinic(lat,nop,ind,lat_new,Tnew) ) return 0;
       else
       if( chk_triclinic(lat,nop,ind,Ttrcl,lat_new,Tnew) ) return 0;
       fprintf(stderr,"Error lat==MONOCLINIC in det_pgrp()\n");
       return 1;
   }

   if( lat==TRICLINIC ) {
       if( chk_triclinic(lat,nop,ind,Ttrcl,lat_new,Tnew) ) return 0;
       fprintf(stderr,"Error lat==MONOCLINIC in det_pgrp()\n");
       return 1;
   }

   return 0;
}

int is_pgrp_eq(int nop_sym, double sym_op[][4][3],
               int nop, double pgrp[][3][3])
{
   int i,j;
   double r[48][3];

   if( nop_sym != nop ) return(0);
   for(i=0; i<nop_sym; i++) {
     for(j=0; j<nop; j++)
     if( ddvec_n(sym_op[i][0],pgrp[j][0],9) < TOL ) goto FOUND;
     return(0); /*  not found  */
     FOUND:
     asgn(r[j],sym_op[i][3]);
   }
   for(i=0; i<nop; i++) {
      asgn_n(sym_op[i][0],pgrp[i][0],9);
      asgn(sym_op[i][3],r[i]);
   }
   return(1);
}

int det_pgrp(int nop, double sym_op[][4][3])
{
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC1],C1_pgrp) ) return(NC1);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NCi],Ci_pgrp) ) return(NCi);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC2],C2_pgrp) ) return(NC2);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NCs],Cs_pgrp) ) return(NCs);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC2h],C2h_pgrp) ) return(NC2h);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND2],D2_pgrp) ) return(ND2);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC2v],C2v_pgrp) ) return(NC2v);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND2h],D2h_pgrp) ) return(ND2h);
   else

   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC4],C4_pgrp) ) return(NC4);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NS4],S4_pgrp) ) return(NS4);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC4h],C4h_pgrp) ) return(NC4h);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND4],D4_pgrp) ) return(ND4);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC4v],C4v_pgrp) ) return(NC4v);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND2d],D2d_pgrp) ) return(ND2d);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND2d2],D2d2_pgrp) ) return(ND2d2);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND4h],D4h_pgrp) ) return(ND4h);
   else

   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC3],C3_pgrp) ) return(NC3);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC3i],C3i_pgrp) ) return(NC3i);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND3],D3_pgrp) ) return(ND3);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND32],D32_pgrp) ) return(ND32);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC3v],C3v_pgrp) ) return(NC3v);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC3v2],C3v2_pgrp) ) return(NC3v2);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND3d],D3d_pgrp) ) return(ND3d);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND3d2],D3d2_pgrp) ) return(ND3d2);
   else

   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC6],C6_pgrp) ) return(NC6);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC3h],C3h_pgrp) ) return(NC3h);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC6h],C6h_pgrp) ) return(NC6h);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND6],D6_pgrp) ) return(ND6);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NC6v],C6v_pgrp) ) return(NC6v);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND3h],D3h_pgrp) ) return(ND3h);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND3h2],D3h2_pgrp) ) return(ND3h2);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[ND6h],D6h_pgrp) ) return(ND6h);
   else

   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NT],T_pgrp) ) return(NT);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NTh],Th_pgrp) ) return(NTh);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NO],O_pgrp) ) return(NO);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NTd],Td_pgrp) ) return(NTd);
   else
   if( is_pgrp_eq(nop,sym_op,nop_pgrp[NOh],Oh_pgrp) ) return(NOh);

   return(0);
}


