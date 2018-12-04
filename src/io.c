#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "sto.h"
#include "type_sg.h"
#include "math_sg.h"
#include "pgrp_dat.h"
#include "wien.h"
#include "sgroup.h"
#include "io.h"

#define ERROR { fprintf(stderr,"error: Not completed input!\n"); fclose(fl); return 1; }  


/* read data from WIEN input file
 * 
 * to enforce additional nonequivalency between atoms the scheme 
 * proposed by Peter Blaha is implemented:
 *  1. atoms with different Z (nuclear charge) values can not be equivalent
 *  2. optional third character in the name of atom may be a digit
 *     which acts like a switch to distinguish between atoms with equal Z */

/* maximal number of digits behind of atom name, 3-10 positions */
#define MAX_DGT 8

int read_WIEN_data(char *fname, t_cell *cl_in, 
              char name_srt[MAX_ATOMS][MAX_CHARS],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
              int *lat, t_WIEN *WIEN)
{
   FILE *fl;
   double d,d1,s11,s22,s33,s12,s13,s23;
   char s[256],st[64];
   int i,j,k,l,nat_old=0, is_wien_rhomb=0, MULT_W;
   
   if( (fl=fopen(fname,"rt")) == NULL ) {
      fprintf(stderr,"error: Can't open file: %s.\n",fname);
      return 1;
   }
   
  is_wien_rhomb=0;
  /*  read title  */
  if( fgets(WIEN->title,256,fl)==NULL ) ERROR
     n_to_zero(WIEN->title);
  /*  read lattice type  */
  if( fgets(s,256,fl)==NULL ) ERROR
  strncpy(st,s,4); st[4]='\0';
  if( !strcmp(st,"P   ") || !strcmp(st," P  ") || 
      !strcmp(st,"  P ") || !strcmp(st,"   P") ) *lat=ORTHOROMBIC_P;
  else
  if( !strcmp(st,"F   ") || !strcmp(st," F  ") || 
      !strcmp(st,"  F ") || !strcmp(st,"   F") )  *lat=ORTHOROMBIC_F;
  else
  if( !strcmp(st,"B   ") || !strcmp(st," B  ") || 
      !strcmp(st,"  B ") || !strcmp(st,"   B") ) *lat=ORTHOROMBIC_I;
  else
  if( !strcmp(st,"CXY ") || !strcmp(st," CXY") ) *lat=ORTHOROMBIC_C;
  else
  if( !strcmp(st,"CYZ ") || !strcmp(st," CYZ") ) *lat=ORTHOROMBIC_A;
  else
  if( !strcmp(st,"CXZ ") || !strcmp(st," CXZ") ) *lat=ORTHOROMBIC_B;
  else
  if( !strcmp(st,"R   ") || !strcmp(st," R  ") || 
      !strcmp(st,"  R ") || !strcmp(st,"   R") )
      { *lat=RHOMBOHEDRAL; is_wien_rhomb=1; }
  else
  if( !strcmp(st,"H   ") || !strcmp(st," H  ") || 
      !strcmp(st,"  H ") || !strcmp(st,"   H") ) *lat=HEXAGONAL;
  else {
     fprintf(stderr, "error: Can't resolve the lattice type given"
                     " in WIEN struct file.\n");
     return 1;
  }
  /*  WIEN->cl.nsrt=NATO param. of WIEN  */
  strncpy(st,s+4+23,3); st[3]='\0';
  WIEN->cl.nsrt=atoi(st);
  /*  read a,b,c alpha, beta, gamma  */
  if( fgets(s,256,fl)==NULL ) ERROR   /*   mode of calc  */
  n_to_zero(s);
  strncpy(WIEN->mode,s,256);
     
  if( fgets(s,256,fl)==NULL ) ERROR
  n_to_zero(s);
  j=strlen(s);
  strncpy(st,s,10);    st[10]='\0';  a[0][0]=atof(st);
  strncpy(st,s+10,10); st[10]='\0';  a[0][1]=atof(st);
  strncpy(st,s+20,10); st[10]='\0';  a[0][2]=atof(st);
  st[0]='\0'; 
  if( j > 30 ) strncpy(st,s+30,10); st[10]='\0';  
  if( stod(st,&(a[1][0]),3) != 1 )  a[1][0]=0.; /* empty */
  st[0]='\0';
  if( j > 40 ) strncpy(st,s+40,10); st[10]='\0';
  if( stod(st,&(a[1][1]),2) != 1 )  a[1][1]=0.; /* empty */
  st[0]='\0'; 
  if( j > 50 ) strncpy(st,s+50,10); st[10]='\0';
  if( stod(st,&(a[1][2]),1) != 1 )  a[1][2]=0.; /* empty */
  
  if( eq(a[1][0],0.) ) a[1][0]=90.; /* 0 -> 90  */
  if( eq(a[1][1],0.) ) a[1][1]=90.; /* 0 -> 90 */
  if( eq(a[1][2],0.) ) a[1][2]=90.; /* 0 -> 90 */
  if( *lat == RHOMBOHEDRAL || *lat == HEXAGONAL ) a[1][2]=120.;
   
  /* check the lattice and angles conformance
   * F,I,C,A-centered : alpha=beta=gamma=90
   * B-centered       : alpha = beta = 90
   * H,R              : alpha=beta=90; gamma=120
   * P                : none */
  switch( *lat ) {
   case ORTHOROMBIC_F:
   case ORTHOROMBIC_I:
   case ORTHOROMBIC_C:
   case ORTHOROMBIC_A:     
     if( ne(a[1][0],90.) ) {
       fprintf(stderr,"error: alpha = %f  and not equal 90. "
                      "Exiting now.\n",a[1][0]);
       return 1;
     } else
     if( ne(a[1][1],90.) ) {
       fprintf(stderr,"error: beta = %f  and not equal 90. "
                      "Exiting now.\n",a[1][1]);
       return 1;
     } else
     if( ne(a[1][2],90.) ) {
       fprintf(stderr,"error: gamma = %f  and not equal 90. "
                       "Exiting now.\n",a[1][2]);
       return 1;
     }
     break;

   case ORTHOROMBIC_B:
     if( ne(a[1][0],90.) ) {
       fprintf(stderr,"error: alpha = %f  and not equal 90. "
                      "Exiting now.\n",a[1][0]);
       return 1;
     } else
     if( ne(a[1][1],90.) ) {
       fprintf(stderr,"error: beta = %f  and not equal 90. "
                      "Exiting now.\n",a[1][1]);
       return 1;
     }
     break;
     
   case RHOMBOHEDRAL:
   case HEXAGONAL:
     if( ne(a[1][0],90.) ) {
        fprintf(stderr,"error: alpha = %f  and not equal 90. "
                      "Exiting now.\n",a[1][0]);
        return 1;
     } else
     if( ne(a[1][1],90.) ) {
       fprintf(stderr,"error: beta = %f  and not equal 90. "
                      "Exiting now.\n",a[1][1]);
       return 1;
     } else
     if( ne(a[1][2],120.) ) {
       fprintf(stderr,"error: gamma = %f  and not equal 120. "
                      "Exiting now.\n",a[1][2]);
       return 1;
     }
     break;
  }

  asgn_n(WIEN->a[0],a[0],2*3);
  WIEN->lat = *lat;
  
  /* s11,s22,... only for centered modes */ 
  s11=a[0][0]*a[0][0]; s22=a[0][1]*a[0][1]; s33=a[0][2]*a[0][2];
  s12=a[0][0]*a[0][1]*cos(PI/180.*a[1][2]);
  s13=0; /* beta=90 a[0][0]*a[0][2]*cos(PI/180.*a[1][1]);*/
  s23=0; /* alpha=90 */
   
  switch( *lat ) {
   case ORTHOROMBIC_I:
      a[0][0]=sqrt( s11 + s22 + s33 - 2.*s12 - 2.*s13 + 2.*s23 )/2.;
      a[0][1]=sqrt( s11 + s22 + s33 - 2.*s12 + 2.*s13 - 2.*s23 )/2.;
      a[0][2]=sqrt( s11 + s22 + s33 + 2.*s12 - 2.*s13 - 2.*s23 )/2.;
      a[1][0]=180./PI*acos(( s11-s22-s33+2.*s23)/(4.*a[0][1]*a[0][2]));
      a[1][1]=180./PI*acos((-s11+s22-s33+2.*s13)/(4.*a[0][0]*a[0][2]));
      a[1][2]=180./PI*acos((-s11-s22+s33+2.*s12)/(4.*a[0][0]*a[0][1]));
      break;
   case ORTHOROMBIC_F:
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1]; a[0][2]*=a[0][2];
      a[1][0]=180./PI*acos( a[0][0]/sqrt((a[0][0]+a[0][1])*(a[0][0]+a[0][2])) );
      a[1][1]=180./PI*acos( a[0][1]/sqrt((a[0][0]+a[0][1])*(a[0][1]+a[0][2])) );
      a[1][2]=180./PI*acos( a[0][2]/sqrt((a[0][0]+a[0][2])*(a[0][1]+a[0][2])) );
      d      =sqrt(a[0][1]+a[0][2])/2;
      d1     =sqrt(a[0][0]+a[0][2])/2;
      a[0][2]=sqrt(a[0][0]+a[0][1])/2;
      a[0][0]=d; a[0][1]=d1;
      break;
   case ORTHOROMBIC_C:
      a[0][0]=sqrt( s11 + s22 - 2.*s12 )/2.;
      a[0][1]=sqrt( s11 + s22 + 2.*s12 )/2.;
      /* a[0][2]=a[0][2]; */
      a[1][0]=180./PI*acos((s13+s23)/(2.*a[0][1]*a[0][2]));
      a[1][1]=180./PI*acos((s13-s23)/(2.*a[0][0]*a[0][2]));
      a[1][2]=180./PI*acos((s11-s22)/(4.*a[0][0]*a[0][1]));
      break;
   case ORTHOROMBIC_A:
      /* a[0][0]=a[0][0]; */
      a[0][1]=sqrt( s22 + s33 + 2.*s23 )/2.;
      a[0][2]=sqrt( s22 + s33 - 2.*s23 )/2.;
      a[1][0]=180./PI*acos((-s22+s33)/(4.*a[0][1]*a[0][2]));
      a[1][1]=180./PI*acos((-s12+s13)/(2.*a[0][0]*a[0][2]));
      a[1][2]=180./PI*acos(( s12+s13)/(2.*a[0][0]*a[0][1]));
      break;
   case ORTHOROMBIC_B:
      a[0][0]=sqrt( s11 + s33 - 2.*s13 )/2.;
      /* a[0][1]=a[0][1]; */
      a[0][2]=sqrt( s11 + s33 + 2.*s13 )/2.;
      a[1][0]=180./PI*acos((s12+s23)/(2.*a[0][1]*a[0][2]));
      a[1][1]=180./PI*acos((s11-s33)/(4.*a[0][0]*a[0][2]));     
      a[1][2]=180./PI*acos((s12-s23)/(2.*a[0][0]*a[0][1]));     
      break;
   case RHOMBOHEDRAL:
      /*  A_rhomb^2 = A_hex^2 /3 + C_hex^2 /9  */
      d=sqrt( a[0][0]*a[0][0]/3 + a[0][2]*a[0][2]/9 );
     /*  Cos[alpha]_rhomb = (-A_hex^2 /6 + C_hex^2 /9)/A_rhomb^2;  */
      d1=180./PI*acos((-a[0][0]*a[0][0]/6 + a[0][2]*a[0][2]/9)/(d*d));
      a[0][0]=a[0][1]=a[0][2]=d;
      a[1][0]=a[1][1]=a[1][2]=d1;
     /*  WIEN input supposes the primitive (rhombohedral) coordinates of atoms  */
      *lat=ORTHOROMBIC_P;
      break;
   case HEXAGONAL:
      a[1][0]=a[1][1]=90.; a[1][2]=120.;
      *lat=ORTHOROMBIC_P;
   }
     
  cl_in->nat=cl_in->nsrt=0;
  for(i=0; i < WIEN->cl.nsrt; i++) {
    if( fgets(s,256,fl)==NULL ) ERROR
    strncpy(st,s+5+3+4,10);           st[10]='\0';  cl_in->r[cl_in->nat][0]=atof(st);
    strncpy(st,s+5+3+4+10+3,10);      st[10]='\0';  cl_in->r[cl_in->nat][1]=atof(st);
    strncpy(st,s+5+3+4+10+3+10+3,10); st[10]='\0';  cl_in->r[cl_in->nat][2]=atof(st);
    asgn(WIEN->cl.r[cl_in->nat], cl_in->r[cl_in->nat]);
    WIEN->cl.srt[cl_in->nat] = i;
    cl_in->nat++;
    if( fgets(s,256,fl)==NULL ) ERROR /*  read MULT  ISPLIT */
    strncpy(st,s+15,2);  st[2]='\0';  MULT_W=atoi(st);
    strncpy(st,s+15+2+17,2);  st[2]='\0';  WIEN->isplit[cl_in->nsrt]=atoi(st);
    WIEN->cl.natsrt[i]=MULT_W;
    for(j=0; j<MULT_W-1; j++) {
      if( fgets(s,256,fl)==NULL ) ERROR
      strncpy(st,s+5+3+4,10);           st[10]='\0';  cl_in->r[cl_in->nat][0]=atof(st);
      strncpy(st,s+5+3+4+10+3,10);      st[10]='\0';  cl_in->r[cl_in->nat][1]=atof(st);
      strncpy(st,s+5+3+4+10+3+10+3,10); st[10]='\0';  cl_in->r[cl_in->nat][2]=atof(st);
      asgn(WIEN->cl.r[cl_in->nat], cl_in->r[cl_in->nat]);
      WIEN->cl.srt[cl_in->nat] = i;
      cl_in->nat++;
    }
    if( fgets(s,256,fl)==NULL ) ERROR /*  read Z  */ 
    n_to_zero(s);
    strncpy(WIEN->atom_name[cl_in->nsrt],s,10+5+5+5+10+5+10+5+5);
    WIEN->atom_name[cl_in->nsrt][10+5+5+5+10+5+10+5+5]='\0';
    strncpy(st,s+10+5+5+5+10+5+10+5,5); st[5]='\0';  WIEN->Z[cl_in->nsrt]=atof(st);
    /* generate sort name */
    /* check if atom name is presented */ 
    for(l=k=0; l<2; l++ )
    if( ! isspace((int)WIEN->atom_name[cl_in->nsrt][l]) ) {
        st[k++]=WIEN->atom_name[cl_in->nsrt][l]; /* atom name detected */
    }

    if ( k == 0 )
    fprintf(stderr, "warning: No literal name specified for atom %d.\n"
                    "         Using Z value as identifier.\n",i);
    for(l=0; l < MAX_DGT; l++ ) /* check 3-... characters to be digits*/
      if( isdigit((int)WIEN->atom_name[cl_in->nsrt][2+l]) ) 
          st[k+l]=WIEN->atom_name[cl_in->nsrt][2+l];
      else break;
    st[k+l]='\0';
    if( k== 0 )
       sprintf(name_srt[cl_in->nsrt],"z%s=%.2f",st,WIEN->Z[cl_in->nsrt]);
     else
       sprintf(name_srt[cl_in->nsrt],"%s",st);
     
    /*  find this sort  */
    for(k=0; k<cl_in->nsrt; k++)
    if( !strcmp(name_srt[cl_in->nsrt],name_srt[k]) ) { /*  this sort already exists  */
         for(l=nat_old; l < cl_in->nat; l++ ) cl_in->srt[l]=k;
         goto NEXT;
    }
    /*  new sort  */
    for(l=nat_old; l < cl_in->nat; l++ ) cl_in->srt[l]=cl_in->nsrt;
    cl_in->nsrt++;
    NEXT:
    nat_old=cl_in->nat;
    if( fgets(s,256,fl)==NULL ) ERROR /*  read ROT. matrix  */
    if( fgets(s,256,fl)==NULL ) ERROR /*  read ROT. matrix  */
    if( fgets(s,256,fl)==NULL ) ERROR /*  read ROT. matrix  */
  } /*  for(i=0; i<WIEN->nato; i++) {  */
  
  WIEN->cl.nat = cl_in->nat;
  fclose(fl);
  return 0;
}

int read_data(char *fname, t_cell *cl_in, char name_srt[MAX_ATOMS][MAX_CHARS],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
              int *lat)
/* "lat" here means only the type of centring mode (P,F,I,C,A),
 * therefore "lat" will take only values: ORTHOROMBIC_P, ORTHOROMBIC_F, ORTHOROMBIC_I,
 * ORTHOROMBIC_A, ORTHOROMBIC_C */
{
   FILE *fl;
   double d,d1;
   char s[256],st[64];
   int i,j;

   if( (fl=fopen(fname,"rt")) == NULL ) {
      fprintf(stderr,"error: Can't open file: %s.\n",fname);
      return 1;
   }

   *lat=ORTHOROMBIC_P;
/*  read type of cell  */
   do
     if( fgets(s,256,fl)==NULL ) ERROR
   while( stos(s,st)==0 );
/*  read translation vectors  */
   do
     if( fgets(s,256,fl)==NULL ) ERROR
   while( stod(s,&a[0][0],6)==0 );
   if( st[0]=='i' || st[0]=='I') {
      *lat=ORTHOROMBIC_I;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL || fabs(a[1][2]-90.)>TOL)
         fprintf(stderr,"warning: for I-centred type angles have to be: "
                "alpha=beta=gamma=90.\nThis relation will be assumed.\n" );
      d=a[0][0]*a[0][0]+a[0][1]*a[0][1]+a[0][2]*a[0][2];
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1]; a[0][2]*=a[0][2];
      a[1][0]=180./PI*acos(( a[0][0]-a[0][1]-a[0][2])/d);
      a[1][1]=180./PI*acos((-a[0][0]+a[0][1]-a[0][2])/d);
      a[1][2]=180./PI*acos((-a[0][0]-a[0][1]+a[0][2])/d);
      a[0][0]=a[0][1]=a[0][2]=sqrt(d)/2;
   } else
   if( st[0]=='f'|| st[0]=='F') {
      *lat=ORTHOROMBIC_F;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL || fabs(a[1][2]-90.)>TOL)
         fprintf(stderr,"warning: for F-centred type angles have to be: "
                "alpha=beta=gamma=90.\nThis relation will be assumed.\n" );
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1]; a[0][2]*=a[0][2];
      a[1][0]=180./PI*acos( a[0][0]/sqrt((a[0][0]+a[0][1])*(a[0][0]+a[0][2])) );
      a[1][1]=180./PI*acos( a[0][1]/sqrt((a[0][1]+a[0][0])*(a[0][1]+a[0][2])) );
      a[1][2]=180./PI*acos( a[0][2]/sqrt((a[0][2]+a[0][0])*(a[0][2]+a[0][1])) );
      d      =sqrt(a[0][1]+a[0][2])/2;
      d1     =sqrt(a[0][0]+a[0][2])/2;
      a[0][2]=sqrt(a[0][0]+a[0][1])/2;
      a[0][0]=d; a[0][1]=d1;
   } else
   if( st[0]=='c'|| st[0]=='C') { /*  C-centred: only orthorombic, cubic, tetragonal  */
      *lat=ORTHOROMBIC_C;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL || fabs(a[1][2]-90.)>TOL)
         fprintf(stderr,"warning: for C-centred type angles have to be: "
                "alpha=beta=gamma=90.\nThis relation will be assumed.\n" );
      a[1][0]=a[1][1]=90.;
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1];
      a[1][2]=180./PI*acos((a[0][0]-a[0][1])/(a[0][0]+a[0][1]));
      d=sqrt(a[0][0]+a[0][1])/2;
      a[0][0]=a[0][1]=d;
   } else
   if( st[0]=='a'|| st[0]=='A') { /*  A-centred: orthorombic and monoclinic  */
      *lat=ORTHOROMBIC_A;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL)
          fprintf(stderr,"warning: for A-centred type angles have to be: "
                  "alpha=beta=90.\nThis relation will be assumed.\n");
      a[1][0]=180./PI*acos((-a[0][1]*a[0][1]+a[0][2]*a[0][2])/
                           ( a[0][1]*a[0][1]+a[0][2]*a[0][2]));
      a[1][1]=180./PI*acos(-a[0][1]*cos(PI/180.*a[1][2])/
                           sqrt( (a[0][1]*a[0][1]+a[0][2]*a[0][2]) ) );
      a[1][2]=180./PI*acos( a[0][1]*cos(PI/180.*a[1][2])/
                           sqrt( (a[0][1]*a[0][1]+a[0][2]*a[0][2]) ) );
      d=sqrt(a[0][1]*a[0][1] + a[0][2]*a[0][2])/2;
      a[0][1]=a[0][2]=d;
   }
   
/*  read number of atoms  */
    do
     if( fgets(s,256,fl)==NULL ) ERROR
    while( stoi(s,&cl_in->nat,1)==0 );
    if( cl_in->nat > MAX_ATOMS ) {
       fprintf(stderr,"error: Parameter MAX_ATOMS=%d is too small!"
               "Increase it (type_sg.h file).\n", MAX_ATOMS);
       fclose(fl);
       return 1;
    }
/*  read coordinates and types  */
    cl_in->nsrt=0;
    for(i=0; i < cl_in->nat; i++)
    {
       do
         if( fgets(s,256,fl)==NULL ) ERROR      
      while( stod(s,cl_in->r[i],3)==0 );

       do
         if( fgets(s,256,fl)==NULL ) ERROR      
       while( stos(s,st)==0 );
       for(j=0; j <= cl_in->nsrt; j++){
          if(j == cl_in->nsrt) {
             if( cl_in->nsrt == MAX_ATOMS ) {
                fprintf(stderr,"error: Parameter MAX_ATOMS=%d is too small!"
                               "Increase it (type_sg.h file).\n", MAX_ATOMS);
                fclose(fl);
                return 1;
             }
             strcpy(name_srt[j],st);
             cl_in->nsrt++;
             cl_in->srt[i]=j;
             break;
          } else
          if( !strcmp(st,name_srt[j]) ) { /*  this sort already exists  */
             cl_in->srt[i]=j;
             break;
          }
       }
    }
    fclose(fl);
    return 0;
}

int read_data_cell(int ibrav, double a[2][3], /*  a,b,c, alpha,beta,gamma  */
              int *lat)
/* "lat" here means only the type of centring mode (P,F,I,C,A),
 * therefore "lat" will take only values: ORTHOROMBIC_P, ORTHOROMBIC_F, ORTHOROMBIC_I,
 * ORTHOROMBIC_A, ORTHOROMBIC_C */
{
   double d,d1;
   char s[256],st[64];
   int  i,j;

/*  read type of cell  */
/*  read translation vectors  */
   if(ibrav == 1) {
      *lat=CUBIC_P;
   } else
   if(ibrav == 2 || ibrav == 6 || ibrav == 9 ) {
	   if (ibrav == 2) *lat=CUBIC_I;
	   if (ibrav == 6) *lat=TETRAGONAL_I;
	   if (ibrav == 9) *lat=ORTHOROMBIC_I;
	   *lat=ORTHOROMBIC_I;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL || fabs(a[1][2]-90.)>TOL)
         fprintf(stderr,"warning: for I-centred type angles have to be: "
                "alpha=beta=gamma=90.\nThis relation will be assumed.\n" );
      d=a[0][0]*a[0][0]+a[0][1]*a[0][1]+a[0][2]*a[0][2];
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1]; a[0][2]*=a[0][2];
      a[1][0]=180./PI*acos(( a[0][0]-a[0][1]-a[0][2])/d);
      a[1][1]=180./PI*acos((-a[0][0]+a[0][1]-a[0][2])/d);
      a[1][2]=180./PI*acos((-a[0][0]-a[0][1]+a[0][2])/d);
      a[0][0]=a[0][1]=a[0][2]=sqrt(d)/2;
   } else
   if(ibrav == 3 || ibrav == 10 ) {
	   if (ibrav == 3) *lat=CUBIC_F;
	   if (ibrav == 10) *lat=ORTHOROMBIC_F;
	   *lat=ORTHOROMBIC_F;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL || fabs(a[1][2]-90.)>TOL)
         fprintf(stderr,"warning: for F-centred type angles have to be: "
                "alpha=beta=gamma=90.\nThis relation will be assumed.\n" );
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1]; a[0][2]*=a[0][2];
      a[1][0]=180./PI*acos( a[0][0]/sqrt((a[0][0]+a[0][1])*(a[0][0]+a[0][2])) );
      a[1][1]=180./PI*acos( a[0][1]/sqrt((a[0][1]+a[0][0])*(a[0][1]+a[0][2])) );
      a[1][2]=180./PI*acos( a[0][2]/sqrt((a[0][2]+a[0][0])*(a[0][2]+a[0][1])) );
      d      =sqrt(a[0][1]+a[0][2])/2;
      d1     =sqrt(a[0][0]+a[0][2])/2;
      a[0][2]=sqrt(a[0][0]+a[0][1])/2;
      a[0][0]=d; a[0][1]=d1;
   } else
/*   if(ibrav == ) {   C-centred: only orthorombic, cubic, tetragonal 
      *lat=ORTHOROMBIC_C;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL || fabs(a[1][2]-90.)>TOL)
         fprintf(stderr,"warning: for C-centred type angles have to be: "
                "alpha=beta=gamma=90.\nThis relation will be assumed.\n" );
      a[1][0]=a[1][1]=90.;
      a[0][0]*=a[0][0]; a[0][1]*=a[0][1];
      a[1][2]=180./PI*acos((a[0][0]-a[0][1])/(a[0][0]+a[0][1]));
      d=sqrt(a[0][0]+a[0][1])/2;
      a[0][0]=a[0][1]=d; 
   } else */
   if(ibrav == 11 || ibrav == 13 ) { /*  A-centred: orthorombic and monoclinic  */
      *lat=ORTHOROMBIC_A;
      if( fabs(a[1][0]-90.)>TOL || fabs(a[1][1]-90.)>TOL)
          fprintf(stderr,"warning: for A-centred type angles have to be: "
                  "alpha=beta=90.\nThis relation will be assumed.\n");
      a[1][0]=180./PI*acos((-a[0][1]*a[0][1]+a[0][2]*a[0][2])/
                           ( a[0][1]*a[0][1]+a[0][2]*a[0][2]));
      a[1][1]=180./PI*acos(-a[0][1]*cos(PI/180.*a[1][2])/
                           sqrt( (a[0][1]*a[0][1]+a[0][2]*a[0][2]) ) );
      a[1][2]=180./PI*acos( a[0][1]*cos(PI/180.*a[1][2])/
                           sqrt( (a[0][1]*a[0][1]+a[0][2]*a[0][2]) ) );
      d=sqrt(a[0][1]*a[0][1] + a[0][2]*a[0][2])/2;
      a[0][1]=a[0][2]=d;
   }else
   if(ibrav == 4){ /* Hexagonal Cell*/
      a[1][0]=a[1][1]=90.; a[1][2]=120.;
      *lat=ORTHOROMBIC_P;
/*      *lat=HEXAGONAL;*/
   }else
   if(ibrav == 7){ /* Rhombohedral Cell*/
      /*  A_rhomb^2 = A_hex^2 /3 + C_hex^2 /9  */
      d=sqrt( a[0][0]*a[0][0]/3 + a[0][2]*a[0][2]/9 );
     /*  Cos[alpha]_rhomb = (-A_hex^2 /6 + C_hex^2 /9)/A_rhomb^2;  */
      d1=180./PI*acos((-a[0][0]*a[0][0]/6 + a[0][2]*a[0][2]/9)/(d*d));
      a[0][0]=a[0][1]=a[0][2]=d;
      a[1][0]=a[1][1]=a[1][2]=d1;
     /*  WIEN input supposes the primitive (rhombohedral) coordinates of atoms  */
      *lat=ORTHOROMBIC_P;
   }else
   if(ibrav == 8){
      *lat=ORTHOROMBIC_P;
   }
    return 0;
}



/* check correspondence of input WIEN structure with space group found  */
#define ERROR_lat { \
 fprintf(fl,"warning: !!! Bravais lattice has changed.\n\n");\
 return 1; \
}

int chk_WIEN_struct( FILE *fl, int lat, t_cell *cl, t_WIEN *WIEN,
                     int npgrp, int nsgrp, int nop, double op[][4][3] )
{
   
   int i,j,atr,atl,isrt,ir,iop,npass;
   int icase_shft=0; /* shift index for T_mon[][3][3] */
   double Telpr[3][3],Tprel[3][3],r[3],r2[3];
   double Tab_orth[3][3]={   
               { 0.0,  -1.0,  0.0 },
               { 1.0,   0.0,  0.0 },
               { 0.0,   0.0,  1.0 }
   };
   double T_mon[5][3][3]={ 
      /* for 7,13,14 MONOCLINIC_P space groups, check cell choices 1,2,3 */
              {{ 1.0,   0.0,  0.0 },
               { 0.0,   1.0,  0.0 },
               { 0.0,   0.0,  1.0 }},

              {{ 0.0,  -1.0,  0.0 }, /* cell choice 2 over cell choice 1 */
               { 1.0,  -1.0,  0.0 },
               { 0.0,   0.0,  1.0 }},

              {{-1.0,   1.0,  0.0 }, /* cell choice 3 over cell choice 1 */
               {-1.0,   0.0,  0.0 },
               { 0.0,   0.0,  1.0 }},

      /* to check unique axis c and b */
          /* unique c over unique c */
              {{ 1.0,   0.0,  0.0 },
               { 0.0,   1.0,  0.0 },
               { 0.0,   0.0,  1.0 }},
          /* monoclinic unique axis b over unique axis c */      
              {{ 0.0,   0.0,  1.0 },
               { 1.0,   0.0,  0.0 },
               { 0.0,   1.0,  0.0 }}
   };
   
   npass=0; /* 1 pass totally */
   if( (npgrp == NCs  && nsgrp == 1) || 
       (npgrp == NC2h && nsgrp == 3) ||
       (npgrp == NC2h && nsgrp == 4) 
     ) npass = 2; /* 3 passes totally */
   
   get_Telpr(lat,Telpr);  inv_matr(Telpr,Tprel);
   /* ORTHOROMBIC_B -> ORTHOROMBIC_A */
   if( WIEN->lat == ORTHOROMBIC_B && lat == ORTHOROMBIC_A ) {
       /* WIEN->lat != MONOCLINIC_B */
       transform_coor( &WIEN->cl,Tab_orth );
       for(i=0; i < WIEN->cl.nat; i++) reduce(WIEN->cl.r[i]);
       WIEN->lat = ORTHOROMBIC_A;
   }
   
   /* unique axis b -> unique axis c */
   if( WIEN->lat == ORTHOROMBIC_P && lat == MONOCLINIC_P ) {
     if( eq(WIEN->a[1][0],90.) && !eq(WIEN->a[1][1],90.) &&
	 eq(WIEN->a[1][2],90.) ) {
       /* unique axis b case */
       transform_coor( &WIEN->cl,T_mon[4] );
       for(i=0; i < WIEN->cl.nat; i++) reduce(WIEN->cl.r[i]);
     } else
     if( eq(WIEN->a[1][0],90.) && eq(WIEN->a[1][1],90.) &&
	 eq(WIEN->a[1][2],90.) ) {
       /* unique axis b or c, check both cases */
	icase_shft=3;
	npass=1; /* 2 passes have to be performed */
     }
   }
   
   if( WIEN->cl.nsrt != cl->nsrt ) {
      fprintf(fl,"warning: !!! Number of inequivalent atoms has changed.\n"
                 "         !!! Old value= %d     New value= %d\n\n",
                 WIEN->cl.nsrt, cl->nsrt );
   }
   
   /* check Bravais lattice WIEN->lat = centering mode */
   switch( lat ) {
     
    case CUBIC_P:
    case CUBIC_I:
    case CUBIC_F:
      if( !( eq(WIEN->a[0][0],WIEN->a[0][1]) && /* a=b */
             eq(WIEN->a[0][0],WIEN->a[0][2]) && /* a=c */
             eq(WIEN->a[1][0],90.) &&           /* alpha=90 */
             eq(WIEN->a[1][1],90.) &&           /* beta=90  */
             eq(WIEN->a[1][2],90.) )            /* gamma=90 */
      ) ERROR_lat
      /* check B,F flags */
      if( (lat == CUBIC_P && WIEN->lat != ORTHOROMBIC_P) || 
          (lat == CUBIC_I && WIEN->lat != ORTHOROMBIC_I) || 
          (lat == CUBIC_F && WIEN->lat != ORTHOROMBIC_F) ) ERROR_lat
      break;
   
    case HEXAGONAL:
    case RHOMBOHEDRAL:
      if( !( eq(WIEN->a[0][0],WIEN->a[0][1]) && /* a=b */
             eq(WIEN->a[1][0],90.) &&           /* alpha=90 */
             eq(WIEN->a[1][1],90.) &&           /* beta=90  */
             eq(WIEN->a[1][2],120.) )           /* gamma=120 */
      ) ERROR_lat
      /* check H,R flags */
      if( (lat == HEXAGONAL && WIEN->lat != HEXAGONAL) || 
          (lat == RHOMBOHEDRAL && WIEN->lat != RHOMBOHEDRAL) ) ERROR_lat     
      break;

    case TETRAGONAL_P:
    case TETRAGONAL_I:
      if( !( eq(WIEN->a[0][0],WIEN->a[0][1]) && /* a=b */
             eq(WIEN->a[1][0],90.) &&           /* alpha=90 */
             eq(WIEN->a[1][1],90.) &&           /* beta=90  */
             eq(WIEN->a[1][2],90.) )            /* gamma=90 */
      ) ERROR_lat
      /* check B flag */
      if( (lat == TETRAGONAL_P && WIEN->lat != ORTHOROMBIC_P) ||
          (lat == TETRAGONAL_I && WIEN->lat != ORTHOROMBIC_I) ) ERROR_lat
      break;

    case ORTHOROMBIC_P:
    case ORTHOROMBIC_A:
    case ORTHOROMBIC_C:
    case ORTHOROMBIC_I:
    case ORTHOROMBIC_F:
      if( !( eq(WIEN->a[1][0],90.) &&          /* alpha=90 */
             eq(WIEN->a[1][1],90.) &&          /* beta=90  */
             eq(WIEN->a[1][2],90.) )           /* gamma=90 */
      ) ERROR_lat
      if( (lat == ORTHOROMBIC_P && WIEN->lat != ORTHOROMBIC_P) ||
          (lat == ORTHOROMBIC_A && WIEN->lat != ORTHOROMBIC_A) ||
          (lat == ORTHOROMBIC_C && WIEN->lat != ORTHOROMBIC_C) ||        
          (lat == ORTHOROMBIC_I && WIEN->lat != ORTHOROMBIC_I) || 
          (lat == ORTHOROMBIC_F && WIEN->lat != ORTHOROMBIC_F) ) ERROR_lat
      break;
   
    case MONOCLINIC_B:
      if( !( eq(WIEN->a[1][0],90.) &&           /* alpha=90 */
             eq(WIEN->a[1][1],90.) )            /* beta=90  */
      ) ERROR_lat
      /* check centering mode */
      if( lat == MONOCLINIC_B && WIEN->lat != ORTHOROMBIC_B ) ERROR_lat
      break;

    case MONOCLINIC_P:
      /* check unique axis b,c */
      if( !(eq(WIEN->a[1][0],90.) && eq(WIEN->a[1][2],90.)) &&
          !(eq(WIEN->a[1][0],90.) && eq(WIEN->a[1][1],90.)) ) ERROR_lat
      /* check centering mode */
      if( lat == MONOCLINIC_P && WIEN->lat != ORTHOROMBIC_P) ERROR_lat
      break;
      
    case TRICLINIC:
      if( WIEN->lat != ORTHOROMBIC_P) ERROR_lat
      break;
   } /* end switch */

   if( WIEN->cl.nat != cl->nat ) {
      fprintf(fl,"warning: !!! Unit cell has been reduced.\n\n");
      return 1;
   }
   
   /* apply symmetry operations to input cell */
   for( ; npass >=0; npass--) {
     for(isrt=atl=0; isrt < WIEN->cl.nsrt; isrt++ ) {
        atr=atl+WIEN->cl.natsrt[isrt]-1;
        for(i=atl; i <= atr; i=ir) {
           ir=i+1;
           for(iop=0; iop < nop; iop++ ) {
              mul_mv(r,T_mon[icase_shft+npass],WIEN->cl.r[i]);
              mul_mv(r,op[iop],r);
              add_vv(r,r,op[iop][3]);
              mul_mv(r,Tprel,r); reduce(r);
              /* search for r position at whole range */
              for(j=atl; j <= atr; j++) {
                 mul_mv(r2,T_mon[icase_shft+npass],WIEN->cl.r[j]);
                 mul_mv(r2,Tprel,r2);
                 reduce(r2);
                 if( eq(ddvec(r,r2),0) ) { /* found */
                    if( j >= ir ) {
                       if( ir > atr ) {
                          if( npass == 0 ) {
                             fprintf(fl,"warning: !!! Wrong multiplicity for atom %d\n"
                                     "         !!! Struct file is not consistent with "
                                     "space group found\n\n",
                                     isrt+1);
                             return 1;
                          }
                          goto NEXT_PASS;
                       }
                       swap(WIEN->cl.r[j], WIEN->cl.r[ir]);                    
                       ir++;
                    }
                    goto FOUND;
                 }
              }
              if( npass == 0 ) {              
                 fprintf(fl,"warning: !!! Struct file is not "
                         "consistent with space group found.\n\n");
                 return 1;
              }
              goto NEXT_PASS;
              FOUND:;
           } /* iop */
           if( ir != (atr+1) ) {
              if( npass == 0 ) {
                 fprintf(fl,"warning: !!! Wrong multiplicity for atom %d\n"
                         "         !!! Struct file is not consistent with "
                         "space group found.\n\n",isrt+1);
                 return 1;
              }
              goto NEXT_PASS;
           }
           atl+=WIEN->cl.natsrt[isrt];
        } /* atl <= i <= atr */
     }
     return 0; /* this pass was successful */
     NEXT_PASS:;
   }
   
   return 1;
} /* end of chk_WIEN_struct() */
#undef ERROR_lat
   
/* cl_new -- cell with full subdivision into nonequivalent atoms */
/* cl     -- cell with sorts given in the input file */
/* swithes:  if is_wien_input=1 -> is_prim=0, if lat==RHOMBOHEDRAL is_prim=1 */
/*           if is_wien_input=1 -> is_noeq=1 for centered lattices  */
int write_res(char *fname, t_cell *cl, char name_srt[MAX_ATOMS][MAX_CHARS],
              t_cell *cl_new, char name_srt_new[MAX_ATOMS][MAX_CHARS],
              t_switch *sw, int fix_org,
              int nsh, double Rsh[][3],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
              int lat, char *lat_name[], char * sgrp_name,
              int npgrp, int nsgrp, int nop, double sym_op[][4][3],
              char name_pgrp[][32],
              double Tin[3][3], t_WIEN *WIEN)
{
  FILE *fl;

  double r[3],Tel[3][3],Telpr[3][3],Tprel[3][3],
         Tr[4][3],pgrp_at[48][4][3], Tba_monoc[3][3];
  double op[48][4][3]; /* to find atom point group if is_prim==1 */
  double Tab_monoc[3][3]={   
               { 0.0,  -1.0,  0.0 },
               { 1.0,  -1.0,  0.0 },
               { 0.0,   0.0,  1.0 },
  };
   
  int i, j, npgrp_at, at, isrt, itr, ntr, lat_pgrp, lat_supp;

  if( fname == NULL) fl=stdout;  else
  if( (fl=fopen(fname,"at")) == NULL ) {
     fprintf(stderr,"error: Can not open file %s for writing.\n", fname);
     return(1);
  }
   
  inv_matr(Tab_monoc,Tba_monoc); 
   
  if( sw->is_prim || sw->is_noeq ) {
      ntr=1;  Tr[0][0]=Tr[0][1]=Tr[0][2]=0.;
  }  else get_centering_translations(lat,&ntr,Tr);

  asgn_n(op[0][0], sym_op[0][0], nop*4*3);
  /* MONOCLINIC_B -> MONOCLINIC_A to find atom point groups */
  if( sw->is_wien_input && lat==MONOCLINIC_B) {
      rotate_group_transl(nop,op,Tba_monoc);
  }
   
  if( sw->is_prim ) { /* transf. space group to conv. basis */
      if( sw->is_wien_input && lat == MONOCLINIC_B ) 
           get_Telpr(MONOCLINIC_A,Telpr);
      else get_Telpr(lat,Telpr);
     
      inv_matr(Telpr,Tprel);
      transform_group_opr(nop,op,Tprel);
      rotate_group_transl(nop,op,Tprel);
  }
     
   /* check if output file for WIEN was requested */
  if( sw->is_wien_input ) {
    if( chk_WIEN_struct(fl,lat,cl_new,WIEN,npgrp,nsgrp,nop,sym_op) )
    fprintf(fl,"----------------------------------------------------------------------\n\n");
  }
   
  if( sw->is_wien_input && lat == RHOMBOHEDRAL) {
     fprintf(fl,"NOTE: atom positions and space group operations\n"
                "      are given in rhombohedral primitive basis.\n"
                "      a,b,c, alpha,beta, gamma - parameters of hexagonal cell.\n\n");
      fprintf(fl,"Bravais lattice: %s\n\n",lat_name[lat]);     
  } else {
     if( sw->is_prim ) {
          fprintf(fl,"NOTE: atom positions and space group operations\n"
                     "      are given in primitive basis.\n\n");
          fprintf(fl,"Bravais lattice : %s\n\nParameters of primitive cell:\n",
                     lat_name[lat]);
     }
     else {
       if( lat == RHOMBOHEDRAL )
           fprintf(fl,"Bravais lattice: %s [hexagonal setting\n\n",lat_name[lat]);
       else
         fprintf(fl,"Bravais lattice: %s\n\n",lat_name[lat]);
     }
  }
   
  fprintf(fl,"     a             b            c\n");
  fprintf(fl,"% .8f   % .8f  % .8f\n",a[0][0],a[0][1],a[0][2]);
  fprintf(fl,"     alpha         beta         gamma\n");
  fprintf(fl,"% .8f   % .8f  % .8f\n\n",a[1][0],a[1][1],a[1][2]);

  if(sw->is_wien_input && WIEN->lat == RHOMBOHEDRAL && lat == RHOMBOHEDRAL)
    fprintf(fl,"\n===== Decomposition of new RHOMBOHEDRAL basis vectors "
               "\n===== over input RHOMBOHEDRAL (not HEXAGONAL!) basis\n");
  else
    fprintf(fl,"\n===== Decomposition of new basis vectors over input basis =====\n");
   
  fprintf(fl,"% .6f  % .6f % .6f  <--- 1\n",Tin[0][0],Tin[1][0],Tin[2][0]);
  fprintf(fl,"% .6f  % .6f % .6f  <--- 2\n",Tin[0][1],Tin[1][1],Tin[2][1]);
  fprintf(fl,"% .6f  % .6f % .6f  <--- 3\n",Tin[0][2],Tin[1][2],Tin[2][2]);
   
  if( sw->is_noeq ) {
    fprintf(fl,"\n==== Number of atoms in cell "
               "(only atoms of primitive cell): %d \n",cl->nat*ntr);
    fprintf(fl,"==== Atom positions (only atoms of primitive cell):\n");
  } else {
    fprintf(fl,"\n==== Number of atoms in cell: %d\n",cl->nat*ntr);
    fprintf(fl,"==== Atom positions:\n");
  }
  for(isrt=i=0; isrt < cl->nsrt; isrt++) {
     fprintf(fl,"\n");
     for(itr=0; itr < ntr; itr++)
     for(at=0; at < cl->natsrt[isrt]; at++) {
        add_vv(r,cl->r[i+at],Tr[itr]); reduce(r);
        fprintf(fl,"% .8f  % .8f % .8f\n",r[0],r[1],r[2]);
        fprintf(fl," %s\n",name_srt[cl->srt[i]]);
     }
     i+=cl->natsrt[isrt];
  }
  fprintf(fl,"\n==== Number of nonequivalent sorts: %d\n", cl_new->nsrt);
  if( ! sw->is_prim)
  fprintf(fl,"==== Nonequivalent atoms, point group for each sort: ====\n");
  else
  fprintf(fl,"\n==== Nonequivalent atoms for each sort: ====\n");
  for(isrt=i=0; isrt < cl_new->nsrt; isrt++) {
     fprintf(fl, "\nSort number: %d\n",isrt+1);
     
 /* some job for find_atom_pgroup() */
     lat_supp=lat;
     asgn(r,cl_new->r[i]);
     if( lat == MONOCLINIC_B ) { /* wien_input requested */
        mul_mv(r,Tab_monoc,cl_new->r[i]); reduce(r);
        lat_supp=MONOCLINIC_A;
        /* mul_mm(Tel,Tba_monoc,Tel); */
     }
     if( sw->is_prim ) {
        mul_mv(r,Telpr,r); reduce(r);
     }
     
     find_atom_pgroup(lat_supp,nop,npgrp,op,r,&lat_pgrp,&npgrp_at,pgrp_at,Tel);
     
     round_to_zero(Tel[0],3*3);
     fprintf(fl,"  Names of point group: %s\n", name_pgrp[npgrp_at]);
     fprintf(fl,"  New basis vectors for this point group:\n");
     fprintf(fl,"   % .4f  % .4f % .4f  <--- 1\n",Tel[0][0],Tel[1][0],Tel[2][0]);
     fprintf(fl,"   % .4f  % .4f % .4f  <--- 2\n",Tel[0][1],Tel[1][1],Tel[2][1]);
     fprintf(fl,"   % .4f  % .4f % .4f  <--- 3\n",Tel[0][2],Tel[1][2],Tel[2][2]);
     fprintf(fl,"\n");
     
     fprintf(fl,"  Atom positions: %d\n", cl_new->natsrt[isrt]*ntr);
     for(itr=0; itr < ntr; itr++)
     for(at=0; at < cl_new->natsrt[isrt]; at++) {
        add_vv(r,cl_new->r[i+at],Tr[itr]); reduce(r);
        fprintf(fl,"  % .8f  % .8f % .8f\n",r[0],r[1],r[2]);
        fprintf(fl,"   %s\n",name_srt_new[cl_new->srt[i]]);
     }
     i+=cl_new->natsrt[isrt];
  }

  fprintf(fl,"\n=========================================================================");
  fprintf(fl,"\n========= Bravais Lattice ====== point group ====== Space group =========");
  fprintf(fl,"\n=========================================================================\n\n");
   
  if( sw->is_wien_input && lat == RHOMBOHEDRAL) {
     fprintf(fl,"NOTE: atom positions and space group operations\n"
                "      are given in rhombohedral primitive basis.\n"
                "      a,b,c, alpha,beta, gamma - parameters of hexagonal cell.\n\n");
      fprintf(fl,"Bravais lattice: %s\n\n",lat_name[lat]);     
  } else {
     if( sw->is_prim ) {
          fprintf(fl,"NOTE: atom positions and space group operations\n"
                     "      are given in primitive basis.\n\n");
          fprintf(fl,"Bravais lattice : %s\n\nParameters of primitive cell:\n",
                     lat_name[lat]);
     }
     else {
       if( lat == RHOMBOHEDRAL )
           fprintf(fl,"Bravais lattice: %s [hexagonal setting\n\n",lat_name[lat]);
       else
         fprintf(fl,"Bravais lattice: %s\n\n",lat_name[lat]);
     }
  }
   
  fprintf(fl,"Names of point group [Short-Full-Schoenflies]: %s\n",
          name_pgrp[npgrp]);

  if( lat==RHOMBOHEDRAL && ! sw->is_prim )
     fprintf(fl,"\nNumber and name of space group: %s [h axes]\n",sgrp_name);
   else
     fprintf(fl,"\nNumber and name of space group: %s\n",sgrp_name);
  fprintf(fl,"\n=========================================================================\n");

  fprintf(fl,"\nNumber of symmetry operations: %d\n",nop);
  for(i=0; i < nop; i++) {
       fprintf(fl,"Operation: %d\n",i+1);
       for(j=0; j<3; j++)
           if( fabs(sym_op[i][3][j]-1./6) < 1.e-5 )
           fprintf(fl,"% .1f  % .1f  % .1f  1/6\n",
                   sym_op[i][j][0],sym_op[i][j][1],sym_op[i][j][2]);
           else
           if( fabs(sym_op[i][3][j]-1./3) < 1.e-5 )
           fprintf(fl,"% .1f  % .1f  % .1f  1/3\n",
                   sym_op[i][j][0],sym_op[i][j][1],sym_op[i][j][2]);
           else
           if( fabs(sym_op[i][3][j]-2./3) < 1.e-5 )
           fprintf(fl,"% .1f  % .1f  % .1f  2/3\n",
                   sym_op[i][j][0],sym_op[i][j][1],sym_op[i][j][2]);
           else
           if( fabs(sym_op[i][3][j]-5./6) < 1.e-5 )
           fprintf(fl,"% .1f  % .1f  % .1f  5/6\n",
                   sym_op[i][j][0],sym_op[i][j][1],sym_op[i][j][2]);
           else
           fprintf(fl,"% .1f  % .1f  % .1f % .3f\n",
                   sym_op[i][j][0],sym_op[i][j][1],sym_op[i][j][2],sym_op[i][3][j]);
       fprintf(fl,"\n");
   }

     fprintf(fl,"============================================================");
     if( fix_org == 001 ) {
        fprintf(fl,"\nNote that shift vectors for this space group are defined\n"
               "only up to the vector ");
        if( sw->is_prim ) {
           if( lat==CUBIC_P || lat==HEXAGONAL ||
               lat==TETRAGONAL_P || lat==ORTHOROMBIC_P || lat==MONOCLINIC_P)
               fprintf(fl,"{ 0, 0, Z }."); else
           if( lat==RHOMBOHEDRAL )
               fprintf(fl,"{ Z, Z, Z }."); else 
           if( lat==CUBIC_I || lat==TETRAGONAL_I || lat== ORTHOROMBIC_I)
               fprintf(fl,"{ Z, Z, 0 }."); else
           if( lat==CUBIC_F || lat==ORTHOROMBIC_F )
               fprintf(fl,"{ Z, Z, -Z }."); else
           if( lat==ORTHOROMBIC_C )
               fprintf(fl,"{ 0, 0, Z }."); else
           if( lat==ORTHOROMBIC_A || lat==MONOCLINIC_A )
               fprintf(fl,"{ 0, Z, Z }.");
        }
        else fprintf(fl,"{ 0, 0, Z }.");        
        fprintf(fl,"\nHere Z can take any value.");     
        
     } else
     if( fix_org == 110 ) {
        fprintf(fl,"\nNote that shift vectors for this space group are defined\n"
               "only up to the vector ");
        if( sw->is_prim ) {
           if( lat==CUBIC_P || lat==HEXAGONAL ||
               lat==TETRAGONAL_P || lat==ORTHOROMBIC_P || lat==MONOCLINIC_P)
               fprintf(fl,"{ X, Y, 0 }."); else
           if( lat==RHOMBOHEDRAL )
               fprintf(fl,"{ -Y, X, -X+Y }."); else     
           if( lat==CUBIC_I || lat==TETRAGONAL_I || lat== ORTHOROMBIC_I)
               fprintf(fl,"{ Y, X, X+Y }."); else
           if( lat==CUBIC_F || lat==ORTHOROMBIC_F )
               fprintf(fl,"{ -X+Y, X-Y, X+Y }."); else
           if( lat==ORTHOROMBIC_C )
               fprintf(fl,"{ X-Y, X+Y, 0 }."); else
           if( lat==ORTHOROMBIC_A || lat==MONOCLINIC_A )
               fprintf(fl,"{ X, Y, -Y }.");
        }
        else fprintf(fl,"{ X, Y, 0 }.");        
        fprintf(fl,"\nHere X and Y can take any values.");
     } else
     if( fix_org == 111 )
        fprintf(fl,"\nNote that shift vector for this space group is defined\n"
                    "only up to the vector { X, Y, Z }.\nHere X, Y and Z can take any values.");

     fprintf(fl,"\n\n===== Number of possible shift vectors: %d =====\n",nsh);
     fprintf(fl,"===== List of shift vectors:\n");
     for(i=0; i<nsh; i++)
        fprintf(fl,"% .4f  % .4f  % .4f\n",Rsh[i][0], Rsh[i][1], Rsh[i][2]);
   if( nsh > 1 ) {
      fprintf(fl,"\nList of shifted cells:\n");
      for(j=1; j<nsh; j++) {
         fprintf(fl,"\nCell #%d \n",j);
         fprintf(fl,"  Shift vector: % .4f  % .4f  % .4f\n",Rsh[j][0], Rsh[j][1], Rsh[j][2]);
         fprintf(fl,"  Nonequivalent atoms :\n");
         for(isrt=i=0; isrt < cl_new->nsrt; isrt++) {
            fprintf(fl,"\n");
            for(itr=0; itr < ntr; itr++)            
            for(at=0; at < cl_new->natsrt[isrt]; at++) {
               sub_vv(r,cl_new->r[i+at],Rsh[j]);
               add_vv(r,r,Tr[itr]);
               reduce(r);
               fprintf(fl,"% .8f  % .8f % .8f\n",r[0], r[1],r[2]);
               fprintf(fl," %s\n",name_srt_new[cl_new->srt[i]]);
            }
            i+=cl_new->natsrt[isrt];
         }
      }
   }
  fclose(fl);
  return 0;
}

/* obtain number of subsort with given Z value */
int get_Z_subsrt( double Z, t_cell *cl, t_WIEN *WIEN)
{
   static int first_call=1;
   int i;
   
   if( first_call ) {
       for(i=0; i < cl->nsrt; i++) WIEN->indZ[i]=0;
       first_call=0;
   }
   /* find which index should be used to refer to
    *  WIEN->indZ[] array, this index is the
    * first occurence of a given Z in WIEN->Z[] array*/
   for(i=0; i < cl->nsrt; i++)
     if( eq(Z,WIEN->Z[i]) ) return ++WIEN->indZ[i];

   fprintf(stderr,"error in get_Z_subsrt(): for Z=%f index not found",Z);
   return -1;
}
/* write out results to WIEN struct file */
/* cl_new -- cell with full subdivision into nonequivalent atoms */
/* cl     -- cell with sorts given in the input file */
int write_WIEN_struct(
    char *fname, t_cell *cl, char name_srt[MAX_ATOMS][MAX_CHARS],
    t_cell *cl_new, char name_srt_new[MAX_ATOMS][MAX_CHARS],
    int fix_org, int nsh, double Rsh[][3],
    double a[2][3], /*  a,b,c, alpha,beta,gamma  */
    int lat, char *lat_name[], char * sgrp_name,
    int npgrp, int nop, double sym_op[][4][3],
    char name_pgrp[][32],
    double Tin[3][3], t_WIEN *WIEN)
{
  FILE *fl;

  char *sptr, s[256];
  char dksgrp_name[100];
  int dki,dkj;
  int i, j, isrt, isub_srt;

  if( fname == NULL) fl=stdout;  else
  if( (fl=fopen(fname,"at")) == NULL ) {
     fprintf(stderr,"Can not open file: %s.\n", fname);
     return 1;
  }

  fprintf(fl,"%s\n", WIEN->title); /* 1 line */
  
/* 2 line
 * define centering mode and lattice type: P,F,B,CXY,CYZ,CXZ,R,H
 * WIEN treats both ORTHOROMBIC_A and ORTHOROMBIC_B,
 * but here only for ORTHOROMBIC_A output is given  */
  switch( lat ) {
    case CUBIC_P:
    case TETRAGONAL_P:
    case ORTHOROMBIC_P:
    case MONOCLINIC_P:
    case TRICLINIC:
                    sptr="P   ";  break;
    case CUBIC_F:
    case ORTHOROMBIC_F:
                    sptr="F   ";  break;
    case CUBIC_I:
    case TETRAGONAL_I:
    case ORTHOROMBIC_I:
                    sptr="B   ";  break;
    case ORTHOROMBIC_C:
                    sptr="CXY ";  break;
    case ORTHOROMBIC_A:
                    sptr="CYZ ";  break;
    case MONOCLINIC_B:
                    sptr="CXZ ";  break;
    case RHOMBOHEDRAL:
                    sptr="R   ";  break;   
    case HEXAGONAL:
                    sptr="H   ";  break;
    default:
                    fprintf(stderr,"Error in write_WIEN_struct(): "
                                   "not supported lattice type.\n");
                    return 1;
     
  }
  dki=dkj=0;
  dksgrp_name[dki++]=sgrp_name[dkj++]; /* first character */
  if(sgrp_name[dkj] != ' ') { /* second character */
    dksgrp_name[dki++]=sgrp_name[dkj++];
    if(sgrp_name[dkj] != ' ') /* third character */
      dksgrp_name[dki++]=sgrp_name[dkj++]; }
  dksgrp_name[dki++]=sgrp_name[dkj++]; /* first blank */
  for(;sgrp_name[dkj]!=')';dkj++)
    if(sgrp_name[dkj] != ' ' && sgrp_name[dkj] != '(' )
      dksgrp_name[dki++]=sgrp_name[dkj];
  dksgrp_name[dki]='\0';
  fprintf(fl,"%sLATTICE,NONEQUIV.ATOMS:%3d" " %s\n",
              sptr,cl_new->nsrt,dksgrp_name);
   
  /* 3 line  */
  fprintf(fl,"%s\n",WIEN->mode);
   
 /* 4 line  a,b,c, alpha, beta, gamma */
  sprintFwd(s,a[0][0],10,6); fprintf(fl,"%s",s);
  sprintFwd(s,a[0][1],10,6); fprintf(fl,"%s",s);
  sprintFwd(s,a[0][2],10,6); fprintf(fl,"%s",s);   
  sprintFwd(s,a[1][0],10,6); fprintf(fl,"%s",s);
  sprintFwd(s,a[1][1],10,6); fprintf(fl,"%s",s);
  sprintFwd(s,a[1][2],10,6); fprintf(fl,"%s",s);   fprintf(fl,"\n");   

  /* 5 line -->  loop over sorts */
  for(isrt=i=0; isrt < cl_new->nsrt; isrt++) {
     fprintf(fl,"ATOM%4d: X=",cl_new->srt[i]+1);
     sprintFwd(s,cl_new->r[i][0],10,8); fprintf(fl,"%s Y=",s);
     sprintFwd(s,cl_new->r[i][1],10,8); fprintf(fl,"%s Z=",s);
     sprintFwd(s,cl_new->r[i][2],10,8); fprintf(fl,"%s\n",s);
     i++;
     /* write MULT ISPLIT line */
     fprintf(fl,"          MULT=%2d          ISPLIT=%2d\n",
             cl_new->natsrt[isrt],WIEN->isplit[cl->srt[i-1]]);
     for(j=1; j < cl_new->natsrt[isrt]; j++) {
       fprintf(fl,"    %4d: X=",cl_new->srt[i]+1);
       sprintFwd(s,cl_new->r[i][0],10,8); fprintf(fl,"%s Y=",s);
       sprintFwd(s,cl_new->r[i][1],10,8); fprintf(fl,"%s Z=",s);
       sprintFwd(s,cl_new->r[i][2],10,8); fprintf(fl,"%s\n",s);
       i++;
     }
     /* get subsort number for this Z */
     if( (isub_srt = get_Z_subsrt(WIEN->Z[cl->srt[i-1]], cl, WIEN)) == -1 )
       return 1;
     /* print atom name */
     sprintf(s,"%d",isub_srt);
     for(j=0; j < MAX_DGT; j++) if( ! isdigit((int)s[j]) ) break;
     for(; j < MAX_DGT; j++) s[j]=' ';
     s[j]='\0';
     fprintf(fl,"%.2s%s%s\n", WIEN->atom_name[cl->srt[i-1]],s,
                              &(WIEN->atom_name[cl->srt[i-1]][2+MAX_DGT]));

     fprintf(fl,"LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000\n"
                "                     0.0000000 1.0000000 0.0000000\n"
                "                     0.0000000 0.0000000 1.0000000\n");
  } /* end loop over sorts */
  
  /* 11 line */
   fprintf(fl,"%4d      NUMBER OF SYMMETRY OPERATIONS\n",nop);
  
  /* 12-14 lines */
  for(i=0; i < nop; i++) {
     for(j=0; j<3; j++) {
       fprintf(fl,"%2d%2d%2d",
               (int)rint(sym_op[i][j][0]),(int)rint(sym_op[i][j][1]),
               (int)rint(sym_op[i][j][2]));
       sprintFwd(s,sym_op[i][3][j],11,8); fprintf(fl,"%s\n",s);
     }
     fprintf(fl,"%8d\n",i+1);
  }
  fclose(fl);
  return 0; 
}
