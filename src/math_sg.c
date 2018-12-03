#include <math.h>
#include "type_sg.h"

#ifndef PI
#define PI 3.14159265358979323844
#endif

double*  asgn(double a[3], double b[3])
{
   a[0]=b[0]; a[1]=b[1]; a[2]=b[2];
   return(a);
}

double*  asgn_n(double *a, double *b, int n)
{
   int i;
   for(i=0; i<n; i++ ) a[i]=b[i];
   return(a);
}

int is_eq(double *a, double *b, int n)
{
   int i;
   for(i=0; i<n; i++ ) if( fabs(a[i]-b[i]) > TOL ) return(0);
   return(1);
}

double*  mul_mv(double p[3], double a[3][3], double r[3])
{  /*  p=a*r  */
   int i,j;
   double r1[3];
   for(i=0;i<3;i++) {
      r1[i]=0.;
      for(j=0;j<3;j++) r1[i]+=a[i][j]*r[j];
   }
   p[0]=r1[0]; p[1]=r1[1]; p[2]=r1[2];
   return(p);
}

double*  mul_vm(double p[3], double a[3][3], double r[3])
{  /*  p=r*a  */
   int i,j;
   double r1[3];
   for(i=0;i<3;i++) {
      r1[i]=0.;
      for(j=0;j<3;j++) r1[i]+=a[j][i]*r[j];
   }
   p[0]=r1[0]; p[1]=r1[1]; p[2]=r1[2];
   return(p);
}

/*  c=a*b  */
void  mul_mm(double c[3][3], double a[3][3], double b[3][3])
{
   double d[3][3]; int i,j,k;
   for(i=0;i<3;i++)
   for(j=0;j<3;j++)
   for(d[i][j]=0.,k=0; k<3; k++) d[i][j]+=a[i][k]*b[k][j];

   for(i=0;i<3;i++)
   for(j=0;j<3;j++) c[i][j]=d[i][j];
}

void transp(double a[3][3],double at[3][3])
{
   int i,j;
   double b[3][3];
   for(i=0; i<3; i++)
   for(j=0; j<3; j++)
       b[i][j]=a[j][i];
   
   for(i=0; i<3; i++)
   for(j=0; j<3; j++)
       at[i][j]=b[i][j];
   
/*   double t;
   t=a[0][1]; a[0][1]=a[1][0]; a[1][0]=t;
   t=a[0][2]; a[0][2]=a[2][0]; a[2][0]=t;
   t=a[1][2]; a[1][2]=a[2][1]; a[2][1]=t;*/
   return;
}

int norm_vec(double r[3]) /*  set length of vec=1  */
{
  double d;
  d=sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
  if( d < TOL ) return(0);
  r[0]/=d; r[1]/=d; r[2]/=d;
  return(1);
}

double dvec(double r[3])
{
   return(sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

double ddvec(double r1[3], double r2[3])
{
  return( sqrt( (r1[0]-r2[0])*(r1[0]-r2[0]) + (r1[1]-r2[1])*(r1[1]-r2[1]) +
                (r1[2]-r2[2])*(r1[2]-r2[2])) );
}

double dvec_n(double r[], int n)
{
  int i;
  double s=0.;
  for(i=0; i<n; i++) s+=r[i]*r[i];
  return( sqrt(s));
}

double ddvec_n(double r1[], double r2[], int n)
{
  int i;
  double s=0.;
  for(i=0; i<n; i++) s+=(r1[i]-r2[i])*(r1[i]-r2[i]);
  return( sqrt(s));
}

double*  add_vv(double c[3], double a[3], double b[3])
{
   double t[3];
   t[0]=a[0]+b[0]; c[0]=t[0];
   t[1]=a[1]+b[1]; c[1]=t[1];
   t[2]=a[2]+b[2]; c[2]=t[2];
   return(c);
}

double*  sub_vv(double c[3], double a[3], double b[3])
{
   double t[3];
   t[0]=a[0]-b[0]; c[0]=t[0];
   t[1]=a[1]-b[1]; c[1]=t[1];
   t[2]=a[2]-b[2]; c[2]=t[2];
   return(c);
}

void sub_mm(double c[3][3], double a[3][3], double b[3][3])
{
   int i; double t[3][3];
   for(i=0; i<9; i++ ) t[0][i]=a[0][i]-b[0][i];
   for(i=0; i<9; i++ ) c[0][i]=t[0][i];
}

double* mul_cv(double b[3], double c, double a[3])
{ /*  b=c*a  */
   double t[3];
   t[0]=c*a[0]; b[0]=t[0];
   t[1]=c*a[1]; b[1]=t[1];
   t[2]=c*a[2]; b[2]=t[2];
   return(b);
}

double mul_vv(double r1[3], double r2[3])
{
   return(r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2]);
}

double*  reduce(double r[3])
{
     r[0]=r[0]-rint(r[0]);
     if( fabs(r[0]) < TOL ) r[0]=0.;
     else  
     if( r[0] < 0. ) r[0]+=1.;

     r[1]=r[1]-rint(r[1]);
     if( fabs(r[1]) < TOL ) r[1]=0.;
     else  
     if( r[1] < 0. ) r[1]+=1.;

     r[2]=r[2]-rint(r[2]);
     if( fabs(r[2]) < TOL ) r[2]=0.;
     else  
     if( r[2] < 0. ) r[2]+=1.;

     return(r);
}

double*  put_in_cell(double r[3], double T[3][3], double T_1[3][3])
{
   mul_mv(r,T_1,r);    /*  transform from cubic to lattice basis  */
   reduce(r);
   mul_mv(r,T,r);      /*  back transformation  */
   return(r);
}

double vol( double r0[3],double r1[3],double r2[3])
{
  return(fabs(
         r0[0]*r1[1]*r2[2] + r0[1]*r1[2]*r2[0] +
         r1[0]*r2[1]*r0[2] - r2[0]*r1[1]*r0[2] -
         r0[0]*r1[2]*r2[1] - r1[0]*r0[1]*r2[2]
              ));
}

double detr(double A[3][3]) /*  determinant of matr A[3x3]  */
{
  double det;
  det= A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] +
       A[1][0]*A[2][1]*A[0][2] -
       A[2][0]*A[1][1]*A[0][2] - A[2][1]*A[1][2]*A[0][0] -
       A[0][1]*A[1][0]*A[2][2];
  return( det );
}
  
int inv_matr(double A[3][3], double B[3][3]) /*  B=inverse matr. of A  */
{
  /* return 1 if det(A)=0 */
  int i;
  double det,C[3][3];

  det=detr(A);
  if( fabs(det) < 1e-28 ) return(1);
  
  C[0][0]=(A[1][1]*A[2][2]-A[2][1]*A[1][2])/det;
  C[1][1]=(A[0][0]*A[2][2]-A[2][0]*A[0][2])/det;
  C[2][2]=(A[0][0]*A[1][1]-A[1][0]*A[0][1])/det;

  C[0][1]=-(A[0][1]*A[2][2]-A[2][1]*A[0][2])/det;
  C[1][0]=-(A[1][0]*A[2][2]-A[2][0]*A[1][2])/det;
  C[0][2]= (A[0][1]*A[1][2]-A[1][1]*A[0][2])/det;
  C[2][0]= (A[1][0]*A[2][1]-A[2][0]*A[1][1])/det;
  C[1][2]=-(A[0][0]*A[1][2]-A[1][0]*A[0][2])/det;
  C[2][1]=-(A[0][0]*A[2][1]-A[2][0]*A[0][1])/det;

  for(i=0; i<9; i++) B[0][i]=C[0][i];
  return(0);
}

void swap( double a[3],double b[3])
{
   double t;
   t=a[0]; a[0]=b[0]; b[0]=t;
   t=a[1]; a[1]=b[1]; b[1]=t;
   t=a[2]; a[2]=b[2]; b[2]=t;
   return;
}

/*  vector product of two vectors c=[axb]  */
double*  mulv_vv( double c[3], double a[3], double b[3])
{
   double t[3];
   t[0]=a[1]*b[2]-a[2]*b[1];
   t[1]=a[2]*b[0]-a[0]*b[2];
   t[2]=a[0]*b[1]-a[1]*b[0];
   c[0]=t[0]; c[1]=t[1]; c[2]=t[2];
   return(c);
}

void right_basis(double T1[3], double T2[3], double T3[3])
{
   double a[3];
   mulv_vv(a,T1,T2);
   if( mul_vv(a,T3) < 0 ) swap(T1,T2);
   return;
}

double*  colmn(double r[3], double a[3][3], int n)
{
   r[0]=a[0][n]; r[1]=a[1][n]; r[2]=a[2][n];
   return(r);
}

void sorti(int *a, int n)
{
   int i,j,k;
   for(i=0; i<n; i++)
   for(j=i+1; j<n; j++)
     if( a[i] > a[j] ) {
        k=a[i]; a[i]=a[j]; a[j]=k;
     }
}

void round_to_zero( double *d, int n)
{
   int i;
   for(i=0; i<n; i++)  if( fabs(d[i]) < TOL ) d[i]=0.;
   return;
}

int find_vec_ivec(double p[3], double (*vecs)[3], int n)
{
   int i;
   double ip[3];
   
   ip[0]=-p[0]; ip[1]=-p[1]; ip[2]=-p[2];
   for(i=0; i<n; i++)
     if( ddvec(p,vecs[i]) < TOL || ddvec(ip,vecs[i]) < TOL ) return(i);
   return(-1);
}
   
  /*** reduce T1,T2 to minimal length ***/
int shortest2( double T1[3], double T2[3])
{
#define MAXITER 10000
   
   int iter,i1,fl_ch;
   double d,d1,d2,p[3];
   
   d1=T1[0]*T1[0]+T1[1]*T1[1]+T1[2]*T1[2];
   d2=T2[0]*T2[0]+T2[1]*T2[1]+T2[2]*T2[2];
   for( iter=0; iter<MAXITER; iter++) {
      fl_ch=0;
      for(i1=-1; i1<=1; i1+=2) {
         p[0]=T1[0]+i1*T2[0];
         p[1]=T1[1]+i1*T2[1];
         p[2]=T1[2]+i1*T2[2];
         d=p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
         if( d<d1 && fabs(d-d1) > TOL ) {
            T1[0]=p[0]; T1[1]=p[1]; T1[2]=p[2];
            d1=d; fl_ch=1;
         } else
         if( d<d2 && fabs(d-d2) > TOL ) {
            T2[0]=p[0]; T2[1]=p[1]; T2[2]=p[2];
            d2=d; fl_ch=1;
        }
      }
      if( !fl_ch ) return(0);
  }
/*    printf(" Error in shortest2(): shortest vectors not found!\n");  */
  return(1);
}
