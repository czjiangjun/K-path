#ifndef MATH_SG_INC
#define MATH_SG_INC


#define eq(x,y) (fabs((x)-(y)) < TOL )
#define lt(x,y) ( (x) < (y)  && fabs((x)-(y)) > TOL )
#define gt(x,y) ( (x) > (y)  && fabs((x)-(y)) > TOL )
#define le(x,y) ( (x) < (y)  || fabs((x)-(y)) < TOL )
#define ne(x,y) (fabs((x)-(y)) > TOL )

extern double*  asgn(double a[3], double b[3]);
extern double*  asgn_n(double *a, double *b, int n);
extern double*  mul_mv(double p[3], double a[3][3], double r[3]);
extern double*  mul_vm(double p[3], double a[3][3], double r[3]);
extern void     mul_mm(double c[3][3], double a[3][3], double b[3][3]);
extern void     transp(double a[3][3], double at[3][3]);
extern int      norm_vec(double r[3]); /*  set length of vec=1  */
extern double   dvec(double r[3]);
extern double   dvec_n(double r[],int n);
extern double   ddvec(double r1[3], double r2[3]);
extern double   ddvec_n(double r1[], double r2[], int n);
extern double*  reduce(double r[3]);
extern double*  put_in_cell(double r[3], double T[3][3], double T_1[3][3]);
extern double   mul_vv(double r1[3], double r2[3]);
extern double   vol( double r0[3],double r1[3],double r2[3]);
extern double   detr(double A[3][3]); /*  determinant of matr A[3x3]  */
extern int      inv_matr(double A[3][3], double B[3][3]); /*  B=inverse matr. of A  */
extern void     swap( double a[3],double b[3]);
extern double*  mulv_vv( double c[3], double a[3], double b[3]);
extern void     right_basis(double T1[3], double T2[3], double T3[3]);
extern double*  add_vv(double c[3], double a[3], double b[3]);
extern double*  sub_vv(double c[3], double a[3], double b[3]);
extern void     sub_mm(double c[3][3], double a[3][3], double b[3][3]);
extern double*  mul_cv(double b[3], double c, double a[3]);
extern double*  colmn(double r[3], double a[3][3], int n);
extern int      is_eq(double *a, double *b, int n);
extern void     sorti(int *a, int n);
extern void round_to_zero( double *d, int n);
extern int find_vec_ivec(double p[3], double (*vecs)[3], int n);
extern int shortest2( double T1[3], double T2[3]);

#endif
