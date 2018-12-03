/*****************************************************************
      FORTRAN interface fuer ERRF
*****************************************************************/

#include "math.h" 

double errf_(x)
double *x;
{   
     return erf(*x);
} 
double errfc_(x)
double *x;
{   
     return erfc(*x);
} 
