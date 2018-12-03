#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "sto.h"
#include "type_sg.h"

int stoi(char *s, int *iar, int max_n)
{
  int j=0,k,n=0;
  char str[256];

  while (1)
  {
     if( n == max_n ) return n;
     while (1) {
       if( s[j] == '\0' || s[j] == '/' || s[j] == '\n' ) return n; 
       if( isdigit((int) s[j]) || s[j]=='+' || s[j]=='-' ) break;
       j++;
    }
    k=0;
    while ( isdigit((int) s[j]) || s[j]=='+' || s[j]=='-')  str[k++]=s[j++];
    str[k]='\0';
    iar[n++]=atoi(str);
  }
}

int stod(char *s, double *dar, int max_n)
{
  int j=0,k,n=0;
  char str[256];

  while (1)
  {
    if( n == max_n ) return n;
    while (1) {
       if( s[j] == '\0' || s[j] == '/' || s[j] == '\n' ) return n; 
       if( isdigit((int) s[j]) || s[j]=='+' || s[j]=='-' ||
	   s[j]=='.' || tolower((int) s[j])=='e') break;
       j++;
    }
    k=0;
    while ( isdigit((int) s[j]) || s[j]=='.' || tolower((int) s[j])=='e' ||
            s[j]=='+' || s[j]=='-')  str[k++]=s[j++];
    str[k]='\0';
    dar[n++]=atof(str);
  }
}

int stos(char *s, char *s_out)
{
  int j=0,k;

  while ( isspace((int) s[j]) ) j++;
  if( s[j] == '\0' || s[j] == '/' ) return 0;
  k=0;
  while ( !( isspace((int) s[j]) || s[j] == '\0' || s[j] == '/') )  
     s_out[k++]=s[j++];
  s_out[k]='\0';
  return 1;
}


/* print a floating value with FORTRAN Fw.d specification,
 * not exactly the same as Fortran behaviour */
char *sprintFwd( char *s, double f, int w, int d )
{
   int d_new;
   char fmt[64];
   
   d_new=d;
   while(1) {
     if( d_new < 0 ) return NULL;
     sprintf(fmt,"%c%d%c%d%c", '%',w,'.',d_new,'f');
     sprintf(s,fmt,f);
     /* check the length of the string */
     if( strlen(s) > w ) d_new--;
     else return s;
   }
}


void print_dbl(char *text, int n, double *d)
{
   int i;
   if( text !=NULL ) printf("%s",text);
   for(i=0; i<n; i++)
     printf(" % .8f",d[i]);
   printf("\n");
   return;
}

/* '\n' -> 0 */
int n_to_zero( char *s )
{
   int i;
   
   for(i=0; s[i] != '\0'; i++)
   if( s[i] == '\n' ) {
       s[i] = '\0';
       return i;
   }
   
   return -1;
}
