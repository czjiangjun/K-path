#ifndef    WIEN_H
#define    WIEN_H

#include  "type_sg.h"

/* definitions for WIEN package; 
 * hide all Wien parameters into single WIEN variable to avoid a mixture 
 * of WIEN and SGROUP parameters */
typedef struct{

   char   title[256];     /* first line of struct file */
   int    nato;           /* from second line of struct file */
   char   mode[256];        /*  RELA, NREL */
   
/* line 7 of struct file: atom name + NPT + ...  */
   char   atom_name[MAX_ATOMS][10+5+5+5+10+5+10+5+5+1]; /* Rmt, NPT */
   int    isplit[MAX_ATOMS];
   double Z[MAX_ATOMS];
   /* for output only, enumerating subsorts with equal Z */
   int    indZ[MAX_ATOMS];
   
   /* for identification of struct files with wrong unit cell setting */
   double a[2][3];
   int lat; /* for centering mode */   
   t_cell cl; /* input cell from struct file */

   
} t_WIEN;

#endif
