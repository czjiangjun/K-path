#ifndef    TYPE_SG_H
#define    TYPE_SG_H

#define  MAX_ATOMS    4096
#define  MAX_CHARS    32

/* default tolerance for floating arithmetic */
#define  TOLER        1.e-4

#ifndef NULL
#define NULL ((void *)0)
#endif

#ifndef PI
#define PI 3.14159265358979323844
#endif

extern  double TOL;
/* type definition for cell */
typedef struct {
   int     nat;                    /*  number of atoms  */
   int     nsrt;                   /*  number of species  */
   double  r[MAX_ATOMS][3];        /*  coordinates  */
   int     srt[MAX_ATOMS];         /*  sort for each atom  */
   int     natsrt[MAX_ATOMS];      /*  number of atoms for each sort */

}  t_cell;

/* command line switches */
typedef struct {
   int  is_noeq; /* print only atoms of primitive cell in con. basis*/
   int  is_prim; /* print atoms of pritive cell in prim. basis */
   
   int  is_wien_input;  /* read data from WIEN struct file*/
   int  is_wien_output; /* write data to WIEN struct file*/
}  t_switch;


#define DEBUG 0

#endif
