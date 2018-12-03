#ifndef IO_H
#define IO_H

#include "wien.h"

/* read data from WIEN input file */
extern int read_WIEN_data(char *fname, t_cell *cl_in,
              char name_srt[MAX_ATOMS][MAX_CHARS],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
              int *lat, t_WIEN *WIEN);

extern int read_data(char *fname, t_cell *cl_in, char name_srt[MAX_ATOMS][MAX_CHARS],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
	      int *lat);

extern int write_res(char *fname, t_cell *cl, char name_srt[MAX_ATOMS][MAX_CHARS],
  	      t_cell *cl_new, char name_srt_new[MAX_ATOMS][MAX_CHARS],
	      t_switch *sw, int fix_org,
	      int nsh, double Rsh[][3],
              double a[2][3], /*  a,b,c, alpha,beta,gamma  */
              int lat, char *lat_name[], char * sgrp_name,
              int npgrp, int nsgrp, int nop, double sym_op[][4][3],
              char name_pgrp[][32],
	      double Tin[3][3], t_WIEN *WIEN);
  
/* write out results to WIEN struct file */
int write_WIEN_struct(
    char *fname, t_cell *cl, char name_srt[MAX_ATOMS][MAX_CHARS],
    t_cell *cl_new, char name_srt_new[MAX_ATOMS][MAX_CHARS],
    int fix_org, int nsh, double Rsh[][3],
    double a[2][3], /*  a,b,c, alpha,beta,gamma  */
    int lat, char *lat_name[], char * sgrp_name,
    int npgrp, int nop, double sym_op[][4][3],
    char name_pgrp[][32],
    double Tin[3][3], t_WIEN *WIEN);

#endif
