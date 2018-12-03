#ifndef STO_INC
#define STO_INC

extern int  stoi(char *s, int *iar, int max_n);
extern int  stod(char *s, double *dar, int max_n);
extern int  stos(char *s, char *s_out);
extern char *sprintFwd( char *s, double f, int w, int d );
extern void print_dbl(char *text, int n, double *d);
extern int  n_to_zero( char *s );
#endif
