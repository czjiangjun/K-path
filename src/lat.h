#ifndef LAT_H
#define LAT_H

extern int is_transl( t_cell *cl, double tr[3],
                      double T[3][3], double T_1[3][3]);

extern int det_lat_NSM(int *lat, double T1[3], double T2[3], double T3[3],
                      double T1p[3], double T2p[3], double T3p[3],
	              double T1t[3], double T2t[3], double T3t[3] );

extern int basis_transl(t_cell *cl_in, t_cell *cl, int *lat, double Tin[3][3],
                         double Tp[3][3], double T[3][3], double Ttrcl[3][3]
			);

#endif
