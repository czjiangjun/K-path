/* 
 *
 * Rot_one[] is rotation for:
 *
 *  CUBIC-I centred
 *  CUBIC-F centred

 *  HEXAGONAL
 *  TRIGONAL

 *  TETRAGONAL-P
 *  TETRAGONAL-I centred

 *  ORTHOROMBIC-A centred
 *  ORTHOROMBIC-F centred
 *
 *  MONOCLINIC-A centred:
 */

double Rot_one[1*3*3]={
                1.0,   0.0,  0.0,
                0.0,   1.0,  0.0,
                0.0,   0.0,  1.0

};

/* Rotations for CUBIC-P  */
double Rot_cub_p[2*3*3]={
                1.0,   0.0,  0.0,
                0.0,   1.0,  0.0,
                0.0,   0.0,  1.0,

                0.0,   1.0,  0.0, /*  Th6 group only changed  */
               -1.0,   0.0,  0.0,
                0.0,   0.0,  1.0
};

/* Rotations for ORTHOROMBIC-P centred */
double Rot_ort_p[6*3*3]={
                1.0,   0.0,  0.0, /*    1; 1  */
                0.0,   1.0,  0.0,
                0.0,   0.0,  1.0,

                0.0,   1.0,  0.0,
                0.0,   0.0,  1.0,
                1.0,   0.0,  0.0,

                0.0,   0.0,  1.0,
                1.0,   0.0,  0.0,
                0.0,   1.0,  0.0,

                0.0,   1.0,  0.0,
               -1.0,   0.0,  0.0,
                0.0,   0.0,  1.0,

                1.0,   0.0,  0.0,
                0.0,   0.0,  1.0,
                0.0,  -1.0,  0.0,

                0.0,   0.0,  1.0,
                0.0,   1.0,  0.0,
               -1.0,   0.0,  0.0

};
/* Rotations for ORTHOROMBIC-I centred */
double Rot_ort_i[4*3*3]={
                1.0,   0.0,  0.0,
                0.0,   1.0,  0.0,
                0.0,   0.0,  1.0,

                0.0,   1.0,  0.0,
                0.0,   0.0,  1.0,
                1.0,   0.0,  0.0,

                0.0,   0.0,  1.0,
                1.0,   0.0,  0.0,
                0.0,   1.0,  0.0,
 
                0.0,   1.0,  0.0,
               -1.0,   0.0,  0.0,
                0.0,   0.0,  1.0
};
/* Rotations for ORTHOROMBIC-C centred */
double Rot_ort_c[2*3*3]={
                1.0,   0.0,  0.0,
                0.0,   1.0,  0.0,
                0.0,   0.0,  1.0,

                0.0,   1.0,  0.0,
               -1.0,   0.0,  0.0,
                0.0,   0.0,  1.0
};

   
/* Rotations for MONOCLINIC-P centred: */
double Rot_mon_p[3*3*3]={
                1.0,  0.0,  0.0, /*    1; 1  */
                0.0,  1.0,  0.0,
                0.0,  0.0,  1.0,
   
                0.0,  1.0,  0.0, /*   15; 4- [ 0 0 1]  */
               -1.0,  0.0,  0.0,
                0.0,  0.0,  1.0,
     
               -1.0,  1.0,  0.0,
               -1.0,  0.0,  0.0,
                0.0,  0.0,  1.0
     
};
