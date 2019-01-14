#include "symbol.inc"
        PROGRAM K_Points_Path

      USE lattice
      USE poscar
      USE base
      USE msymmetry

!      USE fileio
        
      IMPLICIT NONE
      
      TYPE (latt)        LATT_CUR
!      TYPE (latt)        LATT_INI
      TYPE (type_info)   T_INFO
      TYPE (dynamics)    DYN
      TYPE (symmetry)    SYMM

      INTEGER :: IPTYP,NIOND,NIONPD,NTYPD, NTYPPD
      INTEGER :: IPTYP_A
      INTEGER :: IPTYP_B
      INTEGER :: DIR_LEN=0
      CHARACTER(LEN=10) :: DIR_APP
      REAL(q) SAXIS(3), CELDIM(6), CELDIMK(6), PCELDIM(6), PCELDIMK(6)
      INTEGER NCDIJ,NATOMS               !dimension of arrays like CDIJ, CQIJ
      INTEGER :: I,IU0,IU6,IU8

      IU0   = 6
      IU6   = 8
      IU8   = 12

      OPEN(UNIT=IU6,FILE=DIR_APP(1:DIR_LEN)//'OUTKPATH',STATUS='UNKNOWN')
!-----------------------------------------------------------------------      
! read header of POSCAR file to get NTYPD, NTYPDD, NIOND and NIONPD
!-----------------------------------------------------------------------      
      CALL RD_POSCAR_HEAD(LATT_CUR, T_INFO,  &
     &           NIOND,NIONPD, NTYPD,NTYPPD, IU0, IU6)

!      WRITE (IU6,*) 'LATT_CUR=', LATT_CUR
!      WRITE (IU6,*) 'T_INFO=', T_INFO
!      WRITE (IU6,*) 'NIOND=', NIOND,'NIONPD=',NIONPD,'NTYPD=',NTYPD,'NTYPPD=',NTYPPD

!-----------------------------------------------------------------------
! read UNIT=15: POSCAR Startjob and Continuation-job
!-----------------------------------------------------------------------
      CALL RD_POSCAR(LATT_CUR, T_INFO, DYN, &
     &           NIOND,NIONPD, NTYPD,NTYPPD, &
                 IU0, IU6)

!      WRITE (IU6,*) 'LATT_CUR=', LATT_CUR
!      WRITE (IU6,*) 'T_INFO=', T_INFO%SZNAM2
!      WRITE (IU6,*) 'NIOND=', NIOND,'NIONPD=',NIONPD,'NTYPD=',NTYPD,'NTYPPD=',NTYPPD

!-----------------------------------------------------------------------
!  initialize the symmetry stuff
!-----------------------------------------------------------------------
      ALLOCATE(SYMM%ROTMAP(NIOND,1,1), &
               SYMM%TAU(NIOND,3), &
     &         SYMM%TAUROT(NIOND,3),SYMM%WRKROT(3*(NIOND+2)), &
     &         SYMM%PTRANS(NIOND+2,3),SYMM%INDROT(NIOND+2))
      ALLOCATE(DYN%POSION_STD(3,NIOND))

          CALL INISYM(LATT_CUR%B,DYN%POSION,DYN%POSION_STD,DYN%VEL,T_INFO%LSFOR, &
             T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,NIOND,IPTYP, CELDIMK, PCELDIMK, &
             SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
             SYMM%INDROT,T_INFO%ATOMOM,SAXIS,SYMM%MAGROT,NCDIJ,IU6)
           IPTYP_B=IPTYP
!
!           WRITE (IU6,*) 'IPTYP_B=', IPTYP_B
!           WRITE (IU6,*) 'NITYP=', T_INFO%NITYP
!           WRITE (IU6,*) CELDIMK(1:6)
!           WRITE (IU6,*) 'TAU='
!           WRITE (IU6,'(3F12.5)')(SYMM%TAU(I,1:3),I=1,NIOND)

          CALL INISYM(LATT_CUR%A,DYN%POSION,DYN%POSION_STD,DYN%VEL,T_INFO%LSFOR, &
             T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,NIOND,IPTYP, CELDIM, PCELDIM, &
             SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
             SYMM%INDROT,T_INFO%ATOMOM,SAXIS,SYMM%MAGROT,NCDIJ,IU6)
           IPTYP_A=IPTYP

!=======================================================================
! Find SPACE GROUP FROM PRIMITIVE CELL WIEN2k
!=======================================================================
!           WRITE(IU6,*) 'NTYP=', T_INFO%NTYP
!           WRITE(IU6,*) 'NITYP=', T_INFO%NITYP
!           WRITE(IU6,*) 'NIONS=', T_INFO%NIONS
          IF (SYMM%NPTRANS .EQ. 1) THEN
                WRITE(IU6,'(A,/)') "************  SPACE-GROUP ONLY FOR PRIMITIVE CELL  ************"
!                CALL SGROUP(IPTYP_A, CELDIM, T_INFO%NTYP, T_INFO%NITYP, T_INFO%NIONS, DYN%POSION)
          ELSE
                WRITE(IU6,'(A,/)') "!!! PLEASE CHECK: SPACE-GROUP FOR SUPER CELL (MORE THAN 1 PRIMITVIVE CELL) !!!"
!                CALL SGROUP(IPTYP_A, CELDIM, T_INFO%NTYP, T_INFO%NITYP, T_INFO%NIONS, DYN%POSION)
          ENDIF
                CALL SGROUP(IPTYP_A, CELDIM, T_INFO%NTYP, T_INFO%NITYP, T_INFO%NIONS, DYN%POSION_STD)

!           WRITE (IU6,*) 'IPTYP_A=', IPTYP_A
!           WRITE (IU6,*) 'NITYP=', T_INFO%NITYP
!           WRITE (IU6,*) CELDIM(1:6)
!           WRITE (IU6,*) 'TAU='
!           WRITE (IU6,'(3F12.5)')(SYMM%TAU(I,1:3),I=1,NIOND)
           CALL KPATH(IPTYP_B, CELDIM, CELDIMK, IU8)

      DEALLOCATE(SYMM%ROTMAP, SYMM%TAU, SYMM%TAUROT,SYMM%WRKROT, &
     &       SYMM%PTRANS,SYMM%INDROT)
      DEALLOCATE(T_INFO%NITYP, T_INFO%TYPE, T_INFO%LSFOR, &
     &       T_INFO%ITYP, DYN%POSION, DYN%POSION_STD, DYN%POSIOC, DYN%D2C, DYN%VEL, DYN%D2, &
     &       DYN%D3)

      CLOSE (IU6)

        END
