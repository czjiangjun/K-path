               SUBROUTINE KPATH(IPTYP, CELDIM, CELDIMK, IU8)
                       
      IMPLICIT NONE

      INTEGER :: DIR_LEN=0
      CHARACTER(LEN=10) :: DIR_APP
      CHARACTER(LEN=40) :: NAM_BRAVAIS
      INTEGER :: NUMK,IPTYP
      INTEGER :: IU8,GAMMA(3)
      REAL :: CELDIM(6), CELDIMK(8)
      REAL :: alpha,beta,gamma2
      REAL :: min_gamma,max_gamma
      REAL :: v1,v2,v3,v_media
      REAL :: inv_v1,inv_v2,inv_v3
      REAL :: eta,zeta,phi,delta,mu,nu,omega,rho,psi

      REAL  :: K_X(3), K_M(3), K_R(3)
      REAL  :: K_P(3), K_H(3), K_N(3)
      REAL  :: K_U(3), K_K(3), K_W(3), K_L(3)
      REAL  :: K_A(3)
      REAL  :: K_Z(3)
      REAL  :: K_Z1(3)
      REAL  :: K_Y(3), K_Y1(3), SIGMA(3), SIGMA1(3)
      REAL  :: K_B(3), K_P1(3), K_B1(3), K_P2(3), K_F(3), K_Q(3), K_L1(3)
      REAL  :: K_Q1(3)
      REAL  :: K_S(3), K_T(3)
      REAL  :: K_X1(3), K_L2(3)
      REAL  :: K_A1(3), K_C(3), K_H1(3), K_C1(3), K_D(3), K_D1(3)
      REAL  :: K_H2(3), K_M1(3), K_M2(3), K_E(3)
      REAL  :: K_F1(3), K_F2(3), K_F3(3), K_X2(3), K_I(3), K_I1(3), K_Y2(3), &
     &          K_Y3(3), K_N1(3)
      
!      IPTYP=4
      NUMK=10
      GAMMA(1)=0
      GAMMA(2)=0
      GAMMA(3)=0

      OPEN(UNIT=IU8,FILE=DIR_APP(1:DIR_LEN)//'KPOINTS_BAND',STATUS='UNKNOWN')

      IF (IPTYP==1) THEN           ! LATTYP: Found a simple cubic cell.
              NAM_BRAVAIS='Found a simple cubic cell.'
!              NUMK=60

              K_X(1)=0.0; K_X(2)=0.5; K_X(3)=0.0
              K_M(1)=0.5; K_M(2)=0.5; K_M(3)=0.0
              K_R(1)=0.5; K_R(2)=0.5; K_R(3)=0.5

              CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,*)
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,402) K_M(1:3),'!   M'
              WRITE(IU8,*)
              WRITE(IU8,402) K_M(1:3),'!   M'
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,*)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,*)
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,*)
              WRITE(IU8,402) K_M(1:3),'!   M'
              WRITE(IU8,402) K_R(1:3),'!   R'

      ELSE IF (IPTYP==2) THEN       ! LATTYP: Found a body centered cubic cell.
              NAM_BRAVAIS='Found a body centered cubic cell.'
!              NUMK=60

              K_P(1)=0.25; K_P(2)=0.25; K_P(3)=0.0
              K_H(1)=0.5; K_H(2)=-0.5; K_H(3)=0.5
              K_N(1)=0.0; K_N(2)=0.0; K_N(3)=0.5

              CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_H(1:3),'!   H'
              WRITE(IU8,*)
              WRITE(IU8,402) K_H(1:3),'!   H'
              WRITE(IU8,402) K_N(1:3),'!   N'
              WRITE(IU8,*)
              WRITE(IU8,402) K_N(1:3),'!   N'
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,*)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_P(1:3),'!   P'
              WRITE(IU8,*)
              WRITE(IU8,402) K_P(1:3),'!   P'
              WRITE(IU8,402) K_H(1:3),'!   H'
              WRITE(IU8,*)
              WRITE(IU8,402) K_P(1:3),'!   P'
              WRITE(IU8,402) K_N(1:3),'!   N'

      ELSE IF (IPTYP==3) THEN       ! LATTYP: Found a face centered cubic cell.
              NAM_BRAVAIS='Found a face centered cubic cell.'
!              NUMK=100

              K_U(1)=0.625; K_U(2)=0.25; K_U(3)=0.625
              K_K(1)=0.375; K_K(2)=0.375; K_K(3)=0.75
              K_W(1)=0.5; K_W(2)=0.25; K_W(3)=0.75
              K_L(1)=0.5; K_L(2)=0.5; K_L(3)=0.5
              K_X(1)=0.5; K_X(2)=0.0; K_X(3)=0.5

              CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,*)
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,402) K_W(1:3),'!   W'
              WRITE(IU8,*)
              WRITE(IU8,402) K_W(1:3),'!   W'
              WRITE(IU8,402) K_K(1:3),'!   K'
              WRITE(IU8,*)
              WRITE(IU8,402) K_K(1:3),'!   K'
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,*)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_L(1:3),'!   L'
              WRITE(IU8,*)
              WRITE(IU8,402) K_L(1:3),'!   L'
              WRITE(IU8,402) K_U(1:3),'!   U'
              WRITE(IU8,*)
              WRITE(IU8,402) K_U(1:3),'!   U'
              WRITE(IU8,402) K_U(1:3),'!   W'
              WRITE(IU8,*)
              WRITE(IU8,402) K_U(1:3),'!   W'
              WRITE(IU8,402) K_L(1:3),'!   L'
              WRITE(IU8,*)
              WRITE(IU8,402) K_L(1:3),'!   L'
              WRITE(IU8,402) K_K(1:3),'!   K'
              WRITE(IU8,*)
              WRITE(IU8,402) K_U(1:3),'!   U'
              WRITE(IU8,402) K_X(1:3),'!   X'

      ELSE IF (IPTYP==4) THEN       ! LATTYP: Found a hexagonal cell.
              NAM_BRAVAIS='Found a hexagonal cell.'
!              NUMK=90

              K_K(1)=1.0/3.0; K_K(2)=1.0/3.0; K_K(3)=0.0
              K_A(1)=0.0; K_A(2)=0.0; K_A(3)=0.5
              K_L(1)=0.5; K_L(2)=0.0; K_L(3)=0.5
              K_H(1)=1.0/3.0; K_H(2)=1.0/3.0; K_H(3)=0.5
              K_M(1)=0.5; K_M(2)=0.0; K_M(3)=0.0

              CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_M(1:3),'!   M'
              WRITE(IU8,*)
              WRITE(IU8,402) K_M(1:3),'!   M'
              WRITE(IU8,402) K_K(1:3),'!   K'
              WRITE(IU8,*)
              WRITE(IU8,402) K_K(1:3),'!   K'
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,*)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_A(1:3),'!   A'
              WRITE(IU8,*)
              WRITE(IU8,402) K_A(1:3),'!   A'
              WRITE(IU8,402) K_L(1:3),'!   L'
              WRITE(IU8,*)
              WRITE(IU8,402) K_L(1:3),'!   L'
              WRITE(IU8,402) K_H(1:3),'!   H'
              WRITE(IU8,*)
              WRITE(IU8,402) K_H(1:3),'!   H'
              WRITE(IU8,402) K_A(1:3),'!   A'
              WRITE(IU8,*)
              WRITE(IU8,402) K_L(1:3),'!   L'
              WRITE(IU8,402) K_M(1:3),'!   M'
              WRITE(IU8,*)
              WRITE(IU8,402) K_K(1:3),'!   K'
              WRITE(IU8,402) K_H(1:3),'!   H'

      ELSE IF (IPTYP==5) THEN       ! LATTYP: Found a simple tetragonal cell.
              NAM_BRAVAIS='Found a simple tetragonal cell.'
!              NUMK=90

              K_R(1)=0.0; K_R(2)=0.5; K_R(3)=0.5
              K_A(1)=0.5; K_A(2)=0.5; K_A(3)=0.5
              K_L(1)=0.0; K_L(2)=0.5; K_L(3)=0.0
              K_M(1)=0.5; K_M(2)=0.5; K_M(3)=0.0
              K_Z(1)=0.0; K_Z(2)=0.0; K_Z(3)=0.5

              CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,*)
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,402) K_M(1:3),'!   M'
              WRITE(IU8,*)
              WRITE(IU8,402) K_M(1:3),'!   M'
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,*)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,*)
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,402) K_A(1:3),'!   A'
              WRITE(IU8,*)
              WRITE(IU8,402) K_A(1:3),'!   A'
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,*)
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,*)
              WRITE(IU8,402) K_M(1:3),'!   M'
              WRITE(IU8,402) K_A(1:3),'!   A'

      ELSE IF (IPTYP==6) THEN       ! LATTYP: Found a body centered tetragonal cell.
              NAM_BRAVAIS='Found a body centered tetragonal cell.'
              IF (CELDIM(3) .LT. 1.0) THEN     ! C/A-ratio < 1.0 
!                 NUMK=90
                 eta=(1.0+CELDIM(3)**2)/4.0

                 K_X(1)=0.0; K_X(2)=0.0; K_X(3)=0.5
                 K_M(1)=-0.5; K_M(2)=0.5; K_M(3)=0.5
                 K_Z(1)=eta; K_Z(2)=eta; K_Z(3)=-eta
                 K_N(1)=0.0; K_N(2)=0.5; K_N(3)=0.0
                 K_Z1(1)=-eta; K_Z1(2)=1-eta; K_Z1(3)=eta
                 K_P(1)=0.25; K_P(2)=0.25;  K_P(3)=0.25

                 CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,402) K_X(1:3),'!   X'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_X(1:3),'!   X'
                 WRITE(IU8,402) K_M(1:3),'!   M'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_M(1:3),'!   M'
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,*)
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
                 WRITE(IU8,402) K_P(1:3),'!   P'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_P(1:3),'!   P'
                 WRITE(IU8,402) K_N(1:3),'!   N'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_N(1:3),'!   N'
                 WRITE(IU8,402) K_Z1(1:3),'!   Z_1'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_Z1(1:3),'!   Z_1'
                 WRITE(IU8,402) K_M(1:3),'!   M'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_X(1:3),'!   X'
                 WRITE(IU8,402) K_P(1:3),'!   P'
              ELSE                             ! C/A-ratio > 1.0
!                 NUMK=110
                 eta=(1.0+(1.0/CELDIM(3))**2)/4.0
                 zeta=(1.0/CELDIM(3))**2/2.0

                 K_X(1)=0.0; K_X(2)=0.0; K_X(3)=0.5
                 K_N(1)=0.0; K_N(2)=0.5; K_N(3)=0.0
                 K_Y(1)=-zeta; K_Y(2)=zeta; K_Y(3)=0.5
                 K_P(1)=0.25; K_P(2)=0.25;  K_P(3)=0.25
                 K_Y1(1)=-0.5; K_Y1(2)=0.5; K_Y1(3)=-zeta
                 SIGMA(1)=-eta; SIGMA(2)=eta; SIGMA(3)=eta
                 K_Z(1)=0.5; K_Z(2)=0.5; K_X(3)=-0.5
                 SIGMA1(1)=eta; SIGMA1(2)=1-eta; SIGMA1(3)=-eta

                 CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,402) K_X(1:3),'!   X'
                 WRITE(IU8,*)
                 WRITE(IU8,401) K_X(1:3),'!   X'
                 WRITE(IU8,402) K_Y(1:3),'!   Y'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_Y(1:3),'!   Y'
                 WRITE(IU8,402) SIGMA(1:3),'!   sigma'
                 WRITE(IU8,*)
                 WRITE(IU8,402) SIGMA(1:3),'!   sigma'
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,*)
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
                 WRITE(IU8,402) SIGMA1(1:3),'!   sigma_1'
                 WRITE(IU8,*)
                 WRITE(IU8,402) SIGMA1(1:3),'!   sigma_1'
                 WRITE(IU8,402) K_N(1:3),'!   N'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_N(1:3),'!   N'
                 WRITE(IU8,402) K_P(1:3),'!   P'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_P(1:3),'!   P'
                 WRITE(IU8,402) K_Y1(1:3),'!   Y_1'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_Y1(1:3),'!   Y_1'
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_X(1:3),'!   X'
                 WRITE(IU8,402) K_P(1:3),'!   P'
              ENDIF

      ELSE IF (IPTYP==7) THEN       ! LATTYP: Found a trigonal (rhoboedric) cell.
              NAM_BRAVAIS='Found a trigonal (rhoboedric) cell.'
!              NUMK=90
              IF (CELDIM(4) .GT. 0.0) THEN      ! COS(alpha) > 0.0
                 eta=(1.0+4.0*CELDIM(4))/(2.0+4.0*CELDIM(4))
                 nu=3.0/4.0-eta/2.0

                 K_P(1)=eta; K_P(2)=nu;  K_P(3)=nu
                 K_B(1)=eta; K_B(2)=0.5;  K_B(3)=1.0-eta
                 K_P1(1)=1-nu; K_P1(2)=1-nu; K_P1(3)=1.0-eta
                 K_B1(1)=0.5; K_B1(2)=1-eta;  K_B1(3)=eta-1.0
                 K_P2(1)=nu; K_P2(2)=nu; K_P2(3)=eta-1.0
                 K_F(1)=0.5; K_F(2)=0.5; K_F(3)=0.0
                 K_Q(1)=1-nu; K_Q(2)=nu; K_Q(3)=0.0
                 K_L(1)=0.5; K_L(2)=0.0; K_L(3)=0.0
                 K_X(1)=nu; K_X(2)=0.0; K_Q(3)=-nu
                 K_L1(1)=0.0; K_L1(2)=0.0; K_L1(3)=-0.5
                 K_Z(1)=0.5; K_Z(2)=0.5; K_Z(3)=0.5

                 CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,402) K_L(1:3),'!   L'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_L(1:3),'!   L'
                 WRITE(IU8,402) K_B1(1:3),'!   B_1'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_B(1:3),'!   B'
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,*)
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,402) K_X(1:3),'!   X'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_Q(1:3),'!   Q'
                 WRITE(IU8,402) K_F(1:3),'!   F'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_F(1:3),'!   F'
                 WRITE(IU8,402) K_P1(1:3),'!   P_1'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_P1(1:3),'!   P_1'
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_L(1:3),'!   L'
                 WRITE(IU8,402) K_P(1:3),'!   P'
              ELSE                             ! COS(alpha) < 0.0
                 alpha=ACOS(CELDIM(4))
                 eta=1.0/(2.0*(TAN(alpha/2.0))**2)
                 nu=3.0/4.0-eta/2.0

                 K_P1(1)=nu; K_P1(2)=nu-1.0; K_P1(3)=nu-1.0
                 K_F(1)=0.5; K_F(2)=-0.5; K_F(3)=0.0
                 K_Q(1)=eta; K_Q(2)=eta; K_Q(3)=eta
                 K_L(1)=0.5; K_L(2)=0.0; K_L(3)=0.0
                 K_Q1(1)=1.0-eta; K_Q1(2)=-eta; K_Q1(3)=-eta
                 K_P(1)=1.0-nu; K_P(2)=-nu;  K_P(3)=1.0-nu
                 K_Z(1)=0.5; K_Z(2)=-0.5; K_Z(3)=0.5

                 CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,402) K_P(1:3),'!   P'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_P(1:3),'!   P'
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
                 WRITE(IU8,402) K_Q(1:3),'!   Q'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_Q(1:3),'!   Q'
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,*)
                 WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                 WRITE(IU8,402) K_F(1:3),'!   F'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_F(1:3),'!   F'
                 WRITE(IU8,402) K_P1(1:3),'!   P_1'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_P1(1:3),'!   P_1'
                 WRITE(IU8,402) K_Q1(1:3),'!   Q_1'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_Q1(1:3),'!   Q_1'
                 WRITE(IU8,402) K_L(1:3),'!   L'
                 WRITE(IU8,*)
                 WRITE(IU8,402) K_L(1:3),'!   L'
                 WRITE(IU8,402) K_Z(1:3),'!   Z'
              ENDIF

      ELSE IF (IPTYP==8) THEN       ! LATTYP: Found a simple orthorhombic cell.
              !!!! should be checked !!!!
              NAM_BRAVAIS='Found a orthorhombic cell.'

              v1=CELDIM(1); v2=CELDIM(2)*v1; v3=CELDIM(3)*v1

!              NUMK=120

              K_U(1)=0.5; K_U(2)=0.0; K_U(3)=0.5
              K_R(1)=0.5; K_R(2)=0.5; K_R(3)=0.5
              K_X(1)=0.5; K_X(2)=0.0; K_X(3)=0.0
              K_S(1)=0.5; K_S(2)=0.5; K_S(3)=0.0
              K_Y(1)=0.0; K_Y(2)=0.5; K_Y(3)=0.0
              K_T(1)=0.0; K_T(2)=0.5; K_T(3)=0.5
              K_Z(1)=0.0; K_Z(2)=0.0; K_Z(3)=0.5

              CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,*)
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,402) K_S(1:3),'!   S'
              WRITE(IU8,*)
              WRITE(IU8,402) K_S(1:3),'!   S'
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,*)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,402) K_U(1:3),'!   U'
              WRITE(IU8,*)
              WRITE(IU8,402) K_U(1:3),'!   U'
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,*)
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,402) K_T(1:3),'!   T'
              WRITE(IU8,*)
              WRITE(IU8,402) K_T(1:3),'!   T'
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,402) K_T(1:3),'!   T'
              WRITE(IU8,*)
              WRITE(IU8,402) K_U(1:3),'!   U'
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,*)
              WRITE(IU8,402) K_S(1:3),'!   S'
              WRITE(IU8,402) K_R(1:3),'!   R'

      ELSE IF (IPTYP==9) THEN       ! LATTYP: Found a body centered orthorhombic cell.
              !!!! should be checked !!!!
              NAM_BRAVAIS='Found a body centered orthorhombic cell.'

              v1=CELDIM(1); v2=CELDIM(2)*v1; v3=CELDIM(3)*v1

!              NUMK=130
!              IF (CELDIM(2) .GE. 1.0)
              zeta=(1.0+CELDIM(3)**2)/4.0
              delta=(CELDIM(2)**2-CELDIM(3)**2)/4.0
              eta=(1.0+CELDIM(2)**2)/4.0
              mu=(CELDIM(2)**2+CELDIM(3)**2)/4.0

              K_W(1)=0.25; K_W(2)=0.25; K_W(3)=0.25
              K_L(1)=-mu; K_L(2)=mu; K_L(3)=0.5-delta
              K_X(1)=-zeta; K_X(2)=zeta; K_X(3)=zeta
              K_L1(1)=mu; K_L1(2)=-mu; K_L1(3)=0.5+delta
              K_X1(1)=zeta; K_X1(2)=1-zeta; K_X1(3)=-zeta
              K_L2(1)=0.5-delta; K_L2(2)=0.5+delta; K_L2(3)=-mu
              K_Y(1)=eta; K_Y(2)=-eta; K_Y(3)=eta
              K_R(1)=0.0; K_R(2)=0.5; K_R(3)=0.0
              K_Y1(1)=1-eta; K_Y1(2)=eta; K_Y1(3)=-eta
              K_S(1)=0.5; K_S(2)=0.0; K_S(3)=0.0
              K_Z(1)=0.5; K_Z(2)=0.5; K_Z(3)=-0.5
              K_T(1)=0.0; K_T(2)=0.0; K_T(3)=0.5

              CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,*)
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,402) K_L(1:3),'!   L'
              WRITE(IU8,*)
              WRITE(IU8,402) K_L(1:3),'!   L'
              WRITE(IU8,402) K_T(1:3),'!   T'
              WRITE(IU8,*)
              WRITE(IU8,402) K_T(1:3),'!   T'
              WRITE(IU8,402) K_W(1:3),'!   W'
              WRITE(IU8,*)
              WRITE(IU8,402) K_W(1:3),'!   W'
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,*)
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,402) K_X(1:3),'!   X_1'
              WRITE(IU8,*)
              WRITE(IU8,402) K_X(1:3),'!   X_1'
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,*)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,402) K_S(1:3),'!   S'
              WRITE(IU8,*)
              WRITE(IU8,402) K_S(1:3),'!   S'
              WRITE(IU8,402) K_W(1:3),'!   W'
              WRITE(IU8,*)
              WRITE(IU8,402) K_L1(1:3),'!   L_1'
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Y1(1:3),'!   Y_1'
              WRITE(IU8,402) K_Z(1:3),'!   Z'

      ELSE IF (IPTYP==10) THEN       ! LATTYP: Found a face centered orthorhombic cell.
              !!!! should be checked !!!!
              NAM_BRAVAIS='Found a face centered orthorhombic cell.'

              v1=CELDIM(1); v2=CELDIM(2)*v1; v3=CELDIM(3)*v1
              inv_v1=1.0/(v1**2); inv_v2=1.0/(v2**2); inv_v3=1.0/(v3**2);

              IF (inv_v1 .GT. (inv_v2+inv_v3)) THEN
!                NUMK=110

                zeta=(1.0+(v1/v2)**2-(v1/v3)**2)/4.0
                eta=(1.0+(v1/v2)**2+(v1/v3)**2)/4.0

                K_X(1)=0.0; K_X(2)=eta; K_X(3)=eta
                K_A(1)=0.5; K_X(2)=0.5+zeta; K_X(3)=zeta
                K_X1(1)=1.0; K_X1(2)=1.0-eta; K_X1(3)=1.0-eta
                K_A1(1)=0.5; K_A1(2)=0.5-zeta; K_A1(3)=1.0-zeta
                K_Y(1)=0.5; K_Y(2)=0.0; K_Y(3)=0.5
                K_L(1)=0.5; K_L(2)=0.5; K_L(3)=0.5
                K_Z(1)=0.5; K_Z(2)=0.5; K_Z(3)=0.0
                K_T(1)=1.0; K_T(2)=0.5; K_T(3)=0.0

                CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,402) K_T(1:3),'!   T'
                WRITE(IU8,*)
                WRITE(IU8,402) K_T(1:3),'!   T'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_X(1:3),'!   Z'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Z(1:3),'!   X'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Z(1:3),'!   X'
                WRITE(IU8,402) K_A1(1:3),'!   A_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_A1(1:3),'!   A_1'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_T(1:3),'!   T'
                WRITE(IU8,402) K_X1(1:3),'!   X_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,402) K_A(1:3),'!   A'
                WRITE(IU8,*)
                WRITE(IU8,402) K_A(1:3),'!   A'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_L(1:3),'!   L'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              ELSEIF (inv_v1 .LT. (inv_v2+inv_v3)) THEN
!                NUMK=130

                eta=(1.0+(v1/v2)**2-(v1/v3)**2)/4.0
                delta=(1.0+(v1/v2)**2-(v2/v3)**2)/4.0
                phi=(1.0+(v3/v2)**2-(v3/v1)**2)/4.0

                K_H(1)=1.0-phi; K_H(2)=0.5-phi; K_H(3)=0.5
                K_C(1)=0.5; K_C(2)=0.5-eta; K_C(3)=1.0-eta
                K_H1(1)=phi; K_H1(2)=0.5+phi; K_H1(3)=0.5
                K_C1(1)=0.5; K_C1(2)=0.5+eta; K_C1(3)=eta
                K_X(1)=0.0; K_X(2)=0.5; K_X(3)=0.5
                K_D(1)=0.5-delta; K_D(2)=0.5; K_D(3)=1.0-delta
                K_Y(1)=0.5; K_Y(2)=0.0; K_Y(3)=0.5
                K_D1(1)=0.5+delta; K_D1(2)=0.5; K_D1(3)=delta
                K_Z(1)=0.5; K_Z(2)=0.5; K_Z(3)=0.0
                K_L(1)=0.5; K_L(2)=0.5; K_L(3)=0.5

                CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,402) K_C(1:3),'!   C'
                WRITE(IU8,*)
                WRITE(IU8,402) K_C(1:3),'!   C'
                WRITE(IU8,402) K_D(1:3),'!   D'
                WRITE(IU8,*)
                WRITE(IU8,402) K_D(1:3),'!   D'
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,*)
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,402) K_D1(1:3),'!   D_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_D1(1:3),'!   D_1'
                WRITE(IU8,402) K_H(1:3),'!   H'
                WRITE(IU8,*)
                WRITE(IU8,402) K_H(1:3),'!   H'
                WRITE(IU8,402) K_C(1:3),'!   C'
                WRITE(IU8,*)
                WRITE(IU8,402) K_C1(1:3),'!   C_1'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,402) K_H1(1:3),'!   H_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_H(1:3),'!   H'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_L(1:3),'!   L'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              ELSE
!                NUMK=100

                zeta=(1.0+(v1/v2)**2-(v1/v3)**2)/4.0
                eta=(1.0+(v1/v2)**2+(v1/v3)**2)/4.0

                K_X(1)=0.0; K_X(2)=eta; K_X(3)=eta
                K_A(1)=0.5; K_X(2)=0.5+zeta; K_X(3)=zeta
                K_X1(1)=1.0; K_X1(2)=1.0-eta; K_X1(3)=1.0-eta
                K_A1(1)=0.5; K_A1(2)=0.5-zeta; K_A1(3)=1.0-zeta
                K_Y(1)=0.5; K_Y(2)=0.0; K_Y(3)=0.5
                K_L(1)=0.5; K_L(2)=0.5; K_L(3)=0.5
                K_Z(1)=0.5; K_Z(2)=0.5; K_Z(3)=0.0
                K_T(1)=1.0; K_T(2)=0.5; K_T(3)=0.0

                CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,402) K_T(1:3),'!   T'
                WRITE(IU8,*)
                WRITE(IU8,402) K_T(1:3),'!   T'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_X(1:3),'!   Z'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Z(1:3),'!   X'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Z(1:3),'!   X'
                WRITE(IU8,402) K_A1(1:3),'!   A_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_A1(1:3),'!   A_1'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,402) K_A(1:3),'!   A'
                WRITE(IU8,*)
                WRITE(IU8,402) K_A(1:3),'!   A'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_L(1:3),'!   L'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              ENDIF


      ELSE IF (IPTYP==11) THEN       ! LATTYP: Found a base centered orthorhombic cell.
              !!!! should be checked !!!!
              NAM_BRAVAIS='Found a base centered orthorhombic cell.'

              v1=CELDIM(1); v2=CELDIM(2)*v1; v3=CELDIM(3)*v1

!              NUMK=120

              zeta=(1.0+(v1/v2)**2)/4.0

              K_T(1)=-0.5; K_T(2)=0.5; K_T(3)=0.5
              K_A(1)=zeta; K_A(2)=zeta; K_A(3)=0.5
              K_X(1)=zeta; K_X(2)=zeta; K_X(3)=0.0
              K_A1(1)=-zeta; K_A1(2)=1.0-zeta; K_A1(3)=0.5
              K_X1(1)=-zeta; K_X1(2)=1.0-zeta; K_X1(3)=0.0
              K_R(1)=0.0; K_R(2)=0.5; K_R(3)=0.5
              K_Y(1)=-0.5; K_Y(2)=0.5; K_Y(3)=0.0
              K_S(1)=0.0; K_S(2)=0.5; K_S(3)=0.0
              K_Z(1)=0.0; K_Z(2)=0.0; K_Z(3)=0.5

              CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,*)
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,402) K_S(1:3),'!   S'
              WRITE(IU8,*)
              WRITE(IU8,402) K_S(1:3),'!   S'
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,*)
              WRITE(IU8,402) K_R(1:3),'!   R'
              WRITE(IU8,402) K_A(1:3),'!   A'
              WRITE(IU8,*)
              WRITE(IU8,402) K_A(1:3),'!   A'
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,*)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,402) K_X1(1:3),'!   X_1'
              WRITE(IU8,*)
              WRITE(IU8,402) K_X1(1:3),'!   X_1'
              WRITE(IU8,402) K_A1(1:3),'!   A_1'
              WRITE(IU8,*)
              WRITE(IU8,402) K_A1(1:3),'!   A_1'
              WRITE(IU8,402) K_T(1:3),'!   T'
              WRITE(IU8,*)
              WRITE(IU8,402) K_T(1:3),'!   T'
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,402) K_T(1:3),'!   T'

      ELSE IF (IPTYP==12) THEN       ! LATTYP: Found a simple monoclinic cell.
              !!!! should be checked !!!!
              NAM_BRAVAIS='Found a simple monoclinic cell.'

              v1=CELDIM(1); v2=CELDIM(2)*v1; v3=CELDIM(3)*v1
              beta=ACOS(CELDIM(5))

!              NUMK=110

              eta=(1.0-v2*CELDIM(5)/v3)*(2.0*SIN(beta)**2)
              nu=0.5-eta*v3*CELDIM(5)/v2

              K_H2(1)=0.0; K_H2(2)=eta; K_H2(3)=-nu
              K_A(1)=0.5; K_A(2)=0.5; K_A(3)=0.0
              K_M(1)=0.5; K_M(2)=eta; K_M(3)=1.0-nu
              K_C(1)=0.0; K_C(2)=0.5; K_C(3)=0.5
              K_M1(1)=0.5; K_M1(2)=1.0-eta; K_M1(3)=nu
              K_D(1)=0.5; K_D(2)=0.0; K_D(3)=0.5
              K_M2(1)=0.5; K_M2(2)=eta; K_M2(3)=-nu
              K_D1(1)=0.5; K_D1(2)=0.0; K_D1(3)=-0.5
              K_X(1)=0.0; K_X(2)=0.5; K_X(3)=0.0
              K_E(1)=0.5; K_E(2)=0.5; K_E(3)=0.5
              K_Y(1)=0.0; K_Y(2)=0.0; K_Y(3)=0.5
              K_H(1)=0.0; K_H(2)=eta; K_H(3)=1.0-nu
              K_Y1(1)=0.0; K_Y1(2)=0.0; K_Y1(3)=-0.5
              K_H1(1)=0.0; K_H1(2)=1.0-eta; K_H1(3)=nu
              K_Z(1)=0.5; K_Z(2)=0.0; K_Z(3)=0.0

              CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
              WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,402) K_H(1:3),'!   H'
              WRITE(IU8,*)
              WRITE(IU8,402) K_H(1:3),'!   H'
              WRITE(IU8,402) K_C(1:3),'!   C'
              WRITE(IU8,*)
              WRITE(IU8,402) K_C(1:3),'!   C'
              WRITE(IU8,402) K_E(1:3),'!   E'
              WRITE(IU8,*)
              WRITE(IU8,402) K_E(1:3),'!   E'
              WRITE(IU8,402) K_M1(1:3),'!   M_1'
              WRITE(IU8,*)
              WRITE(IU8,402) K_M1(1:3),'!   M_1'
              WRITE(IU8,401) K_A(1:3),'!   A'
              WRITE(IU8,*)
              WRITE(IU8,401) K_A(1:3),'!   A'
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,*)
              WRITE(IU8,402) K_X(1:3),'!   X'
              WRITE(IU8,402) K_H1(1:3),'!   H_1'
              WRITE(IU8,*)
              WRITE(IU8,402) K_M(1:3),'!   M'
              WRITE(IU8,402) K_D(1:3),'!   D'
              WRITE(IU8,*)
              WRITE(IU8,402) K_D(1:3),'!   D'
              WRITE(IU8,402) K_Z(1:3),'!   Z'
              WRITE(IU8,*)
              WRITE(IU8,402) K_Y(1:3),'!   Y'
              WRITE(IU8,402) K_D(1:3),'!   D'

      ELSE IF (IPTYP==13) THEN       ! LATTYP: Found a base centered monoclinic cell.
              NAM_BRAVAIS='Found a base centered monoclinic cell.'
              !!!! should be checked !!!!

              v1=CELDIM(1); v2=CELDIM(2)*v1; v3=CELDIM(3)*v1
              beta=ACOS(CELDIM(5))
              v_media=v2/v3*CELDIM(5)+(v2*SIN(beta)/v1)**2-1.0

              IF (CELDIMK(6) .LT. 0.0) THEN
!                NUMK=100

                zeta=(2.0-v2/v3*CELDIM(5))/(4.0*SIN(beta)**2)
                eta=0.5+2.0*zeta*v3/v2*CELDIM(5)
                psi=0.75-v1**2/(4.0*v2**2*SIN(beta)**2)
                phi=psi+(0.75-psi)*v2*CELDIM(5)/v3

                K_L(1)=0.5; K_L(2)=0.5; K_L(3)=0.5
                K_N(1)=0.5; K_N(2)=0.0; K_N(3)=0.0
                K_M(1)=0.5; K_M(2)=0.0; K_M(3)=0.5
                K_N1(1)=0.0; K_N1(2)=-0.5; K_N(3)=0.0
                K_X(1)=1.0-psi; K_X(2)=psi-1.0; K_X(3)=0.0
                K_F(1)=1.0-zeta; K_F(2)=1.0-zeta; K_F(3)=1.0-eta
                K_X1(1)=psi; K_X1(2)=1.0-psi; K_X1(3)=0.0
                K_F1(1)=zeta; K_F1(2)=zeta; K_F1(3)=eta
                K_X2(1)=psi-1.0; K_X2(2)=-psi; K_X2(3)=0.0
                K_F2(1)=-zeta; K_F2(2)=-zeta; K_F2(3)=1.0-eta
                K_Y(1)=0.5; K_Y(2)=0.5; K_Y(3)=0.0
                K_F3(1)=1.0-zeta; K_F3(2)=-zeta; K_F3(3)=1.0-eta
                K_Y1(1)=-0.5; K_Y1(2)=-0.5; K_Y1(3)=0.0
                K_I(1)=phi; K_I(2)=1.0-phi; K_I(3)=0.5
                K_Z(1)=0.0; K_Z(2)=0.0; K_Z(3)=0.5
                K_I1(1)=1.0-phi; K_I1(2)=phi-1.0; K_I1(3)=0.5

                CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,402) K_F(1:3),'!   F'
                WRITE(IU8,*)
                WRITE(IU8,402) K_F(1:3),'!   F'
                WRITE(IU8,402) K_L(1:3),'!   L'
                WRITE(IU8,*)
                WRITE(IU8,402) K_L(1:3),'!   L'
                WRITE(IU8,402) K_I(1:3),'!   I'
                WRITE(IU8,*)
                WRITE(IU8,402) K_I1(1:3),'!   I_1'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,402) K_F1(1:3),'!   F_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,402) K_X1(1:3),'!   X_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_N(1:3),'!   N'
                WRITE(IU8,*)
                WRITE(IU8,402) K_M(1:3),'!   M'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
              ELSEIF (ABS(CELDIMK(6)) .LT. 1.0E-7) THEN
!                NUMK=80

                zeta=(2.0-v2/v3*CELDIM(5))/(4.0*SIN(beta)**2)
                eta=0.5+2.0*zeta*v3/v2*CELDIM(5)
                psi=0.75-v1**2/(4.0*v2**2*SIN(beta)**2)
                phi=psi+(0.75-psi)*v2*CELDIM(5)/v3

                K_L(1)=0.5; K_L(2)=0.5; K_L(3)=0.5
                K_N(1)=0.5; K_N(2)=0.0; K_N(3)=0.0
                K_M(1)=0.5; K_M(2)=0.0; K_M(3)=0.5
                K_N1(1)=0.0; K_N1(2)=-0.5; K_N(3)=0.0
                K_X(1)=1.0-psi; K_X(2)=psi-1.0; K_X(3)=0.0
                K_F(1)=1.0-zeta; K_F(2)=1.0-zeta; K_F(3)=1.0-eta
                K_X1(1)=psi; K_X1(2)=1.0-psi; K_X1(3)=0.0
                K_F1(1)=zeta; K_F1(2)=zeta; K_F1(3)=eta
                K_X2(1)=psi-1.0; K_X2(2)=-psi; K_X2(3)=0.0
                K_F2(1)=-zeta; K_F2(2)=-zeta; K_F2(3)=1.0-eta
                K_Y(1)=0.5; K_Y(2)=0.5; K_Y(3)=0.0
                K_F3(1)=1.0-zeta; K_F3(2)=-zeta; K_F3(3)=1.0-eta
                K_Y1(1)=-0.5; K_Y1(2)=-0.5; K_Y1(3)=0.0
                K_I(1)=phi; K_I(2)=1.0-phi; K_I(3)=0.5
                K_Z(1)=0.0; K_Z(2)=0.0; K_Z(3)=0.5
                K_I1(1)=1.0-phi; K_I1(2)=phi-1.0; K_I1(3)=0.5

                CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,402) K_F(1:3),'!   F'
                WRITE(IU8,*)
                WRITE(IU8,402) K_F(1:3),'!   F'
                WRITE(IU8,402) K_L(1:3),'!   L'
                WRITE(IU8,*)
                WRITE(IU8,402) K_L(1:3),'!   L'
                WRITE(IU8,402) K_I(1:3),'!   I'
                WRITE(IU8,*)
                WRITE(IU8,402) K_I1(1:3),'!   I_1'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,402) K_F1(1:3),'!   F_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_N(1:3),'!   N'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_M(1:3),'!   M'

              ELSEIF ((ABS(CELDIMK(6)) .GT. 1.0E-7) .AND. (v_media .LT. 0.0)) THEN
!                NUMK=110

                mu=(1.0+(v2/v1)**2)/4.0
                delta=v2*v3*CELDIM(5)/(2.0*v1**2)
                zeta=mu-0.25+(1.0-v2/v3*CELDIM(5))/(4.0*SIN(beta)**2)
                eta=0.5+2.0*zeta*v3/v2*CELDIM(5)
                phi=1.0+zeta-2.0*mu
                psi=eta-2.0*delta

                K_N(1)=0.5; K_N(2)=0.0; K_N(3)=0.0
                K_F(1)=1.0-phi; K_F(2)=1.0-phi; K_F(3)=1.0-psi
                K_N1(1)=0.0; K_N1(2)=-0.5; K_N1(3)=0.0
                K_F1(1)=phi; K_F1(2)=phi-1.0; K_F1(3)=psi
                K_X(1)=0.5; K_X(2)=-0.5; K_X(3)=0.0
                K_F2(1)=1.0-phi; K_F2(2)=-phi; K_F2(3)=1.0-psi
                K_Y(1)=mu; K_Y(2)=mu; K_Y(3)=delta
                K_H(1)=zeta; K_H(2)=zeta; K_H(3)=eta
                K_Y1(1)=1.0-mu; K_Y1(2)=-mu; K_Y1(3)=-delta
                K_H1(1)=1.0-zeta; K_H1(2)=-zeta; K_H1(3)=1.0-eta
                K_Y2(1)=-mu; K_Y2(2)=-mu; K_Y2(3)=-delta
                K_H2(1)=-zeta; K_H2(2)=-zeta; K_H2(3)=1.0-eta
                K_Y3(1)=mu; K_Y3(2)=mu-1.0; K_Y3(3)=delta
                K_I(1)=0.5; K_I(2)=-0.5; K_I(3)=0.5
                K_Z(1)=0.0; K_Z(2)=0.5; K_Z(3)=0.5
                K_M(1)=0.5; K_M(2)=0.0; K_M(3)=0.5

                CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,402) K_F(1:3),'!   F'
                WRITE(IU8,*)
                WRITE(IU8,402) K_F(1:3),'!   F'
                WRITE(IU8,402) K_H(1:3),'!   H'
                WRITE(IU8,*)
                WRITE(IU8,402) K_H(1:3),'!   H'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,402) K_I(1:3),'!   I'
                WRITE(IU8,*)
                WRITE(IU8,402) K_I(1:3),'!   I'
                WRITE(IU8,402) K_F1(1:3),'!   F_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_H1(1:3),'!   H_1'
                WRITE(IU8,402) K_Y1(1:3),'!   Y_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y1(1:3),'!   Y_1'
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,*)
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_N(1:3),'!   N'
                WRITE(IU8,*)
                WRITE(IU8,402) K_M(1:3),'!   M'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'

              ELSEIF ((ABS(CELDIMK(6)) .GT. 1.0E-7) .AND. (ABS(v_media) .LT. 1.0E-7)) THEN
!                NUMK=100

                mu=(1.0+(v2/v1)**2)/4.0
                delta=v2*v3*CELDIM(5)/(2.0*v1**2)
                zeta=mu-0.25+(1.0-v2/v3*CELDIM(5))/(4.0*SIN(beta)**2)
                eta=0.5+2.0*zeta*v3/v2*CELDIM(5)
                phi=1.0+zeta-2.0*mu
                psi=eta-2.0*delta

                K_N(1)=0.5; K_N(2)=0.0; K_N(3)=0.0
                K_F(1)=1.0-phi; K_F(2)=1.0-phi; K_F(3)=1.0-psi
                K_N1(1)=0.0; K_N1(2)=-0.5; K_N1(3)=0.0
                K_F1(1)=phi; K_F1(2)=phi-1.0; K_F1(3)=psi
                K_X(1)=0.5; K_X(2)=-0.5; K_X(3)=0.0
                K_F2(1)=1.0-phi; K_F2(2)=-phi; K_F2(3)=1.0-psi
                K_Y(1)=mu; K_Y(2)=mu; K_Y(3)=delta
                K_H(1)=zeta; K_H(2)=zeta; K_H(3)=eta
                K_Y1(1)=1.0-mu; K_Y1(2)=-mu; K_Y1(3)=-delta
                K_H1(1)=1.0-zeta; K_H1(2)=-zeta; K_H1(3)=1.0-eta
                K_Y2(1)=-mu; K_Y2(2)=-mu; K_Y2(3)=-delta
                K_H2(1)=-zeta; K_H2(2)=-zeta; K_H2(3)=1.0-eta
                K_Y3(1)=mu; K_Y3(2)=mu-1.0; K_Y3(3)=delta
                K_I(1)=0.5; K_I(2)=-0.5; K_I(3)=0.5
                K_Z(1)=0.0; K_Z(2)=0.5; K_Z(3)=0.5
                K_M(1)=0.5; K_M(2)=0.0; K_M(3)=0.5

                CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,402) K_F(1:3),'!   F'
                WRITE(IU8,*)
                WRITE(IU8,402) K_F(1:3),'!   F'
                WRITE(IU8,402) K_H(1:3),'!   H'
                WRITE(IU8,*)
                WRITE(IU8,402) K_H(1:3),'!   H'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,402) K_I(1:3),'!   I'
                WRITE(IU8,*)
                WRITE(IU8,402) K_H1(1:3),'!   H_1'
                WRITE(IU8,402) K_Y1(1:3),'!   Y_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y1(1:3),'!   Y_1'
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,*)
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_N(1:3),'!   N'
                WRITE(IU8,*)
                WRITE(IU8,402) K_M(1:3),'!   M'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'

              ELSEIF ((ABS(CELDIMK(6)) .GT. 1.0E-7) .AND. (v_media .GT. 0.0)) THEN
!                NUMK=120
                
                zeta=((v2/v1)**2+(1.0-v2/v1*CELDIM(5))/SIN(beta)**2)/4.0
                eta=0.5+2.0*zeta*v3/v2*CELDIM(5)
                mu=eta/2.0+(v2/2.0*v1)**2-v2*v3*CELDIM(5)/(2.0*v1**2)
                nu=2.0*mu-zeta
                rho=1.0-zeta*(v1/v2)**2
                omega=(4.0*nu-1.0-(v2*SIN(beta)/v1)**2)*v3/(2.0*v2*CELDIM(5))
                delta=zeta*v3*CELDIM(5)/v2+omega/2.0-0.25

                K_M(1)=0.5; K_M(2)=0.0; K_M(3)=0.5
                K_F(1)=nu; K_F(2)=nu; K_F(3)=omega
                K_N(1)=0.5; K_N(2)=0.0; K_N(3)=0.0
                K_F1(1)=1.0-nu; K_F1(2)=1.0-nu; K_F1(3)=1.0-omega
                K_N1(1)=0.0; K_N1(2)=-0.5; K_N1(3)=0.0
                K_F2(1)=nu; K_F2(2)=nu-1.0; K_F2(3)=omega
                K_X(1)=0.5; K_X(2)=-0.5; K_X(3)=0.0
                K_H(1)=zeta; K_H(2)=zeta; K_H(3)=eta
                K_Y(1)=mu; K_Y(2)=mu; K_Y(3)=delta
                K_H1(1)=1.0-zeta; K_H1(2)=-zeta; K_H1(3)=1.0-eta
                K_Y1(1)=1.0-mu; K_Y1(2)=-mu; K_Y1(3)=-delta
                K_H2(1)=-zeta; K_H2(2)=-zeta; K_H2(3)=1.0-eta
                K_Y2(1)=-mu; K_Y2(2)=-mu; K_Y2(3)=-delta
                K_I(1)=rho; K_I(2)=1.0-rho; K_I(3)=0.5
                K_Y3(1)=mu; K_Y3(2)=mu-1.0; K_Y3(3)=delta
                K_I1(1)=1.0-rho; K_I1(2)=rho-1.0; K_I1(3)=0.5
                K_Z(1)=0.0; K_Z(2)=0.0; K_Z(3)=0.5
                K_L(1)=0.5; K_L(2)=0.5; K_L(3)=0.5

                CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,402) K_F(1:3),'!   F'
                WRITE(IU8,*)
                WRITE(IU8,402) K_F(1:3),'!   F'
                WRITE(IU8,402) K_L(1:3),'!   L'
                WRITE(IU8,*)
                WRITE(IU8,402) K_L(1:3),'!   L'
                WRITE(IU8,402) K_I(1:3),'!   I'
                WRITE(IU8,*)
                WRITE(IU8,402) K_I1(1:3),'!   I_1'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,402) K_H(1:3),'!   H'
                WRITE(IU8,*)
                WRITE(IU8,402) K_H(1:3),'!   H'
                WRITE(IU8,402) K_F1(1:3),'!   F_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_H1(1:3),'!   H_1'
                WRITE(IU8,402) K_Y1(1:3),'!   Y_1'
                WRITE(IU8,*)
                WRITE(IU8,402) K_Y1(1:3),'!   Y_1'
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,*)
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_N(1:3),'!   N'
                WRITE(IU8,*)
                WRITE(IU8,402) K_M(1:3),'!   M'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'

              ENDIF

      ELSE IF (IPTYP==14) THEN       ! LATTYP: Found a triclinic cell.
              NAM_BRAVAIS='Found a triclinic cell.'
              !!!! should be checked !!!!
!              NUMK=70

              max_gamma=MAX(CELDIMK(4),CELDIMK(5),CELDIMK(6))
              min_gamma=MIN(CELDIMK(4),CELDIMK(5),CELDIMK(6))
              IF (max_gamma .LE. 1.0E-7) THEN
                K_R(1)=0.5; K_R(2)=0.5; K_R(3)=0.5
                K_L(1)=0.5; K_L(2)=0.5; K_L(3)=0.0
                K_X(1)=0.5; K_X(2)=0.0; K_X(3)=0.0
                K_M(1)=0.0; K_M(2)=0.5; K_M(3)=0.5
                K_Y(1)=0.0; K_Y(2)=0.5; K_Y(3)=0.0
                K_N(1)=0.5; K_N(2)=0.0; K_N(3)=0.5
                K_Z(1)=0.0; K_Z(2)=0.0; K_Z(3)=0.5

!                IF (max_gamma .LT. -1.0E-7) THEN
!
!                ELSE
!
!                ENDIF
              ELSEIF (min_gamma .GE. -1.0E-7) THEN
                K_R(1)=0.0; K_R(2)=-0.5; K_R(3)=0.5
                K_L(1)=0.5; K_L(2)=-0.5; K_L(3)=0.0
                K_X(1)=0.0; K_X(2)=-0.5; K_X(3)=0.0
                K_M(1)=0.0; K_M(2)=0.0; K_M(3)=0.5
                K_Y(1)=0.5; K_Y(2)=0.0; K_Y(3)=0.0
                K_N(1)=-0.5; K_N(2)=-0.5; K_N(3)=0.5
                K_Z(1)=-0.5; K_Z(2)=0.0; K_Z(3)=0.5

!                IF (min_gamma .GT. 1.0E-7) THEN
!
!                ELSE
!
!                ENDIF
              ENDIF

                CALL WRITE_HEAD(NUMK,NAM_BRAVAIS,IU8)
                WRITE(IU8,402) K_X(1:3),'!   X'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Y(1:3),'!   Y'
                WRITE(IU8,*)
                WRITE(IU8,402) K_L(1:3),'!   L'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_Z(1:3),'!   Z'
                WRITE(IU8,*)
                WRITE(IU8,402) K_N(1:3),'!   N'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,*)
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'
                WRITE(IU8,402) K_M(1:3),'!   M'
                WRITE(IU8,*)
                WRITE(IU8,402) K_R(1:3),'!   R'
                WRITE(IU8,401) GAMMA(1:3),'!   \Gamma'

      ENDIF

  401 FORMAT(3(I2,6X),A)
  402 FORMAT(3(F6.3,2X),A)

      CLOSE (IU8)

      RETURN
      END SUBROUTINE

           SUBROUTINE WRITE_HEAD(NUMK,NAME_BRA,IU8)
      IMPLICIT NONE

      INTEGER :: NUMK,IU8
      CHARACTER(LEN=40) :: NAME_BRA
      
      WRITE(IU8,'(A37,A32)') 'k-points along high symmetry lines:  ', NAME_BRA
      WRITE(IU8,400) NUMK, '! ', NUMK, 'intersections'
      WRITE(IU8,'(A9)')'Line-mode'
      WRITE(IU8,'(A10)')'reciprocal'

  400 FORMAT(2(I4,4X,A))
      RETURN
      END SUBROUTINE

