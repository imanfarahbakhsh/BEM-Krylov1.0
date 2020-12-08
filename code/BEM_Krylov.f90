!     Last change:  IMAN  9 Dec 2020   11:24 pm
!
! This program solves the two dimensional (LA)place equation
! using the (B)oundary (E)lement method with (CON) stant
! boundary elements
!******************************************************************************
!******************************************************************************
!******************************************************************************
MODULE BEM
IMPLICIT NONE
INTEGER,PARAMETER          :: N    =500            ! N= Number of boundary elements which is equal to the number
                                                   ! of boundary nodes
INTEGER,PARAMETER          :: INN  =((N/4)-1)**2   ! IN= Number of internal points where the function u is calculated
INTEGER,PARAMETER          :: AUXN =N/4
INTEGER,PARAMETER          :: IM   =(N/4)-1
INTEGER,PARAMETER          :: JM   =(N/4)-1
INTEGER,PARAMETER          :: ITMAX=100
REAL(8),PARAMETER          :: PI   =3.1415926535897932384626433832795
REAL(8),PARAMETER          :: LX   =1.D0!2*PI
REAL(8),PARAMETER          :: LY   =1.D0!2*PI
REAL(8),PARAMETER          :: DX   =LX/AUXN
REAL(8),PARAMETER          :: DY   =LY/AUXN
REAL(8),PARAMETER          :: TOL  =1.D-10
END MODULE
!******************************************************************************
!******************************************************************************
!******************************************************************************
PROGRAM LABECON
USE BEM
IMPLICIT NONE
INTEGER                 :: LSING
INTEGER,DIMENSION(N)    :: INDEX
REAL(8),DIMENSION(N)    :: XM,YM,UB,UNB
REAL(8),DIMENSION(N+1)  :: XL,YL
REAL(8),DIMENSION(N,N)  :: G,H,A
REAL(8),DIMENSION(INN)  :: XIN,YIN,UIN
!------------------------------------------------------------------------------
CALL INPUT_PRODUCTION (XL,YL,XIN,YIN,INDEX,UB)  ! Read the data from INPUTFILE
CALL GMATR (XL,YL,XM,YM,G)                      ! Compute the G matrix
CALL HMATR (XL,YL,XM,YM,H)                      ! Compute the H matrix
CALL ABMATR (G,H,A,UNB,UB,INDEX)                ! Form the system of equations AX=B
CALL SOLVEQ (A,UNB,LSING)                       ! Solve the system of equations
CALL REORDER (UB,UNB,INDEX)                     ! Form the vectors U and UN of all the boundary values
CALL UINTER (XL,YL,XIN,YIN,UB,UNB,UIN)          ! Compute the values UIN of u at the internal points
CALL OUTPUT (XIN,YIN,UIN)          ! Print the results in the OUTPUTFILE
END

!******************************************************************************
!******************************************************************************
!******************************************************************************
!
! This subroutine computes the elements of the G matrix
!
SUBROUTINE GMATR (XL,YL,XM, YM,G)
USE BEM
IMPLICIT NONE
INTEGER                 ::I,J,JP1
REAL(8),DIMENSION(N+1)  ::XL,YL
REAL(8),DIMENSION(N)    ::XM,YM
REAL(8),DIMENSION(N,N)  ::G
REAL(8)                 ::RESULT
!
! Compute the nodal coordinates and store them in the arrays
! XM and YM
!
XL(N+1)=XL(1)
YL(N+1)=YL(1)
DO 10 I=1,N
XM(I)=(XL(I)+XL(I+1))/2.D0
10 YM(I)=(YL(I)+YL(I+1))/2.D0
!
! Compute the elements of matrix G
!
DO 20 I=1,N
DO 20 J=1,N
JP1=J+1
IF (I.NE.J) THEN
CALL RLINTC (XM(I),YM(I),XL(J),YL(J),XL(JP1),YL(JP1),RESULT)
G(I,J)=RESULT
ELSEIF (I.EQ.J) THEN
CALL SLINTC (XL(J),YL(J),XL(JP1),YL(JP1),RESULT)
G(I,J)=RESULT
ENDIF
20 CONTINUE
RETURN
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
!
! This subroutine computes the off-diagonal elements of the
! matrix G
!
! RA= The distance of point 0 from the Gauss integration point
! on the boundary element
!
! WG= The weights of the Gauss integration
!
! XI= The coordinates of the Gauss integration points in the
! interval [-1,1]
!
! XC,YC= The global coordinates of the Gauss integration points
!
SUBROUTINE RLINTC (X0,Y0,X1,Y1,X2,Y2,RESULT)
USE BEM
IMPLICIT NONE
INTEGER               ::I
REAL(8),DIMENSION(4)  ::XC,YC,XI,WG
REAL(8)               ::RESULT,X0,Y0,X1,X2,Y1,Y2,AX,AY,BX,BY,RA,SL
DATA XI/-0.86113631,-0.33998104,0.33998104,0.86113631/
DATA WG/0.34785485,0.65214515,0.65214515,0.34785485/
!------------------------------------------------------------------------------
AX=(X2-X1)/2.D0
AY=(Y2-Y1)/2.D0
BX=(X2+X1)/2.D0
BY=(Y2+Y1)/2.D0
!
!Compute the line integral
!
RESULT=0.D0
!
DO 30 I=1,4
XC(I)=AX*XI(I)+BX
YC(I)=AY*XI(I)+BY
RA=SQRT((XC(I)-X0)**2+(YC(I)-Y0)**2)
30 RESULT=RESULT+DLOG(RA)*WG(I)
SL=2.D0*SQRT(AX**2+AY**2)
RESULT=RESULT*SL/(4.D0*PI)
RETURN
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
!This subroutines computes the diagonal elements of the matrix G
!
SUBROUTINE SLINTC (X1,Y1,X2,Y2,RESULT)
USE BEM
IMPLICIT NONE
REAL(8)        ::AX,AY,X2,Y2,X1,Y1,SL,RESULT
!------------------------------------------------------------------------------
AX=(X2-X1)/2.D0
AY=(Y2-Y1)/2.D0
SL=SQRT(AX**2+AY**2)
RESULT=SL*(DLOG(SL)-1.D0)/PI
RETURN
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
!
!This subroutine computes the elements of the H matrix
!
SUBROUTINE HMATR (XL,YL,XM,YM,H)
USE BEM
IMPLICIT NONE
INTEGER                 ::I,J
REAL(8),DIMENSION(N+1)  ::XL,YL
REAL(8),DIMENSION(N)    ::XM,YM
REAL(8),DIMENSION(N,N)  ::H
REAL(8)                 ::RESULT
!------------------------------------------------------------------------------
XL(N+1)=XL(1)
YL(N+1)=YL(1)
DO 10 I=1,N
XM(I)=(XL(I)+XL(I+1))/2.D0
10 YM(I)=(YL(I)+YL(I+1))/2.D0
!
!Compute the elements of the H matrix
!
DO 20 I=1,N
DO 20 J=1,N
IF (I.NE.J) THEN
CALL DALPHA (XM(I),YM(I),XL(J),YL(J),XL(J+1),YL(J+1),RESULT)
H(I,J)=RESULT
ELSEIF (I.EQ.J) THEN
H(I,J)=-0.5D0
ENDIF
20 CONTINUE
RETURN
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
!
!This subroutine computes the off diagonal elements of
!the H matrix
!
SUBROUTINE DALPHA (X0,Y0,X1,Y1,X2,Y2,RESULT)
USE BEM
IMPLICIT NONE
REAL(8)       ::COS1,SIN1,X0,Y0,X1,Y1,DX2R,DY2R
REAL(8)       ::X2,Y2,RESULT,DY1,DX1,DX2,DY2,DL1,DA
!------------------------------------------------------------------------------
DY1=Y1-Y0
DX1=X1-X0
DY2=Y2-Y0
DX2=X2-X0
DL1=SQRT(DX1**2+DY1**2)
COS1=DX1/DL1
SIN1=DY1/DL1
DX2R=DX2*COS1+DY2*SIN1
DY2R=-DX2*SIN1+DY2*COS1
DA=ATAN2(DY2R,DX2R)
RESULT=DA/(2.D0*PI)
RETURN
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
!
!This subroutine rearranges the matrices G and H and produces
!the matrices A and B=UNB
!
SUBROUTINE ABMATR (G,H,A,UNB,UB,INDEX)
USE BEM
IMPLICIT NONE
INTEGER                ::I,J
REAL(8),DIMENSION(N,N) ::G,H,A
REAL(8),DIMENSION(N)   ::UB,UNB
INTEGER,DIMENSION(N)   ::INDEX
!
!Reorder the columns of the system of equations and
!store them in A
!
DO 40 J=1,N
IF (INDEX(J).EQ.0) THEN
DO 20 I=1,N
20 A(I,J)=-G(I,J)
ELSEIF (INDEX(J).NE.0) THEN
DO 30 I=1,N
30 A(I,J)=H(I,J)
END IF
40 CONTINUE
!
!Compute the right hand side vector and store it in UNB
!
DO 50 I=1,N
UNB(I)=0.D0
DO 60 J=1,N
IF (INDEX(J).EQ.0) THEN
UNB(I)=UNB(I)-H(I,J)*UB(J)
ELSEIF (INDEX(J).NE.0) THEN
UNB(I)=UNB(I)+G(I,J)*UB(J)
ENDIF
60 CONTINUE
50 CONTINUE
RETURN
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE SOLVEQ (A,UNB,LSING)
USE BEM
IMPLICIT NONE
INTEGER                ::LSING
REAL(8),DIMENSION(N,N) ::A
REAL(8),DIMENSION(N)   ::UNB
REAL(8)                ::BEGIN_TIME,END_TIME
!------------------------------------------------------------------------------
CALL CPU_TIME(BEGIN_TIME)
!------------------------------------------------------------------------------
!CALL LEQS     (A,UNB,LSING)
!CALL BCGSTAB  (A,UNB)
CALL CGS      (A,UNB)
!CALL BICG     (A,UNB)
!CALL GMRES_STAR (A,UNB)
!------------------------------------------------------------------------------
CALL CPU_TIME(END_TIME)
PRINT*,'Total time which spent in solver is',END_TIME-BEGIN_TIME,'sec.'
!------------------------------------------------------------------------------
IF (LSING.EQ.0) THEN
WRITE (2,150)
150 FORMAT(/,'',69('*')//2X'The system has been solved regularly'/)
ELSEIF (LSING.EQ.1) THEN
WRITE (2,170)
170 FORMAT (/,'',69('*')//2X'The system is singular'/)
ENDIF
RETURN
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
!
! This subroutine uses Gauss elimination to solve
! a system of linear equations, [A]{X}={B}, where
! A : One-dimensional array which contains the occasional row of
! the two-dimensional array of the coefficients of the unknowns
! B : One-dimensional array which contains the known coefficients
! N : Integer denoting the number of the unknowns
! LSING: Integer taking the values:
! LSING = 0, if the system has been solved regularly
! LSING = I, if the system is singular
!
SUBROUTINE LEQS (A,B,LSING)
USE BEM
IMPLICIT NONE
INTEGER              ::I,J,LSING,JJ,JY,IHELP,IJ,IMAX,I1,I2,K
INTEGER              ::NY,IQS,IX,IXJ,JX,IJREF,JJX,NN,I3
REAL(8)              ::AMAX,ATEMP
REAL(8),DIMENSION(1) ::A,B
!
LSING=0
JJ=-N
DO 10 J=1,N
JY=J+1
JJ=JJ+N+1
AMAX=0.D0
IHELP=JJ-J
DO 20 I=J,N
IJ=IHELP+I
IF(ABS(AMAX)-ABS(A(IJ))) 30,20,20
30 AMAX=A(IJ)
IMAX=I
20 CONTINUE
IF (ABS(AMAX).EQ.0.D0) THEN
LSING=1
RETURN
END IF
I1=J+N*(J-2)
IHELP=IMAX-J
DO 40 K=J,N
I1=I1+N
I2=I1+IHELP
ATEMP=A(I1)
A(I1)=A(I2)
A(I2)=ATEMP
40 A(I1)=A(I1)/AMAX
ATEMP=B(IMAX)
B(IMAX)=B(J)
B(J)=ATEMP/AMAX
IF (J-N) 50,70,50
50 IQS=N*(J-1)
DO 10 IX=JY,N
IXJ=IQS+IX
IHELP=J-IX
DO 60 JX=JY,N
IJREF=N*(JX-1)+IX
JJX=IJREF+IHELP
60 A(IJREF)=A(IJREF)-(A(IXJ)*A(JJX))
10 B(IX)=B(IX)-B(J)*A(IXJ)
70 NY=N-1
NN=N*N
DO 80 J=1,NY
I1=NN-J
I2=N-J
I3=N
DO 80 K=1,J
B(I2)=B(I2)-A(I1)*B(I3)
I1=I1-N
80 I3=I3-1
RETURN
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
!
!This subroutine rearranges the arrays UB and UNB in such a way
!that all values of the function u are stored in UB while all
!the values of the normal derivative un are stored in UNB
!
SUBROUTINE REORDER (UB,UNB,INDEX)
USE BEM
IMPLICIT NONE
INTEGER               ::I
REAL(8),DIMENSION(N)  ::UB,UNB
INTEGER,DIMENSION(N)  ::INDEX
REAL(8)               ::CH
!
DO 20 I=1,N
IF (INDEX (I)) 20,20,10
10 CH=UB(I)
UB(I)=UNB(I)
UNB(I)=CH
20 CONTINUE
RETURN
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
!
!This subroutine computes the values of u at the internal points
!
SUBROUTINE UINTER (XL,YL,XIN,YIN,UB,UNB,UIN)
USE BEM
IMPLICIT NONE
INTEGER                  ::J,K,JP1
REAL(8),DIMENSION(N+1)   ::XL,YL
REAL(8),DIMENSION(INN)   ::XIN,YIN,UIN
REAL(8),DIMENSION(N)     ::UB,UNB
REAL(8)                  ::RESG,RESH
!
!Compute the values of u at the internal points
!
DO 10 K=1,INN
UIN(K)=0.D0
DO 20 J=1,N
JP1=J+1
CALL DALPHA (XIN(K),YIN(K),XL(J),YL(J),XL(JP1),YL(JP1),RESH)
CALL RLINTC (XIN(K),YIN(K),XL(J),YL(J),XL(JP1),YL(JP1),RESG)
20 UIN(K)=UIN(K)+RESH*UB(J)-RESG*UNB(J)
10 CONTINUE
RETURN
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
!This subroutine prints the results in the output file.
!
SUBROUTINE OUTPUT (XIN,YIN,UIN)
USE BEM
IMPLICIT NONE
INTEGER                :: K
REAL(8),DIMENSION(INN) :: XIN,YIN,UIN
!------------------------------------------------------------------------------
OPEN (UNIT=2,FILE="/results/OUTPUTFILE.PLT")     ! Open the output file
WRITE (2,*)'VARIABLES=X,Y,U'
WRITE (2,*)'ZONE'
WRITE (2,*)'F=POINT'
WRITE (2,*)'I=',IM
WRITE (2,*)'J=',JM
DO K=1,INN
WRITE(2,400)XIN(K),YIN(K),UIN(K)
400 FORMAT(E14.5,2X,E14.5,2X,E14.5)
END DO
CLOSE(2)
END
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE INPUT_PRODUCTION (XBC,YBC,X,Y,INDEX,U)
USE BEM
IMPLICIT NONE
INTEGER                  :: I,J
INTEGER,DIMENSION(N)     :: INDEX
REAL(8),DIMENSION(N)     :: XBC,YBC,U
REAL(8),DIMENSION(IM,JM) :: X,Y
!------------------------------------------------------------------------------
DO I=1,N
SELECT CASE (I)
CASE (1:AUXN+1)
XBC(I)=(I-1)*DX
YBC(I)=0.D0
CASE (AUXN+2:2*AUXN+1)
XBC(I)=LX
YBC(I)=(I-AUXN-1)*DY
CASE (2*AUXN+2:3*AUXN+1)
XBC(I)=(3*AUXN-I+1)*DX
YBC(I)=LY
CASE (3*AUXN+2:N)
XBC(I)=0.D0
YBC(I)=(N-I+1)*DY
END SELECT
END DO
DO I=1,N
SELECT CASE (I)
    CASE (1:AUXN)
    U(I)=100.D0
    INDEX(I)=0
    CASE (AUXN+1:2*AUXN)
    U(I)=0.D0
    INDEX(I)=1
    CASE (2*AUXN+1:3*AUXN)
    U(I)=0.D0
    INDEX(I)=1
    CASE (3*AUXN+1:N)
    U(I)=50.D0
    INDEX(I)=0
END SELECT
END DO
X(1,:)=DX
DO J=1,JM
DO I=1,IM-1
X(I+1,J)=X(I,J)+DX
END DO
END DO
Y(:,1)=DY
DO I=1,IM
DO J=1,JM-1
Y(I,J+1)=Y(I,J)+DY
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE INPUT_PRODUCTION
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                         BCGSTAB_SUBROUTINE
!______________________________________________________________________________
SUBROUTINE BCGSTAB (A,B)
USE BEM
IMPLICIT NONE
INTEGER                            :: I,ITER
REAL(8),DIMENSION(N)               :: B,X_OLD,X,P_OLD,P,R_OLD,R,AP,AQQ,RS,QQ
REAL(8)                            :: NORM,S,ALPHA,BETA,M,MS,MMS,MM,MN,OM
REAL(8),DIMENSION(N,N)             :: A
!------------------------------------------------------------------------------
R_OLD=B
RS=R_OLD
P_OLD=R_OLD
NORM=1.D0
DO WHILE (NORM.GT.TOL)
ITER=ITER+1
!PRINT*,ITER
!------------------------------------------------------------------------------
DO I=1,N
AP(I)=DOT_PRODUCT(A(I,:),P_OLD)
END DO
!------------------------------------------------------------------------------
M=DOT_PRODUCT(R_OLD,RS)
MM=DOT_PRODUCT(AP,RS)
ALPHA=M/MM
QQ=R_OLD-ALPHA*AP
!------------------------------------------------------------------------------
DO I=1,N
AQQ(I)=DOT_PRODUCT(A(I,:),QQ)
END DO
!------------------------------------------------------------------------------
MS=DOT_PRODUCT(AQQ,QQ)
MMS=DOT_PRODUCT(AQQ,AQQ)
OM=MS/MMS
X=X_OLD+ALPHA*P_OLD+OM*QQ
R=QQ-OM*AQQ
S=DOT_PRODUCT(R,R)
NORM=SQRT(S)/N
MN=DOT_PRODUCT(R,RS)
BETA=(MN/M)*(ALPHA/OM)
P=R+BETA*(P_OLD-OM*AP)
P_OLD=P
X_OLD=X
R_OLD=R
!PRINT*,NORM
CALL SUCCESSIVE_SOL ("BCGSTAB",NORM,ITER)
END DO
B=X
!------------------------------------------------------------------------------
END SUBROUTINE BCGSTAB
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                     SUCCESSIVE_SOLUTION SUBROUTINE
!______________________________________________________________________________
SUBROUTINE SUCCESSIVE_SOL (METHOD,ERR,ITER)
IMPLICIT NONE
INTEGER                             :: ITER
CHARACTER(*)                        :: METHOD
CHARACTER*30                        :: FNAME
REAL(8)                             :: ERR
!------------------------------------------------------------------------------
FNAME="/results/CONV-"//METHOD//".PLT"
!------------------------------------------------------------------------------
OPEN(3,FILE=FNAME)
WRITE(3,*) ITER,ERR
!------------------------------------------------------------------------------
END SUBROUTINE SUCCESSIVE_SOL
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                         CGS_SUBROUTINE
!______________________________________________________________________________
SUBROUTINE CGS (A,B)
USE BEM
IMPLICIT NONE
INTEGER                            :: I,ITER
REAL(8),DIMENSION(N,N)             :: A
REAL(8),DIMENSION(N)               :: B,X_OLD,X,P_OLD,P,R_OLD,U,U_OLD,R,AP,AUQ,RS,Q
REAL(8)                            :: NORM,S,ALPHA,BETA,M,MM,MN
!------------------------------------------------------------------------------
R_OLD=B
RS=R_OLD
P_OLD=R_OLD
U_OLD=R_OLD
NORM=1.D0
DO WHILE (NORM.GT.TOL)
ITER=ITER+1
!PRINT*,ITER
!------------------------------------------------------------------------------
DO I=1,N
AP(I)=DOT_PRODUCT(A(I,:),P_OLD)
END DO
!------------------------------------------------------------------------------
M=DOT_PRODUCT(R_OLD,RS)
MM=DOT_PRODUCT(AP,RS)
ALPHA=M/MM
Q=U_OLD-ALPHA*AP
X=X_OLD+ALPHA*(U_OLD+Q)
!------------------------------------------------------------------------------
DO I=1,N
AUQ(I)=DOT_PRODUCT(A(I,:),(U_OLD+Q))
END DO
!------------------------------------------------------------------------------
R=R_OLD-ALPHA*AUQ
S=DOT_PRODUCT(R,R)
NORM=SQRT(S)/N
MN=DOT_PRODUCT(R,RS)
BETA=MN/M
U=R+BETA*Q
P=U+BETA*(Q+BETA*P_OLD)
P_OLD=P
X_OLD=X
R_OLD=R
U_OLD=U
!PRINT*,NORM
CALL SUCCESSIVE_SOL ("CGS",NORM,ITER)
END DO
B=X
!------------------------------------------------------------------------------
END SUBROUTINE CGS
!______________________________________________________________________________
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                         BCG_SUBROUTINE
!_____________________________________________________________________________
SUBROUTINE BICG (A,B)
USE BEM
IMPLICIT NONE
INTEGER                            :: I,J,ITER
REAL(8),DIMENSION(N)               :: X_OLD,X,P_OLD,P,R_OLD,R,B,&
                                      AP,RS_OLD,PS_OLD,PS,RS,ATPS
REAL(8),DIMENSION(N,N)             :: A
REAL(8)                            :: NORM,S,ALPHA,BETA,M,MM,SN
!-----------------------------------------------------------------------------
R_OLD=B
RS_OLD  = R_OLD
P_OLD   = R_OLD
PS_OLD  = RS_OLD
NORM=1.D0
DO WHILE (NORM.GT.TOL)
ITER=ITER+1
!PRINT*,ITER
!------------------------------------------------------------------------------
DO I=1,N
AP(I)=DOT_PRODUCT(A(I,:),P_OLD)
END DO
!------------------------------------------------------------------------------
M=DOT_PRODUCT(R_OLD,RS_OLD)
MM=DOT_PRODUCT(AP,PS_OLD)
ALPHA=M/MM
X=X_OLD+ALPHA*P_OLD
R=R_OLD-ALPHA*AP
!------------------------------------------------------------------------------
DO J=1,N
ATPS(J)=DOT_PRODUCT(A(:,J),PS_OLD)
END DO
!------------------------------------------------------------------------------
RS=RS_OLD-ALPHA*ATPS
S=DOT_PRODUCT(R,RS)
SN=DOT_PRODUCT(R,R)
NORM=SQRT(SN)/N
BETA=S/M
P=R+BETA*P_OLD
PS=RS+BETA*PS_OLD
PS_OLD   =PS
P_OLD    =P
X_OLD    =X
R_OLD    =R
RS_OLD   =RS
!PRINT*,NORM
CALL SUCCESSIVE_SOL ("BICG",NORM,ITER)
END DO
B=X
!------------------------------------------------------------------------------
END SUBROUTINE BICG
!******************************************************************************
!******************************************************************************
!                      GMRES_STAR Subroutine
!******************************************************************************
!******************************************************************************
SUBROUTINE GMRES_STAR (A,B)
USE BEM
IMPLICIT NONE
INTEGER                           :: I,J,K
REAL(8)                           :: NORMC,NORM,ALPHA,CR,S,SR
REAL(8),DIMENSION(N)              :: B,ZM,X_OLD,XX_OLD,R_OLD,X,R,C
REAL(8),DIMENSION(ITMAX,N)        :: CC,U
REAL(8),DIMENSION(N,N)            :: A
!------------------------------------------------------------------------------
R_OLD=B
DO I=1,ITMAX
!------------------------------------------------------------------------------
CALL BCGSTAB_INNER  (30,A,XX_OLD,R_OLD,ZM)
!CALL CGS_INNER     (30,A,XX_OLD,R_OLD,ZM)
!CALL BICG_INNER    (30,A,XX_OLD,R_OLD,ZM)
!------------------------------------------------------------------------------
DO J=1,N
C(J)=DOT_PRODUCT(A(J,:),ZM)
END DO
!------------------------------------------------------------------------------
DO K=1,I-1
ALPHA=DOT_PRODUCT(CC(K,:),C)
C=C-ALPHA*CC(K,:)
ZM=ZM-ALPHA*U(K,:)
END DO
S=DOT_PRODUCT(C,C)
NORMC=DSQRT(S)
CC(I,:)=C/NORMC
U(I,:)=ZM/NORMC
CR=DOT_PRODUCT(CC(I,:),R_OLD)
X=X_OLD+CR*U(I,:)
R=R_OLD-CR*CC(I,:)
SR=DOT_PRODUCT(R,R)
NORM=DSQRT(SR)/N
R_OLD=R
X_OLD=X
!PRINT*,I,NORM
CALL SUCCESSIVE_SOL ("GMRES",NORM,I)
IF (NORM.LT.TOL) THEN
EXIT
ELSE IF (I.EQ.ITMAX.AND.NORM.GT.TOL) THEN
PRINT*, 'convergency is not met by this number of iteration'
END IF
!------------------------------------------------------------------------------
END DO
B=X
!------------------------------------------------------------------------------
END SUBROUTINE GMRES_STAR
!******************************************************************************
!******************************************************************************
!                          BCGSTAB_INNER Subroutine
!******************************************************************************
!******************************************************************************
SUBROUTINE BCGSTAB_INNER (INNER_ITER,A,X_OLD,B,X)
USE BEM
IMPLICIT NONE
INTEGER                            :: I,J,INNER_ITER
REAL(8),DIMENSION(N)               :: X_OLD,X,P_OLD,P,R_OLD,R,AP,AQQ,RS,QQ,B
REAL(8),DIMENSION(N,N)             :: A
REAL(8)                            :: ALPHA,BETA,M,MS,MMS,MM,MN,OM
!------------------------------------------------------------------------------
R_OLD=B
!------------------------------------------------------------------------------
RS=R_OLD
P_OLD=R_OLD
DO I=1,INNER_ITER
MS=0.D0
MMS=0.D0
MN=0.D0
!------------------------------------------------------------------------------
DO J=1,N
AP(J)=DOT_PRODUCT(A(J,:),P_OLD)
END DO
!------------------------------------------------------------------------------
M=DOT_PRODUCT(R_OLD,RS)
MM=DOT_PRODUCT(AP,RS)
ALPHA=M/MM
QQ=R_OLD-ALPHA*AP
!------------------------------------------------------------------------------
DO J=1,N
AQQ(J)=DOT_PRODUCT(A(J,:),QQ)
END DO
!------------------------------------------------------------------------------
MS=DOT_PRODUCT(AQQ,QQ)
MMS=DOT_PRODUCT(AQQ,AQQ)
OM=MS/MMS
X=X_OLD+ALPHA*P_OLD+OM*QQ
R=QQ-OM*AQQ
MN=DOT_PRODUCT(R,RS)
BETA=(MN/M)*(ALPHA/OM)
P=R+BETA*(P_OLD-OM*AP)
P_OLD=P
X_OLD=X
R_OLD=R
END DO
!------------------------------------------------------------------------------
END SUBROUTINE BCGSTAB_INNER
!******************************************************************************
!******************************************************************************
!                         CGS_INNER Subroutine
!______________________________________________________________________________
SUBROUTINE CGS_INNER (INNER_ITER,A,X_OLD,B,X)
USE BEM
IMPLICIT NONE
INTEGER                 :: I,J,INNER_ITER
REAL(8),DIMENSION(N)    :: X_OLD,X,P_OLD,P,R_OLD,U,U_OLD,R,AP,AUQ,RS,Q,B
REAL(8),DIMENSION(N,N)  :: A
REAL(8)                 :: ALPHA,BETA,M,MM,MN
!------------------------------------------------------------------------------
R_OLD=B
!------------------------------------------------------------------------------
RS=R_OLD
P_OLD=R_OLD
U_OLD=R_OLD
DO I=1,INNER_ITER
!------------------------------------------------------------------------------
DO J=1,N
AP(J)=DOT_PRODUCT(A(J,:),P_OLD)
END DO
!------------------------------------------------------------------------------
M=DOT_PRODUCT(R_OLD,RS)
MM=DOT_PRODUCT(AP,RS)
ALPHA=M/MM
Q=U_OLD-ALPHA*AP
X=X_OLD+ALPHA*(U_OLD+Q)
!------------------------------------------------------------------------------
DO J=1,N
AUQ(J)=DOT_PRODUCT(A(J,:),U_OLD+Q)
END DO
!------------------------------------------------------------------------------
R=R_OLD-ALPHA*AUQ
MN=DOT_PRODUCT(R,RS)
BETA=MN/M
U=R+BETA*Q
P=U+BETA*(Q+BETA*P_OLD)
P_OLD=P
X_OLD=X
R_OLD=R
U_OLD=U
END DO
!------------------------------------------------------------------------------
END SUBROUTINE CGS_INNER
!******************************************************************************
!******************************************************************************
!                         BCG_SUBROUTINE_SAAD
!_____________________________________________________________________________
SUBROUTINE BICG_INNER (INNER_ITER,A,X_OLD,B,X)
USE BEM
IMPLICIT NONE
INTEGER                            :: I,J,K,INNER_ITER
REAL(8),DIMENSION(N)               :: X_OLD,X,P_OLD,P,R_OLD,R,B,&
                                      AP,RS_OLD,PS_OLD,PS,RS,ATPS
REAL(8),DIMENSION(N,N)             :: A
REAL(8)                            :: S,ALPHA,BETA,M,MM
!-----------------------------------------------------------------------------
R_OLD=B
RS_OLD  = R_OLD
P_OLD   = R_OLD
PS_OLD  = RS_OLD
DO I=1,INNER_ITER
!------------------------------------------------------------------------------
DO K=1,N
AP(K)=DOT_PRODUCT(A(K,:),P_OLD)
END DO
!------------------------------------------------------------------------------
M=DOT_PRODUCT(R_OLD,RS_OLD)
MM=DOT_PRODUCT(AP,PS_OLD)
ALPHA=M/MM
X=X_OLD+ALPHA*P_OLD
R=R_OLD-ALPHA*AP
!------------------------------------------------------------------------------
DO J=1,N
ATPS(J)=DOT_PRODUCT(A(:,J),PS_OLD)
END DO
!------------------------------------------------------------------------------
RS=RS_OLD-ALPHA*ATPS
S=DOT_PRODUCT(R,RS)
BETA=S/M
P=R+BETA*P_OLD
PS=RS+BETA*PS_OLD
PS_OLD   =PS
P_OLD    =P
X_OLD    =X
R_OLD    =R
RS_OLD   =RS
END DO
!------------------------------------------------------------------------------
END SUBROUTINE BICG_INNER
!******************************************************************************
!******************************************************************************
!                    INITIAL_RESIDUAL_SUBROUTINE
!______________________________________________________________________________
SUBROUTINE INITIAL_RESIDUAL (A,RHS,X0,R0)
USE BEM
IMPLICIT NONE
INTEGER                     :: I
REAL(8),DIMENSION(N)        :: RHS,X0,R0,Y0
REAL(8),DIMENSION(N,N)      :: A
!------------------------------------------------------------------------------
DO I=1,N
Y0(I)=DOT_PRODUCT(A(I,:),X0)
END DO
!------------------------------------------------------------------------------
R0=RHS-Y0
!------------------------------------------------------------------------------
END SUBROUTINE INITIAL_RESIDUAL
!******************************************************************************

