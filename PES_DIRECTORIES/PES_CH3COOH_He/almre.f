c This is the real part of  Wormer's ALM function. Intermolecular
c vector R assumed along the z axis. Seems much faster than the
c original general one. No spherical harmonics called.
c
C***********************************************************************
      FUNCTION almre(LA,KA,LB,KB,L,AOMEGA,BOMEGA,COMEGA)
C
C     ANGULAR FUNCTION DEPENDING ON TWO SETS OF THREE EULER ANGLES
C     (AOMEGA AND BOMEGA) TWO POLAR ANGLES (COMEGA) are not needed.
c     It is assumed that theta = phi = 0.
C
C     FUNCTION IS A COUPLING OF TWO WIGNER D-MATRICES, DEPENDING
C     ON THE EULER ANGLES OF A AND B,
C     AND A SPHERICAL HARMONIC IN RACAH'S NORMALIZATION, equal to 1
C     if THE INTERMOLECULAR VECTOR is along z axis,
C     TO A SCALAR FUNCTION.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (LMAX=50,THRESH=1.D-10)
      DIMENSION WORK(0:2*LMAX,0:2*LMAX)
      DIMENSION AOMEGA(3),BOMEGA(3),COMEGA(2)
      dimension DLKA(0:2*LMAX),DLKB(0:2*LMAX)

      IF ( (LA.GT.LMAX) .OR. (LB.GT.LMAX) .OR. (L.GT.LMAX) ) THEN
         WRITE(*,'(''0L-VALUES LARGER THAN '',I3,''. INCREASE LMAX IN'',
     1             '' FUNCTION ALM'')') LMAX
         STOP 16
      ENDIF

C COLUMN KA OF  WIGNER D-MATRIX
      CALL DMAT(DBLE(LA),AOMEGA(2),LMAX,WORK)
c      Z = DCMPLX(0.D0, DBLE(KA)*AOMEGA(3))
c      EXPK = CDEXP(Z)
c      DMA = DBLE(-LA-1)
      DO 10 MA = -LA,LA
c       DMA = DMA + 1.D0
c       Z = DCMPLX(0.D0,DMA*AOMEGA(1))
       DLKA(LA+MA) = WORK(LA+MA,LA+KA)
   10 CONTINUE

C COLUMN KB OF  WIGNER D-MATRIX
      CALL DMAT(DBLE(LB),BOMEGA(2),LMAX,WORK)
c      Z = DCMPLX(0.D0, DBLE(KB)*BOMEGA(3))
c      EXPK = CDEXP(Z)
c      DMB = DBLE(-LB-1)
      DO 20 MB = -LB,LB
c       DMB = DMB + 1.D0
c       Z = DCMPLX(0.D0,DMB*BOMEGA(1))
c       DLKB(LB+MB) = CDEXP(Z) * DCMPLX(WORK(LB+MB,LB+KB),0.D0) * EXPK
       DLKB(LB+MB) = WORK(LB+MB,LB+KB)
   20 CONTINUE

C CONTRACT, RESULT IN Z.
      CALL SUM(LA,ka,LB,kb,L,DLKA,DLKB,aomega,bomega,Z)
      IF ( MOD(LA+LB+L,2) .EQ. 0 ) THEN
         almre = Z
      ELSE
         almre = -Z
      ENDIF

      END
C***********************************************************************
      SUBROUTINE SUM(LA,ka,LB,kb,L,DLKA,DLKB,aomega,bomega,z)

C     COUPLE DLKA,DLKB AND CLM  TO SCALAR FUNCTION.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER( ZERO=0.D0, ONE=1.0D0, JMAX=50 )
      dimension DLKA(-LA:LA),DLKB(-LB:LB),aomega(3),bomega(3)
      data nzer /0/

      z = 0.D0
      lmm = min(la,lb)
      DO 30 MA=-lmm,lmm
      mb = -ma
      term = ka*aomega(3) + kb*bomega(3) + ma*(aomega(1) - bomega(1))
      term = dcos(term)
          z = z + THREEJ(LA,LB,L, MA,mb,nzer)
     1                    *DLKA(MA)*DLKB(MB)*term
   30 CONTINUE
      END
C***********************************************************************
      SUBROUTINE DMAT(DJ,BETA,JMAX,D)
C
C     -- COMPUTE WIGNER D-MATRIX ("SMALL D") FOR ANGLE BETA.
C
C     INPUT PARAMETERS:
C        DJ:   DOUBLE PRECISION J-QUANTUM NUMBER, INTEGRAL OR
C              HALF-INTEGRAL, DJ <= JMAX
C        BETA: ROTATION ANGLE (RADIANS), 0 <= BETA <= PI.
C        JMAX: GIVES DIMENSION OF D.
C
C     OUTPUT PARAMETER:
C        D:    SQUARE MATRIX OF  FORTRAN DIMENSIONS 2*JMAX+1 AND
C              FILLED IN FROM (0,0) TO (2*DJ,2*DJ)
C
C     NOTES:
C     ----
C     1. THE ELEMENT D(M1,M2) (-J <= M1,M2 <= +J) IS ACCESSED
C        BY D(M1+DJ,M2+DJ).
C     2. THE ALGORITHM IS LEAST STABLE FOR BETA APPROXIMATELY PI/2.
C        IF HIGH PRECISION IS REQUIRED COMPUTE THE MATRIX FOR PI/4
C        SQUARE  THE MATRIX
C     3. THE SUBROUTINE CHECKS ON STABILITY AND WRITES WARNING
C        MESSAGES TO UNIT IWARN (DEFAULT 6)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical lbeta
      PARAMETER (PI = 3.14159 26535 89793 2D0,TOLER=1.D-15,
     1           MAXT = 50)
      DIMENSION D(0:2*JMAX,0:2*JMAX),T(0:2*MAXT)
      DATA IWARN/6/

      TWOJ  = 2.D0*DJ
      ITWOJ = NINT(TWOJ)
      IF ( ABS(DBLE(ITWOJ)-TWOJ) .GT. TOLER ) THEN
         WRITE(*,'(''0*** ERROR: DJ = '',D23.15,'' WHICH IS NOT'',
     1             '' INTEGRAL OR HALF-INTEGRAL'')') DJ
         STOP 16
      ENDIF
C
      IF (JMAX.GT.MAXT) THEN
         WRITE(*,'( ''0*** ERROR: J TOO LARGE'',
     1              '' INCREASE TMAX IN SUBROUTINE DMAT'')')
         STOP 16
      ENDIF
C
      IF (  (BETA.LT.-TOLER) .OR. (BETA.GT.PI+TOLER) ) THEN
         WRITE(*,'(''0*** ERROR: BETA = '',D23.15,'' DEGREES, WHICH'',
     1             '' IS OUT OF RANGE'')') BETA*180.D0/PI
         STOP 16
      ENDIF
C
      IF ( ABS(DJ) .LT. TOLER ) THEN
         D(0,0) = 1.D0
         RETURN
      ENDIF
C
      IF ( ABS(BETA) .LT. TOLER ) THEN
C        D IS UNIT MATRIX
         DO 20 K=0,ITWOJ
         D(K,K) = 1.D0
         DO 10 L=0,K-1
         D(K,L) = 0.D0
         D(L,K) = 0.D0
   10    CONTINUE
   20    CONTINUE
         RETURN
      ENDIF
C
      IF ( ABS(BETA-PI) .LT. TOLER ) THEN
C        D IS DIAGONAL (ALONG THE MINOR DIAGONAL)
         DO 40 K=0,ITWOJ
         DO 30 L=0,ITWOJ
         D(K,L) = 0.D0
   30    CONTINUE
   40    CONTINUE
         PHASE = 1.D0
         DO 50 K=0,ITWOJ
         D(K,ITWOJ-K) = PHASE
         PHASE = - PHASE
   50    CONTINUE
         RETURN
      ENDIF
C
C     FOR  0 <= BETA <= PI/2  RECURSIVE CALCULATION OF
C     UNDER TRIANGLE IS FAIRLY STABLE IF WE FIRST REMOVE
C     POWERS OF TAN(BETA/2)
C

      LBETA = .FALSE.
      IF ( BETA .GT. 0.5D0*PI ) THEN
          BETA = PI - BETA
          LBETA = .TRUE.
      ENDIF

C     -  FIRST (= -J) COLUMN BY RECURSION RELATION, WHICH
C        CAN EASILY BE DERIVED FROM THE GENERAL EXPRESSION (3.65)
C        OF BIEDENHARN & LOUCK.
C     -  SECOND COLUMN (-J+1) FROM FIRST COLUMN BY
C        EQ. (3.84) OF B. & L.
C
      STWOJ = SQRT(TWOJ)
      TGHB  = TAN(0.5D0*BETA)
      CB    = COS(BETA)
      SB    = SIN(BETA)
      DENOM = 1.D0/(STWOJ*SB)
      D(0,0) = (COS(0.5D0*BETA))**ITWOJ
      DM = -DJ
C     (DM=M, M IS ROW INDEX, M=-J,...,J)
      DO 60 K=1,ITWOJ
      DM =  DM + 1.D0
      D(K,0) = -D(K-1,0)*SQRT( (DJ-DM+1.D0)/(DJ+DM) )
      D(K,1) = -D(K,0)*(TWOJ*CB+2.D0*DM)*DENOM*TGHB
   60 CONTINUE
C
C     THE REMAINING COLUMNS (M'=-J+2,...,M) FOLLOW BY EQ. (3.84)
C     OF B. & L.
C
      DM  = 1.D0-DJ
C     (DM=M)
      DO 80 L=2,ITWOJ
      DM  = DM  + 1.D0
      DMP = 1.D0-DJ
C     (DMP=M', M' IS COLUMN INDEX)
      DO 70 K=2,L
      DMP = DMP + 1.D0
      D(L,K) = 1.D0/SQRT( (DJ-DMP+1.D0)*(DJ+DMP) )*
     1  (  2.D0*D(L,K-1)*((DMP-1.D0)*CB-DM)*TGHB/SB
     2         -D(L,K-2)*SQRT((DJ+DMP-1.D0)*(DJ-DMP+2.D0))*TGHB*TGHB  )
   70 CONTINUE
   80 CONTINUE

C     SET UP POWERS OF TGHB
      T(0) = 1.D0
      DO 90 K=1,ITWOJ
      T(K) = T(K-1)*TGHB
   90 CONTINUE

C     TRANSPOSE AND MULTIPLY BY POWERS OF TAN(BETA/2)

      SIGN1=1.D0
      DO 110 K=1,ITWOJ
      SIGN1 = -SIGN1
      SIGN2 = -SIGN1
      DO 100  L=0,K-1
      D(K,L) = T(K-L)*D(K,L)
      SIGN2 = -SIGN2
      D(L,K) = SIGN2*D(K,L)
  100 CONTINUE
  110 CONTINUE

C
C     CHECK ON STABILITY AND CORRECTNESS
C
      IF ( ABS (D(0,0)-D(ITWOJ,ITWOJ)) .GT. TOLER*1.D+3 ) THEN
         IF (LBETA) BETA = PI -BETA
         WRITE(IWARN,'(''0*** WARNING D(-J,-J) .NE. D(J,J) ''/
     1                 '' BETA     = '',D23.15/
     2                 '' D(-J,-J) = '',D23.15/
     3                 '' D( J, J) = '',D23.15)' )
     4             BETA*180.D0/PI, D(0,0),D(ITWOJ,ITWOJ)
         IF (LBETA) BETA = PI -BETA
      ENDIF
      IF (LBETA) THEN
C        BETA WAS  INITIALLY > PI/2, SET BACK.
         BETA = PI - BETA
         IF ( ABS(DBLE(NINT(DJ))-DJ) .GT. TOLER ) THEN
C           HALF INTEGRAL J, PHASE = (-1)**2J = -1
            PHASE= -1.D0
         ELSE
C          INTEGRAL J,     PHASE = (-1)**2J = +1
           PHASE =  1.D0
         ENDIF
         DO 130 K=0,ITWOJ
         DO 120 L=0,ITWOJ/2
         TEMP = D(ITWOJ-L,K)
         D(ITWOJ-L,K) = PHASE*D(L,K)
         D(L,K) = PHASE*TEMP
  120    CONTINUE
         PHASE = - PHASE
  130    CONTINUE
      ENDIF
      END
C***********************************************************************
      FUNCTION THREEJ(J1,J2,J3,M1,M2,M3)

C     WIGNER 3-J ONLY FOR INTEGER QUANTUM NUMBERS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER( ZERO=0.D0, ONE=1.0D0, JMAX=50 )
      DOUBLE PRECISION  BINOM(0:JMAX,0:JMAX)
      LOGICAL FIRST
      DATA FIRST  /.TRUE./

      THREEJ = ZERO

      IF ( (J1 +J2 +J3 + 1) .GT. JMAX ) THEN
         WRITE(*,'( ''0SUM OF THE J-VALUES LARGER THAN '',I3)') JMAX-1
         WRITE(*,'( ''0INCREASE JMAX IN FUNCTION THREEJ'')')
         GOTO 999
      ENDIF

      IF (FIRST) THEN
C FILL OUT BINOM BY PASCAL'S TRIANGLE RELATION
        BINOM(0,0) = 1.D0
        DO 20 I=1,JMAX
         BINOM(I,0) = 1.D0
         BINOM(I,I) = 1.D0
         DO 10 J=1,I-1
           BINOM(I,J) = BINOM(I-1,J-1) + BINOM(I-1,J)
   10    CONTINUE
   20   CONTINUE
        FIRST = .FALSE.
      ENDIF

      IF ((M1+M2+M3).NE.0) GOTO 999
      I1 = -J1+J2+J3
      IF (I1.LT.0) GOTO 999
      I2 =  J1-J2+J3
      IF (I2.LT.0) GOTO 999
      I3 =  J1+J2-J3
      IF (I3.LT.0) GOTO 999
      K1 =  J1+M1
      IF (K1.LT.0) GOTO 999
      K2 =  J2+M2
      IF (K2.LT.0) GOTO 999
      K3 =  J3+M3
      IF (K3.LT.0) GOTO 999
      L1 =  J1-M1
      IF (L1.LT.0) GOTO 999
      L2 =  J2-M2
      IF (L2.LT.0) GOTO 999
      L3 =  J3-M3
      IF (L3.LT.0) GOTO 999
      N1 = -J1-M2+J3
      N2 =  M1-J2+J3
      N3 =  J1-J2+M3
      IMIN = MAX0(-N1,-N2,0)
      IMAX = MIN0(L1,K2,I3)
      IF (IMIN .GT. IMAX) GOTO 999
      SIGN = - ONE
      DO 30 I=IMIN,IMAX
      SIGN = -SIGN
      THREEJ = THREEJ + SIGN*BINOM(I1,N1+I)*BINOM(I2,N2+I)*BINOM(I3,I)
   30 CONTINUE
      THREEJ = THREEJ * DSQRT(BINOM(J2+J2,I3)*BINOM(J1+J1,I2)
     1       / (BINOM(J1+J2+J3+1,I3)*DBLE(J3+J3+1)
     2       * BINOM(J1+J1,L1)*BINOM(J2+J2,L2)*BINOM(J3+J3,L3)))
      IF (MOD(N3+IMIN,2) .NE. 0) THREEJ = - THREEJ
  999 RETURN
      END
C***********************************************************************
