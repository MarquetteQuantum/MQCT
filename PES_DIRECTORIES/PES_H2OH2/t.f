      function tnormed(p1,q1,p2,p,th,ph,thp,php)
c
c la fonction t de PMG94, cette fois normee a 1
c les angles sont en degres
c
c attention, coquille formule 8 dans PMG94,
c lire 
c
c                            (1 + \delta_q1_0)^{-1} 
c
c et non (...)^{-2} !!!
c
      implicit none
      integer p1,q1,p2,p
      real*8 tnormed, t, th,ph,thp,php, xn, Pi, dq1, one, two
      parameter (Pi=3.14159265358979324d0, one=1.0d0, two=2.0d0)

!      write(6,*) p1,q1,p2,p,th,ph,thp,php

      if (q1.eq.0) then
         dq1=1
      else
         dq1=0
      endif

      xn = two / (one+dq1) / (2*p1+1)

      tnormed = t(p1,q1,p2,p,th,ph,thp,php) / sqrt(xn)
      return
      end

      function tprime(p1,q1,p2,p,th,ph,thp,php)
c
c la fonction t de PMG94, cette fois normee a 4Pi
c les angles sont en degres
c
c attention, coquille formule 8 dans PMG94,
c lire 
c
c                            (1 + \delta_q1_0)^{-1} 
c
c et non (...)^{-2} !!!
c
      implicit none
      integer p1,q1,p2,p
      real*8 tprime, t, th,ph,thp,php, xn, Pi, dq1, one, two
      parameter (Pi=3.14159265358979324d0, one=1.0d0, two=2.0d0)

!      write(6,*) p1,q1,p2,p,th,ph,thp,php

      if (q1.eq.0) then
         dq1=1
      else
         dq1=0
      endif

      xn = two / (one+dq1) / (2*p1+1)

      tprime = t(p1,q1,p2,p,th,ph,thp,php)
     $     * (4*Pi) / sqrt(xn)
      return
      end

      function t(p1,q1,p2,p,th,ph,thp,php)
      implicit none
      integer p1,q1,p2,p, r1,r2,r
      real*8 t, th,ph,thp,php, d1,d2,fact, PLMx,
     $     xp1,xp2,xp,xr1,xr2,xr, xf3j, thrj, pi,
     $     PARITC
      parameter (Pi=3.14159265358979324d0)

      t=0
      do r1=-p1,p1
         do r2=-p2,p2

            r=-r1-r2
            if (abs(r).gt.p) cycle
               
            if (q1.eq.r1) then
               d1=1
            else
               d1=0
            endif
            if (q1+r1.eq.0) then
               d2=1
            else
               d2=0
            endif

            xp1=p1; xp2=p2; xp=p
            xr1=r1; xr2=r2; xr=r

! obeying painfully to SG conventions...

            t = t + thrj(xp1,xp2,xp,xr1,xr2,xr)
     $           * PLMx(p2,r2,cos(thp*pi/180.0d0))
     $           * (d1+PARITC(p1+q1+p2+p)*d2)
     $           * PLMx(p,r,cos(th*pi/180.0d0))
     $           * cos((r2*php+r*ph)*pi/180.0d0) / (2*PI)

         enddo
      enddo
      if (q1.eq.0) t=t*0.5
      end

      function PLMx(L,M,C)
      implicit none
      integer L, M
      real*8 PLMx, PLM, C, PARITC
      PLMx=PLM(L,M,C)
      if (M.lt.0) then
         PLMx=PLMx*PARITC(M)
      endif
      end

      FUNCTION PLM(LIN,MIN,COSTH)
C
C     COMPUTES NORMALIZED ASSOC. LEGENDRE POLYNOMIALS BY RECURSION.
C     THE VALUES RETURNED ARE NORMALIZED FOR INTEGRATION OVER X
C     (I.E. INTEGRATION OVER COS THETA BUT NOT PHI).
C     NOTE THAT THE NORMALIZATION GIVES
C           PLM(L,0,1)=SQRT(L+0.5)
C           PLM(L,0,X)=SQRT(L+0.5) P(L,X)
C     FOR M.NE.0, THE VALUE RETURNED DIFFERS FROM THE USUAL
C           DEFINITION OF THE ASSOCIATED LEGENDRE POLYNOMIAL
C           (E.G. EDMONDS PAGES 23-24)
C           BY A FACTOR OF (-1)**M*SQRT(L+0.5)*SQRT((L-M)!/(L+M)!)
C     THUS THE SPHERICAL HARMONICS ARE
C          CLM = PLM * EXP(I*M*PHI) / SQRT(L+0.5)
C          YLM = PLM * EXP(I*M*PHI) / SQRT(2*PI)
C     THIS ROUTINE ALWAYS RETURNS THE VALUE FOR ABS(MIN); NOTE THAT
C          FOR MIN.LT.0 THIS VALUE SHOULD BE MULTIPLIED BY PARITC(MIN)
C
C     FUNCTION PM1(LIN,MIN,COSTH)
C     This routine appears to be much more stable for large l, m than
C       the routine from Nerf/ modified according to R.T Pack
C     It was obtained:
C     From: Marie-Lise Dubernet <mld@ipp-garching.mpg.de>
C     Date: Mon, 19 Jun 1995 12:48:11 +0200 (MET DST)
C     Some mods 27-28 June 95 by SG for speed and to accord w/ MOLSCAT 
C     Bugs fixed 21 Sept 95 (SG)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     CHECK FOR ABS(COSTH).LE.1.D0 ... 
      IF (ABS(COSTH).GT.1.D0) THEN
        WRITE(6,*) ' *** ILLEGAL ARGUMENT TO PLM. X =',COSTH
        STOP
      ENDIF
C     SAVE ARGUMENTS IN LOCAL VARIABLES
      L=LIN
      M=ABS(MIN)
      X=COSTH
C
C  IF M>L PLM=0 !
      IF(M.GT.L) THEN
        PLM=0.D0
        RETURN
      ENDIF
      LMAX=L
C
      IF (M.GT.0) GO TO 5
C  HERE FOR REGULAR LEGENDRE POLYNOMIALS
      PLM=1.D0
      PM2=0.D0
      XL=0.D0
      DO 2 L=1,LMAX
      XL=XL+1.D0
      PP=((2.D0*XL-1.D0)*X*PLM-(XL-1.D0)*PM2)/XL
      PM2=PLM
    2 PLM=PP
      GO TO 9000
C
C  HERE FOR ALEXANDER-LEGENDRE POLYNOMIALS
C
    5 IMAX=2*M
      RAT=1.D0
      AI=0.D0
      DO 6 I=2,IMAX,2
      AI=AI+2.D0
    6 RAT=RAT*((AI-1.D0)/AI)
C     Y=SIN(THETA)
      Y=SQRT(1.D0-X*X)
      PLM=SQRT(RAT)*(Y**M)
      PM2=0.D0
      LOW=M+1
      XL=LOW-1
      DO 10 L=LOW,LMAX
      XL=XL+1.D0
      AL=DBLE((L+M)*(L-M))
      AL=1.D0/AL
      AL2=(DBLE((L+M-1)*(L-M-1)))*AL
      AL=SQRT(AL)
      AL2=SQRT(AL2)
      PP=(2.D0*XL-1.D0)*X*PLM*AL-PM2*AL2
      PM2=PLM
   10 PLM=PP
      PLM=PLM*PARITC(MIN)
C
C     CONVERT TO MOLSCAT'S IDIOSYNCRATIC NORMALIZATION
9000  PLM=PLM*SQRT(XL+0.5D0)
      RETURN
      END

      FUNCTION THRJ(F1,F2,F3,G1,G2,G3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SMALL CHANGES 31 JUL 95 (SG)
      SAVE MUNG,X,Y
      PARAMETER (MXIX=302)
      DIMENSION X(MXIX),Y(MXIX)
      DATA MUNG/0/
      IF (MUNG.EQ.21) GO TO 69
      MUNG = 21
      X(1) = 0.D0
      DO 100 I = 1, MXIX-1
      A = I
      X(I+1) = LOG(A) +X(I)
      Y(I+1) = LOG(A)
  100 CONTINUE
   69 IF(F1-ABS(G1)) 1,13,13
   13 IF(F2-ABS(G2))1,14,14
   14 IF(F3-ABS(G3))1,15,15
   15 SUM=F1+F2+F3
      NSUM=SUM+.001D0
      IF(SUM-NSUM)2,2,1
    1 THRJ=0.D0
      RETURN
    2 IF(ABS(G1+G2+G3)-1.D-08)3,3,1
    3 IF(F1+F2-F3)1,4,4
    4 IF(F1+F3-F2)1,5,5
    5 IF(F2+F3-F1)1,6,6
    6 J1=2.D0*F3+2.001D0
      J2=F1+F2-F3+1.001D0
      J3=F1-F2+F3+1.001D0
      J4=-F1+F2+F3+1.001D0
      J5=F1+F2+F3+2.001D0
      J6=F1+G1+1.001D0
      J7=F1-G1+1.001D0
      J8=F2+G2+1.001D0
      J9=F2-G2+1.001D0
      J10=F3+G3+1.001D0
      J11=F3-G3+1.001D0
      IF(J5.GT.MXIX) THEN
        WRITE(6,601) J5,MXIX
  601   FORMAT(' *** DIMENSION ERROR IN THRJ - INDEX.GT.MXIX',2I5)
        STOP
      ENDIF
      R=0.5D0*(Y(J1)+X(J2)+X(J3)+X(J4)-X(J5)
     1+X(J6)+X(J7)+X(J8)+X(J9)+X(J10)+X(J11))
      SUM=0.D0
      F=-1
      KZ=-1
    7 KZ=KZ+1
      F=-F
      J1=KZ+1
      J2=F1+F2-F3-KZ+1.001D0
      IF(J2)20,20,8
    8 J3=F1-G1-KZ+1.001D0
      IF(J3)20,20,9
    9 J4=F2+G2-KZ+1.001D0
      IF(J4)20,20,10
   10 J5=F3-F2+G1+KZ+1.001D0
      IF(J5)7,7,11
   11 J6=F3-F1-G2+KZ+1.001D0
      IF(J6)7,7,12
   12 JMAX=MAX(J1,J2,J3,J4,J5,J6)
      IF(JMAX.GT.MXIX) THEN
        WRITE(6,601) JMAX,MXIX
        STOP
      ENDIF
      S=-(X(J1)+X(J2)+X(J3)+X(J4)+X(J5)+X(J6))
      SUM=SUM+F*EXP(R+S)
      GO TO 7
   20 INT=ABS(F1-F2-G3)+0.0001D0
      VAL=((-1.D0)**INT)*SUM/SQRT(2.D0*F3+1.D0)
      IF(ABS(VAL).LE.1.D-6) VAL=0.D0
      THRJ=VAL
      RETURN
      END

      FUNCTION PARITC(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARITC=1.D0
      IF((I/2)*2-I.NE.0) PARITC=-1.D0
      RETURN
      END
