      SUBROUTINE USER_DEFINED_PES(V,R,rvib,rvib2,alpha,beta,gamma,
     & 								aalpha,bbeta,ggamma)
!     INPUT:  R - distance betweenn COMs, rvib - vibrational coordinate
! 		      alpha,beta,gamma - Euler's angles of the first molecule
!   		  aalpha,bbeta, ggamma - Euler's angles of the second molecule
!     OUTPUT: V - value of the potential
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 V,R,alpha,beta,gamma,aalpha,bbeta,ggamma,t,rvib,rvib2
      REAL*8 R1,COG, phi, cosa, cosb
	  real*8, parameter :: PI = 4.d0*datan(1.d0)
	  
	  cosa = dcos(beta)
	  cosb = dcos(bbeta)
	  phi = PI - gamma
	  V = POTFUN(R,cosa,cosb,phi)
	  
      END SUBROUTINE USER_DEFINED_PES


!----------!
! OPTION 2 !
!----------! 
!	  USE KEYWORD "EXPANSION=YES" TO INITIATE THIS OPTION
      SUBROUTINE USER_DEFINED_TERMS(T,I,R)
!     THIS SUBROTUNE COMPUTES RADIAL COEFFICENTS OF THE PES EXPANSION AT A GIVEN DISTANCE
!     INPUT:  R - distance between COMs of particles, I - TERM NUMBER
!     OUTPUT: T - value of coefficent 	  
      IMPLICIT NONE	  
      REAL*8 T,R
      INTEGER I
!     USER MUST INSERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:
      STOP "ERROR: USER_DEFINED_TERMS IS NOT SUPPLIED"
	  END SUBROUTINE USER_DEFINED_TERMS
!----------!
! OPTION 3 !
!----------! 
!	  USE KEYWORDS "EXPANSION=YES, TERMS_FILE=YES" TO INITIATE THIS OPTION
! 	  SIMILAR TO OPTION 2, BUT NO SUBROUTINE IS REQUIRED
!     USER SHOULD PROVIDE THE FILE EXPAN_PES_TERMS.DAT 
!     IN THE MAIN PROGRAM DIRECTORY CONTAINING THE COEFFICEINS 
!     OF POTENTIAL EXPANSION PRECOMPUTED EXTERNALLY.
! 	  SEE EXAMPLE FILES SUPPLIED WITH THE CODE.
!----------!
! OPTION 4 !
!----------! 
!	  USE KEYWORDS "EXPANSION=YES, TERMS_ONFLY=YES" TO INITIATE THIS OPTION
      SUBROUTINE USER_DEFINED_COEFFS(T,DTDR,I,R) 
!     THIS SUBROUTINE COMPUTES RADIAL COEFFICENTS OF THE PES EXPANSION 
!     AND THEIR DERIVATIVES AT A GIVEN DISTANCE R
      IMPLICIT NONE
!     INPUT : R - distance between COMs of particles, I - TERM NUMBER
!     OUTPUT: T - value of coefficent, DTDR - its radial derivative 	  
      REAL*8 T,R,DTDR 
      INTEGER I
!     USER MUST INCERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE IS SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:	
      STOP "ERROR: USER_DEFINED_COEFFS IS NOT SUPPLIED"
      END SUBROUTINE USER_DEFINED_COEFFS 


        FUNCTION POTFUN(r,xcosa,xcosb,xphi)
c r here is in bohr
c POTFUN returns energy in Hartree
C     Program to compute point on fitted CCSD(T)-SAPT PES of CO-CO
C     Units: bohr, degrees, hartree
C
C     Reference:
C     G. W. M. Vissers, P. E. S. Wormer and A. van der Avoird
C     Phys. Chem. Chem. Phys., 2003, 5,  4767

      implicit real*8 (a-h,o-z)
      integer outfile
      parameter( maxang  = 91,
     .           maxpow  = 4,
     .           maxLX   = 5,
     .           outfile = 10,
     .           inCSR   = 11,
     .           inCLR   = 12,
     .           inalpha = 13,
     .           inLLL   = 14)
      real*8  CSR(0:maxpow, maxang),
     .        CLR(2, maxang),
     .        alpha(maxang),
     .        ALM(maxang)
      integer LLL(3,maxang)
      logical first/.true./

      pi = dacos(-1.d0)
      thetaa=dacos(xcosa)
      thetab=dacos(xcosb)
      phi=xphi

C--Get fit parameters from files 11 - 14
!      if (first) then
         open(unit=inCSR,   file='./CSR.dat')
         open(unit=inCLR,   file='./CLR.dat')
         open(unit=inalpha, file='./alpha.dat')
         open(unit=inLLL,   file='./LLL.dat')
         Call getparam(CSR, CLR, alpha, LLL, maxang, maxpow)
         first = .false.
		 close(11)
		 close(12)
		 close(13)
		 close(14)
!      endif

C--Compute necessary angular functions
            Call compALM(ALM, maxLX, maxang, LLL, thetaa, thetab, phi)

C--Compute energy
            ener = energy(ALM, R, alpha, CSR, CLR, maxpow, maxang)
            potfun=ener
            return
      end

      REAL*8 FUNCTION ENERGY(ALM, R, alpha, CSR, CLR, maxpow, maxang)
      implicit real*8 (a-h,o-z)
      real*8  CSR(0:maxpow, maxang),
     .        CLR(2, maxang),
     .        alpha(maxang),
     .        ALM(maxang)

      ener = 0.d0

      do l = 1,maxang
          ESR  = CSR(0,l)
          RR   = 1
          do i = 1,maxpow
             RR  = RR*R
             ESR = ESR + CSR(i,l)*RR
          enddo
          ESR  = exp(-alpha(l)*R)*ESR

          nexp = -idnint(CLR(2,l))
          dmp  = TTdamp(R, alpha(l), nexp)
          ELR  = CLR(1,l)*(R**CLR(2,l))
          ener = ener + (ESR  + dmp*ELR)*ALM(l)
      enddo
      energy = ener

      end
      SUBROUTINE GETPARAM(CSR, CLR, alpha,LLL, maxang, maxpow)
      implicit real*8 (a-h,o-z)
      real*8  CSR(0:maxpow, maxang),
     .        CLR(2, maxang),
     .        alpha(maxang)
      integer LLL(3,maxang)
      logical echo


      do i=1,maxang
         read(11,*) (CSR(j,i),j=0,maxpow)
      enddo

      do i=1,maxang
         read(12,*) (CLR(j,i),j=1,2)
      enddo

      do i=1,maxang
         read(13,*) alpha(i)
      enddo

      do i=1,maxang
         read(14,*) (LLL(j,i),j=1,3)
      enddo

      echo = .false.
      if (echo) then
         write(*,*) ' Short range coefficients'
         do i=1,maxang
            write(*,*) (CSR(j,i),j=0,maxpow)
         enddo

         write(*,*) ' Long range coefficients'
         do i=1,maxang
            write(*,*) (CLR(j,i),j=1,2)
         enddo

         write(*,*) ' alphas'
         do i=1,maxang
            write(*,*) alpha(i)
         enddo

         write(*,*) ' L quantum numbers'
         do i=1,maxang
            write(*,*) (LLL(j,i),j=1,3)
         enddo
      endif

      end
      SUBROUTINE COMPALM(ALM, maxLX, maxang, LLL, thetaA, thetaB, phi)

C     Real invariant angular functions for two Sigma-state diatomics.
C     Connecting vector R is along z-axis.  Angles in radians.
C
C     Use of subroutine PLMS:
C     P(n,m;x) = sqrt(2*(n-m)!/(n+m)!) (1-x^2)^(m/2)*(d/dx)^m{P(n,x)}
C     which is normalized as:
C     int_{th=0..pi,phi=0..2*pi} [P(n,m;cos(th))*cos(m*phi)]^2
C                      * sin(th)dth dphi  = 4*pi/(2*n+1)
C
C     Norm ALMr:
C           4*pi/( (2*LA+1)*(2*LB+1)*(2*L+1) ) = int (ALMr)^2 dv
C               dv = sin(thetaA)sin(thetaB) dthetaA dthetaB dphi
C               0 < thetaA, thetaB, phi < pi
C     The ALMr are orthogonal.
C     This Fortran function gives identical results as PESW's Matlab
C     function ALMr.

      implicit real*8 (a-h,o-z)
      parameter (nA = 150,
     .           nB = 150)
      real*8 spA( (nA+1)*(nA+2)/2 ),
     .       spB( (nB+1)*(nB+2)/2 ),
     .       ALM( maxang )
      integer LLL(3,maxang)


      if ( maxLX .gt. nA)  then
         write(*,*) ' nA or nB not large enough, adapt and recompile'
         stop 16
      endif

C--Get phaseless, semi-normalized, associated Legendre functions:
      Call plm_schm(spA, nA, maxLX, thetaA)
      Call plm_schm(spB, nB, maxLX, thetaB)

      do l = 1,maxang
         LA    = LLL(1,l)
         LB    = LLL(2,l)
         Lsum  = LLL(3,l)
         tj    = threejp(LA, LB, Lsum, 0, 0,0)
         lmA   = LA*(LA+1)/2+1
         lmB   = LB*(LB+1)/2+1
         ALMr  = tj*spA(lmA)*spB(lmB)
      if ( mod(LA+LB+Lsum,2) .ne. 0 ) then
         write(*,*) 'LA+LB+L not even',LA,LB,Lsum
         stop 16
      endif


         mph = 1
         do  m=1,min(LA,LB)
            mph  = - mph
            lmA  = lmA + 1
            lmB  = lmB + 1
            tj   = threejp(LA, LB, Lsum, m, -m, 0)
            ALMr = ALMr + mph*tj*spA(lmA)*spB(lmB)*cos(dble(m)*phi)
         enddo
         ALM(l)  = ALMr
      enddo

      end
      SUBROUTINE PLM_SCHM(PLM, N, LMAX, THETA)

C     Compute Schmidt semi-normalized associated Legendre polynomials:
C     P(l,m;x) = sqrt(2*(l-m)!/(l+m)!)
C                * (1-x^2)^(m/2) * (d/dx)^m { P(l,x) } for  m > 0
C     P(l,0;x) = P(l,0)                                for  m = 0.

C     Notes:
C   1.   The associated Legendre P(L,M;X) is accessed via
C        PLM(L*(L+1)/2+M+1).
C   2.   Results are identical to Matlab's function LEGENDRE with
C        option 'SCH' (Schmidt semi-normalization).

C     Author: P.E.S. Wormer (2003)

      implicit real*8 (a-h,o-z)
      dimension plm( (n+1)*(n+2)/2 )

C--Compute associated unnormalized Legendre functions
      Call assleg(plm, n, dcos(theta))

C--Renormalization loop
      plm(1) = plm(1)                  ! l=m=0
      if ( n .eq. 0) return

      kp    = 1
      dl    = 0.d0
      sqrt2 = sqrt(2.d0)

      do  l=1,n
         kp      = kp + 1
         dl      = dl + 1.d0
         fac     = sqrt2
         plm(kp) = plm(kp)             ! l>0, m = 0
         dlm     = dl
         dlm1    = dl + 1.d0

         do m=1,l
            kp = kp + 1                ! kp   = l(l+1)/2 + 1 + m
            dlm  = dlm  + 1.d0         ! dlm  = l+m
            dlm1 = dlm1 - 1.d0         ! dlm1 = l-m + 1
            fac = fac / dsqrt(dlm*dlm1)
            plm(kp) = fac * plm(kp)
         enddo
      enddo

      end
      SUBROUTINE ASSLEG(P,N,X)

C     Subroutine computes associated Legendre polynomials as defined
C     by A.R. Edmonds (Angular Momentum in Quantum Mechanics):
C              P(l,m;x) =  (1-x^2)^(m/2) * (d/dx)^m { p(l,x) },
C     x is the usual coordinate (-1 < x < +1 ), n is the maximum
C     l-quantum number. No m-dependent phase.

C     The associated legendre function P(l,m;x) may be  accessed via
C     P( l*(l+1)/2 + m + 1 ).

C     Author P.E.S. WORMER (1982/2003)

      implicit real*8 (a-h,o-z)
      dimension p( (n+1)*(n+2)/2 )

      p(1) = 1.d0
      if (n .eq. 0) return
      sint = dsqrt(1.d0 - x*x)
      p(2) = x
      p(3) = sint
      if (n .eq. 1) return

      lm1 = 1
      lm  = 3
      dl  = 1.d0

      do  l=2,n
         dl = dl + 1.d0
         lm1 = lm1 + 1
C        (   lm1 = l*(l-1)/2 + 1 )
         lm  = lm + 1
C        (   lm = l*(l+1)/2 + 1 )

         p(lm) = x*p(lm1) - sint*p(lm1+1)/dl
C        ( p(l,0) = x*p(l-1,0) - dsqrt(1-x*x)*p(l-1,1)/l   )
         mmax = l-1
         dlm = dl
            do m=1,mmax
               lm1 = lm1 + 1
               lm  = lm + 1

               p(lm) = dlm*sint*p(lm1-1) + x*p(lm1)
C              ( p(l,m) = (l+m-1)*dsqrt(1-x*x)*p(l-1,m-1) + x*p(l-1,m))
               dlm = dlm + 1.d0
            enddo
         lm = lm + 1
         p(lm) = (dl+dl-1.d0)*sint*p(lm1)
C        (  p(l,l) = (2*l-1)*dsqrt(1-x*x)*p(l-1,l-1)    )

      enddo

      end
      FUNCTION THREEJP(J1,J2,J3,M1,M2,M3)

C     Wigner 3-j (only for integer quantum numbers).
C
C     Author: J. Tennyson, modified by P.E.S. Wormer
C
C     Important note:
C     ==============
C     Function is not protected against addressing outside BINOM!
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER( ZERO=0.D0, ONE=1.0D0)
      COMMON /BIN/  BINOM(0:200, 0:200), JMAX

C Fill out BINOM by Pascal's triangle relation (if necessary)
      IF (JMAX .EQ. 0) CALL PASCAL

      THREEJP = ZERO

      IF ((M1+M2+M3).NE.0)  RETURN
      I1 = -J1+J2+J3
      IF (I1.LT.0)  RETURN
      I2 =  J1-J2+J3
      IF (I2.LT.0)  RETURN
      I3 =  J1+J2-J3
      IF (I3.LT.0)  RETURN
      K1 =  J1+M1
      IF (K1.LT.0)  RETURN
      K2 =  J2+M2
      IF (K2.LT.0)  RETURN
      K3 =  J3+M3
      IF (K3.LT.0)  RETURN
      L1 =  J1-M1
      IF (L1.LT.0)  RETURN
      L2 =  J2-M2
      IF (L2.LT.0)  RETURN
      L3 =  J3-M3
      IF (L3.LT.0)  RETURN
      N1 = -J1-M2+J3
      N2 =  M1-J2+J3
      N3 =  J1-J2+M3
      IMIN = MAX0(-N1,-N2,0)
      IMAX = MIN0(L1,K2,I3)
      IF (IMIN .GT. IMAX)  RETURN
      SIGN = - ONE
      DO 30 I=IMIN,IMAX
         SIGN = -SIGN
         THREEJP = THREEJP + SIGN*BINOM(I1,N1+I)*
     1                          BINOM(I2,N2+I)*BINOM(I3,I)
   30 CONTINUE

      THREEJP = THREEJP * DSQRT(BINOM(J2+J2,I3)*BINOM(J1+J1,I2)
     1       / (BINOM(J1+J2+J3+1,I3)*DBLE(J3+J3+1)
     2       * BINOM(J1+J1,L1)*BINOM(J2+J2,L2)*BINOM(J3+J3,L3)))
      IF (MOD(N3+IMIN,2) .NE. 0) THREEJP = - THREEJP

      END

      SUBROUTINE PASCAL

C     Fill COMMON BIN with binomial coefficients (Pascal's triangle)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (KMAX = 200)
      COMMON /BIN/  BINOM(0:200, 0:200), JMAX

      IF ( KMAX .EQ. JMAX ) RETURN

C     Set jmax as sign that BINOM has been set:
      JMAX = KMAX

      BINOM(0,0) = 1.D0
      DO 20 I=1,JMAX
         BINOM(I,0) = 1.D0
         BINOM(I,I) = 1.D0
         DO 10 J=1,I-1
            BINOM(I,J) = BINOM(I-1,J-1) + BINOM(I-1,J)
   10    CONTINUE
   20 CONTINUE

      END
      BLOCK DATA BLPAS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /BIN/  BINOM(0:200, 0:200), JMAX
      DATA JMAX/0/
      END
      REAL*8 FUNCTION TTdamp(R, beta, n)

C     Tang-Toennies damping function.

      implicit real*8 (a-h,o-z)
      parameter (nmax = 20)
      real*8 factorial(0:nmax)
      logical first/.true./

      if (n .gt. nmax) then
          write(*, '(//a,i5/a,i1)')
     .   ' *** Error: calling TTdamp with 3rd parameter ', n,
     .   ' maximum for this parameter is ', nmax
         stop 16
      endif

      if (beta .ge. 100.d0) then ! no damping
         TTdamp = 1.d0
         return
      endif

!      if (first) then
         factorial(0) = 1.d0
         do i = 1,nmax
           factorial(i) = factorial(i-1)*dble(i)
         enddo
!         first = .false.
!      endif

      sum = 0.d0
      do k = 0, n
         sum = sum + ((beta*R)**k)/factorial(k)
      enddo

      TTdamp = 1.d0 - exp(-beta*R)*sum

      end
