      module ini_PESH2OH2
      integer mxlmb, i, ityp
      parameter(mxlmb=10000)
      integer icntrl, nlam, lam(4,mxlmb), p1, q1, p2, p
      integer mxlam
      integer kount	  
      character dir*80
      real*8 COEF(mxlmb),
     $     thx, phx, thpx, phpx,r, y
      logical corr_9d_flag, r12_flag
      logical :: ini_PESH2OH2_defined = .FALSE.
      logical :: first = .FALSE.
      REAL *8, PARAMETER :: pi = dacos(-1d0)	  
      end	  
	  
	  
!!!OPTION 2
      SUBROUTINE USER_DEFINED_TERMS(T,I,R)!! SUBROTUNE COMPUTES RADIAL COEFFICENTS OF THE PES EXPANSION  AT GIVEN DISTANCE
!!! INPUT : R - distance between COMs of particles, I - TERM NUMBER
!!! OUTPUT: T - value of coefficent 	  
      IMPLICIT NONE	  
      REAL*8 T,R
      INTEGER I

!!!! HERE WHAT USER PREPARED 	    
	  END SUBROUTINE USER_DEFINED_TERMS

	  
!!! OPTION 3: SIMILAR TO OPTION 2, BUT USER NEEDS TO PROVIDE POTENTIAL EXPANSION COEFFICEINS IN FILE "EXPAN_PES_TERMS.DAT"

!!! OPTION 4 
      SUBROUTINE USER_DEFINED_COEFFS(T,DTDR,I,R) !! SUBROTUNE COMPUTES RADIAL COEFFICENTS OF THE PES EXPANSION AND THEIR DERIVATIVES AT GIVEN DISTANCE
      IMPLICIT NONE
!!! INPUT : R - distance between COMs of particles, I - TERM NUMBER
!!! OUTPUT: T - value of coefficent, DTDR - its radial derivative 	  
      REAL*8 T,R,DTDR 
      INTEGER I
	  
!!! USER CODE	
      END SUBROUTINE USER_DEFINED_COEFFS 
	  
!! OPTION 1	  
      subroutine  USER_DEFINED_PES(V,R_d,rvib,rvib2,alpha,
     & beta,gamma,aalpha,bbeta
     & ,ggamma)
C -- NORMALISATION MOLSCAT POUR LES COEFFICIENTS  ----
C---------------------------------------------------------
      use ini_PESH2OH2
      implicit none
      REAL*8 V,R_d,alpha,beta,gamma,aalpha,bbeta,ggamma,t,rvib,rvib2
      REAL*8 t_body_fixed,tnormed
      REAL*8 r_inp(3),r_out(3),thppx,phppx	  
	  real*8,parameter :: bk_pi=4.0d0*atan(1.0d0)
!      EXTERNAL tnormed,t_body_fixed	
	  
c-------------------------------- INITIALIZATION
      IF(.not.ini_PESH2OH2_defined) THEN
      icntrl=-1
      dir='./H2Ov000_H2v0/'
!      corr_9d_flag=.true.
!      r12_flag=.true.
c      write(6,'(a,$)')
c     $     'Enter corr_9d_flag (t/f) and r12_flag (t/f) ? '
c      read(5,*) corr_9d_flag, r12_flag
      r=0
      do i=1,mxlmb
         coef(i)=0
      enddo
      mxlam = mxlmb
      call potenl(icntrl, 
     $     mxlam, nlam, lam, r, coef,ityp)
      write(6,*) 'NLAM', nlam,'  MXLAM', MXLAM

c-------------------------------- PES EVALUATION
      ini_PESH2OH2_defined = .TRUE.
      ENDIF
      icntrl=0
	  
!         thx = dacos(dcos(beta))*180d0/pi
		 
!         phx = datan2(dsin(beta)*dsin(gamma),
!     & dsin(beta)*dcos(gamma))*180d0/pi
	 
!         thpx = dacos(dcos(aalpha)*dcos(bbeta)*dsin(beta)
!     & -dsin(bbeta)*dcos(beta))*180d0/pi
	 
!         phpx = datan2(dsin(aalpha)*dcos(bbeta)*dcos(gamma)
!     & -dcos(aalpha)*dcos(bbeta)*dcos(beta)*dsin(gamma)
!     & -dsin(bbeta)*dsin(beta)*dsin(gamma),
!     & -dcos(aalpha)*dcos(bbeta)*dcos(beta)*dcos(gamma)
!     & -dsin(bbeta)*dsin(beta)*dcos(gamma)
!     & -dsin(aalpha)*dcos(bbeta)*dsin(gamma))*180d0/pi

!         thpx = aalpha*180d0/pi
	 
!         phpx = bbeta*180d0/pi
	 
!         IF(phx.lt.0d0) phx = 360.0d0+phx		 
!         IF(phpx.lt.0d0) phpx = 360.0d0+phpx	

		call xyz(alpha,beta,gamma,aalpha,bbeta,ggamma,thx,phx,thpx,phpx)  
		thx=thx*180.0d0/pi
		phx=phx*180.0d0/pi
		thpx=thpx*180.0d0/pi
		phpx=phpx*180.0d0/pi
		
         call potenl(icntrl, 
     $        mxlam, nlam, lam, R_d, coef,ityp)
         y=0d0
         kount = 0		 
         do i=1,nlam
            p1=lam(1,i)
            q1=lam(2,i)
            p2=lam(3,i)
            p =lam(4,i)
            y=y+t(p1,q1,p2,p,thx,phx,thpx,phpx)*coef(i)
         enddo
      V = y
 
      end

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
C     cFUNCTION PM1(LIN,MIN,COSTH)
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

      subroutine POTENL(ICNTRL, 
     $     NDIM, NLAM, LAM, RD, P, ITYP)
c
c -------------------------------------------------------------------
c Vibrationally averaged 5D PES expansion for H2O-H2,
c Valiron et al, Version 1.1 january 2004.
c -------------------------------------------------------------------
c
c  ARGUMENT LIST
c
c  ICNTRL      (input)   -1  mandatory initialization
c                         0  to evaluate P(R)
c  DIR         (input)   slash/ terminated directory containing data   (-1)
c  CORR_9D_FLAG (input)   Add expectation value of 9-D PES modulation  (-1)
c  R12_FLAG     (input)   Add R12 calibration data                     (-1)
c  NDIM        (input)   dimension of LAM(4,*) and P(*)
c  NLAM        (output)  number of terms of the expansion 
c  LAM(4,NDIM) (output)  (4,NDIM) p1, q1, p2, p expansion terms        (-1)
c  RD          (input)   intermolecular distance (a.u.)                (0)
c  P(NDIM)     (output)  coefficients for angular expansion in cm-1    (0)
c              properly interpolated or extrapolated at distance R 
c              assuming normalized expansion terms (beware the normalization
c              error in PMG94's paper).
c              These coefficients come in the ordering indicated in LAM. 
c              This ordering reflects the list of terms given in the
c              input file mesh.dat and can thus be easily changed. 
c              Terms can easily be removed or added in mesh.dat provided
c              they form a subset of the original 149 term expansion in
c              files ref.DDDD.dat. 
c
c     (-1)  relevant for initialization only (ICNTRL=-1)
c     (0)   relevant for PES evaluation only (ICNTRL=0)
c
c  RELATED FILES
c
c     mesh.dat, ref.DDDD.dat, cal.DDDD.dat and vib.DDDD.dat
c     in directory DIR.
c
      implicit none
c      integer icntrl, ndim, nlam, lam(4,ndim), np1, np1_r12

c MLD change
      integer idim
      parameter (idim=150)
      integer icntrl, ndim, nlam, lam(4*idim), np1, np1_r12
      integer ityp
c MLD end change

      character*80 dir
      real*8 rd, p(ndim)
      integer ndx, ntx, nttx, ndx_r12, ntx_r12
      parameter (ndx=30, ntx=150, nttx=150, ndx_r12=20, ntx_r12=9)
      integer ndist, id, Ngeom, Nfunc, Ngeom0, Nfunc0, i, it, nt, k, k2,
     $     l0(4,ntx), l(4,ntx), nex,
     $     ll_r12(4,ntx_r12), ndist_r12, nt_r12, Ngeom_r12(ndx_r12)
      integer ndist_c9d, Ngeom_c9d, Nfunc_c9d, l_c9d(8,ntx), ind, ind2
      real*8 fdist, Dist(ndx), Dist_c9d(ndx), Dist_r12(ndx_r12),
     $     Sinv(ndx), RMS(ndx), ERR(ndx)
      real*8 fdist_c9d
      real*8 coef(ndx,ntx),
     $     cc(ndx,nttx), c1(ndx,nttx), c2(ndx,nttx), c3(ndx,nttx),
     $     cc_r12(ndx_r12,ntx_r12), c1_r12(ndx_r12,ntx_r12),
     $     c2_r12(ndx_r12,ntx_r12), c3_r12(ndx_r12,ntx_r12),
     $     cc_r12_old(ndx_r12,ntx_r12), cc_r12_shift(ntx_r12)
      real*8 cc_c9d(ndx,ntx)
      integer ntt, itt, in, ll(4,nttx), indll(nttx), indll_r12(ntx_r12)
      logical chk(ntx), r12_flag, corr_9d_flag, error
      integer n0000, kleft, kright,
     $     kleft_r12, kright_r12
      real*8 a(nttx), alpha(nttx), dleft, dright, hinvl, hinvr
      real*8 c(nttx), b(nttx),
     $     v0000, ccinv0000, fact, 
     $     c_r12(ntx_r12), b_r12(ntx_r12),
     $     dleft_r12, dright_r12, hinvl_r12, hinvr_r12, p_r12(ntx_r12),
     $     pm_r12(ntx_r12)
      character cfit*20, comment*132, name*132
      character cfit_c9d*20
      integer le
      real*8 h, r, rinv, aniso, tres, thres
      real*8 f, u, fu, Pi, Pis2
c MLD change
      real*8 epsil, rm
      integer q1, p1
      real*8 one, two, dq1, xn , z
      parameter (one=1.0d0, two=2.0d0)
      data rm /0.52917706D0/, epsil /1.0d0/
      namelist /POTL/dir,corr_9d_flag,r12_flag
c MLD endchange      

c     funct(u)=(1-cos(u*Pi))/2 ; f(u)=funct(funct(u))
      f(u)=(1-cos((1-cos(u*Pi))*Pis2))*0.5d0
      save

      if (ICNTRL.eq.0) goto 1000   
c-------------------------------------------------------------------
c Initialization (ICNTRL=-1)
c-------------------------------------------------------------------
      if (ICNTRL.ne.-1) stop 'invalid ICNTRL, expected -1'
	  open(unit=6,file='temp.out')
      dir = './PES_H2OH2/H2Ov000_H2v0/'
      corr_9d_flag = .true.
      r12_flag = .true.
!      read(*,POTL)
      Pi=acos(-1.0d0)
      Pis2=Pi/2

c Set directory containing the fitted and mesh data (slash terminated)

      le=len_trim(dir)

c Open mesh.dat

      name=dir(:le)//'mesh.dat'
      open(unit=1, file=name, status='old')

      read(1,*) comment
      write(6,*) '___________________________________________________'
      write(6,*)
      write(6,*)
     $     'H20-H2 5D PES, Version 1.1, Valiron et al, january 2004.'
      write(6,*)
      write(6,*) comment(1:len_trim(comment))
      write(6,*) '___________________________________________________'

c-------------------------------------------------------------------
c  Cheap CCSD(T)/CP rigid-body reference
c-------------------------------------------------------------------

c Read fitted data in ref.DDDD.dat files

      read(1,*) ndist
      write(6,*)
      write(6,*)
     $     '********************** CCSD(T)/CP rigid-body reference'
      write(6,*)
      write(6,*) 'Number of distances', ndist
      if (ndist.gt.ndx) stop 'increase ndx'

      do id=1,ndist

         read(1,*) fdist, cfit
         write(6,'(f8.3,1x,2a)') fdist, dir(:le), cfit

         name=dir(:le)//cfit
         open(unit=2, file=name, status='old')

c Skip two comment lines
         read(2,*)
         read(2,*)

c Read data at dist(id)
         read(2,*) Dist(id), Ngeom, Nfunc, Sinv(id), RMS(id), ERR(id)
         if (dist(id).ne.fdist) then
            write(6,*) 'Inconsistent distance in ref.DDDD.dat'
            write(6,*) fdist, dist(id)
            stop
         endif
         if (id.eq.1) then
            Ngeom0=Ngeom
            Nfunc0=Nfunc
            nt=Nfunc
            write(6,*) 'Ngeom = ', Ngeom, '    Nfunc = ', Nfunc
            if (Nfunc.gt.ntx) stop 'increase ntx'
         else
            if (Ngeom0.ne.Ngeom) stop 'inconsistent Ngeom'
            if (Nfunc0.ne.Nfunc) stop 'inconsistent Nfunc'
            if (Dist(id).le.Dist(id-1)) stop 'non monotonic Dist'
         endif
         
         do it=1,nt
            read(2,*) (l(i,it),i=1,4), coef(id,it)
            p1=l(1,it)
            q1=l(2,it)
            if (q1.eq.0) then
             dq1=1
            else
             dq1=0
            endif
            xn = two / (one+dq1) / (2*p1+1)
            z = sqrt(xn)
            coef(id,it) = coef(id,it)/z
         enddo

         if (id.eq.1) then
            do it=1,nt
               do i=1,4
                  l0(i,it)=l(i,it)
               enddo
            enddo
         else
            k=0
            do it=1,nt
               do i=1,4
                  k=k+abs(l0(i,it)-l(i,it))
               enddo
            enddo
            if (k.ne.0) stop 'inconsistent expansion terms'
         endif

         close(unit=2)

      enddo

      write(6,*)
      write(6,*) ndist, ' Distances'
      write(6,*)
      write(6,*)
     $     '   Dist  Ngeom Nfunc  Sinv       RMS         ERR'
      do id=1,ndist
         write(6,'(f8.3,2i6,f8.3,2f12.6)')
     $        Dist(id), Ngeom, nt, Sinv(id), RMS(id), ERR(id)
      enddo

      write(6,*)
      write(6,*) 'Original terms from ref.DDDD.dat files', nt
      write(6,*)
      write(6,'(3(i4,2x,4i3,3x))') (it, (l(i,it), i=1,4), it=1,nt)
!	  print*,'check1'
c-------------------------------------------------------------------
c  Vibrationally averaged correction
c
c  Present code assumes the radial mesh is the same as for the 
c  reference PES. However the angular expansion could be a subset
c  of the angular expansion for the reference PES. 
c-------------------------------------------------------------------

c read <delta9d> coefficients to add the vibrationally averaged
c correction which is an average of the 9d correction PES over a
c given vibrational wavefunction for H2O and H2

      write(6,*)
      write(6,*) '********************** vibrational corrections'
      write(6,*)

c     corr_9d_flag was moved to argument list

      read(1,*) ndist_c9d
c      read(1,*), ndist_c9d
      if (ndist_c9d.ne.ndist) then
         write(6,*) 'Present version requires the same radial grid',
     $        ' for reference and vibrational corrections'
         write(6,*) 'ndist, ndist_c9d', ndist, ndist_c9d
         stop
      endif

      do id=1,ndist_c9d
c         read(1,*), fdist_c9d, cfit_c9d
         read(1,*) fdist_c9d, cfit_c9d
         name = dir(:le)//cfit_c9d
         write(6,'(f8.3,1x,2a)') fdist_c9d, dir(:le), cfit_c9d
         
         open(unit = 2, file=name, status='old')
c Skip two comment lines
         read(2,*)
         read(2,*)

c Read data at dist(id)
         read(2,*) Dist_c9d(id), Ngeom_c9d, Nfunc_c9d          
         if (dist_c9d(id).ne.fdist_c9d) then
            write(6,*) 'Inconsistent distance in vib.DDDD.dat'
            write(6,*) fdist, dist_c9d(id)
            stop
         endif
         if (dist_c9d(id).ne.dist(id)) then
            write(6,*) 'Present version requires the same radial grid',
     $           ' for reference and vibrational corrections'
            write(6,*) 'id, dist_c9d(id), dist(id)',
     $           id, dist_c9d(id), dist(id)
            stop
         endif
         do it = 1,Nfunc_c9d
            read(2,*) (l_c9d(i,it), i=1,4), cc_c9d(id,it)
            p1=l_c9d(1,it)
            q1=l_c9d(2,it)
            if (q1.eq.0) then
             dq1=1
            else
             dq1=0
            endif
            xn = two / (one+dq1) / (2*p1+1)
            z = sqrt(xn)
            cc_c9d(id,it) = cc_c9d(id,it)/z
         enddo
      enddo
      
c Add averaged C9d correction to pure 5d coeffs
c N^2 loop. Could be clever but negligible cpu time anyway...

      if (corr_9d_flag) then
         write(6,*)
         write(6,*) '-------------------------------'
         write(6,*) 'Vibrational correction applied'
         write(6,*) '-------------------------------'
         write(6,*)
         do id=1,ndist_c9d
            do ind2 = 1, Nfunc_c9d
               error=.true.
               do ind = 1, nt
                  if (l(1,ind) .eq. l_c9d(1,ind2) .and.
     $                 l(2,ind) .eq. l_c9d(2,ind2) .and.
     $                 l(3,ind) .eq. l_c9d(3,ind2) .and.
     $                 l(4,ind) .eq. l_c9d(4,ind2)) then
                     coef(id,ind) = coef(id,ind) + cc_c9d(id,ind2)
                     error=.false.
                  endif
               enddo
               if (error) then
                  write(6,*) 'Vibrational term ind2=', ind2,
     $                 ' not found in reference expansion'
                  write(6,*) (l_c9d(i,ind2), i=1,4)
                  stop
               endif
            enddo
         enddo
         
      else 
         write(6,*)
         write(6,*) '-------------------------------'
         write(6,*) 'Vibrational correction ignored'
         write(6,*) '-------------------------------'
         write(6,*)
      endif                    
c 
c end of vibrationally averaged correction
c____________________________________________________________________
!	  print*,'check2'
c Select final expansion terms and corresponding pointers to original data

      read(1,*) ntt
      write(6,*)
      write(6,*) 'Final expansion terms', ntt
      write(6,*)
      if (ntt.gt.nttx) stop 'increase nttx'
      read(1,*) ((ll(i,itt), i=1,4), itt=1,ntt)
      write(6,'(3(i4,2x,4i3,3x))') (itt, (ll(i,itt), i=1,4), itt=1,ntt)

      do it=1,ntx
         chk(it)=.false.
      enddo
      do itt=1,nttx
         indll(itt)=0
      enddo
      do itt=1,ntt
         do it=1,nt
            if (       l(1,it).eq.ll(1,itt)
     $           .and. l(2,it).eq.ll(2,itt)
     $           .and. l(3,it).eq.ll(3,itt)
     $           .and. l(4,it).eq.ll(4,itt) ) then
               if (chk(it)) then
                  write(6,*) 'l =', (l(i,it),i=1,4)
                  write(6,*) 'll=', (ll(i,itt),i=1,4)
                  write(6,*) 'itt=', itt, '     it=', it
                  stop 'error -- chk(it) already set'
               else
                  indll(itt)=it
                  chk(it)=.true.
               endif
            endif
         enddo
      enddo

      write(6,*)
      write(6,*) 'Pointers to original fit expansion terms'
      write(6,*)
      write(6,'(5(i4,a,i3,3x))') (itt, ' ->', indll(itt), itt=1,ntt)
      
      error=.false.
      do itt=1,ntt
         if (indll(itt).eq.0) then
            error=.true.
            write(6,*)
     $           'invalid pointer for term ll=', (ll(i,itt),i=1,4)
         endif
      enddo
      if (error) stop

      write(6,*)
      write(6,*) 'Copy original expansion coeffs to final cc matrix'
      do itt=1,ntt
         it=indll(itt)
         do id=1,ndist
            cc(id,itt)=coef(id,it)
         enddo
      enddo

      n0000=0
      do itt=1,ntt
         if (       ll(1,itt).eq.0
     $        .and. ll(2,itt).eq.0
     $        .and. ll(3,itt).eq.0
     $        .and. ll(4,itt).eq.0 ) then
            n0000=itt
         endif
      enddo
      if (n0000.eq.0) stop 'n0000 term missing in expansion'

c Retrieve and set up long range terms as C/R**beta

      write(6,*)
      write(6,*) 'Set up long range extrapolation'
      write(6,*) 'All terms above 0.001 at R=14 and 15 au are',
     $     ' extrapolated'
      write(6,*)
      write(6,*) ' #term          term            C                beta'

c All terms above thres=1d-2 at R = 14 and 15 au are extrapolated 
c by a power law obtained from values at 14 and 15 au.

      do itt=1,ntt
         b(itt)=0
         c(itt)=0
      enddo
      thres=1d-2
      nex=0
      do itt=1, ntt
         if (abs(cc(ndist-1,itt)).ge.thres
     $        .and. abs(cc(ndist,itt)).ge.thres) then
            nex=nex+1
            b(itt)=log(cc(ndist-1,itt)/cc(ndist,itt))
     $           /log(dist(ndist)/dist(ndist-1))
            c(itt)=cc(ndist,itt)*dist(ndist)**b(itt)
            write(6,'(i6,4i6,1pe16.6,0pf12.6)')
     $           itt, (ll(i,itt),i=1,4), c(itt), b(itt)
         endif
      enddo

      write(6,*)
      write(6,*) nex, ' terms extrapolated beyond R=15 au'

!		print*,'check3'
c Initialize extrapolation at short range

      write(6,*)
      write(6,*) 'Monitor exponential short range extrapolation'
      tres=20000
      write(6,'(1x,a,1pe12.2,a,0pf8.3)') 'Select terms above TRES =',
     $     tres, '  at R=', dist(1)
      do itt=1,ntt
         alpha(itt)=log(cc(1,itt)/cc(2,itt))/(dist(2)-dist(1))
         a(itt)=cc(1,itt)*exp(alpha(itt)*dist(1))
         if (cc(1,itt).gt.tres) then
            write(6,'(i4,4x,4i3,1pe14.3,0pf12.3)')
     $           itt, (ll(i,itt),i=1,4), a(itt), alpha(itt)
         endif
      enddo
      write(6,*)
      write(6,*)'    Dist      Isotropic   Sigma(Aniso)    Ratio'
      r=0
      do 
         aniso=0
         do itt=1,ntt
            if (itt.ne.n0000 .and. cc(1,itt).gt.tres) then
               aniso=aniso+abs(a(itt))*exp(-alpha(itt)*r)
            endif
         enddo
         v0000=a(n0000)*exp(-alpha(n0000)*r)
         write(6,'(f10.3, 1p2e14.3, 0pf10.2)')
     $        r, v0000, aniso, aniso/v0000
         r=r+0.1d0
         if (r.gt.dist(3)+0.0001) exit
      enddo
      write(6,*)
     $     '  ==> Anisotropic terms grow faster than isotropic one;'
      write(6,*)
     $     '  ==> blind extrapolation of anisotropic terms is DANGEROUS'
      write(6,*)
     $     '  ==> proportional extrapolation more secure'
!	  print*,'check4'
c set up smooth transition domains using step function f
      read(1,*) fact
      dleft=dist(1)
      kleft=2
      hinvl=1/(dist(kleft+1)-dist(1))
      dright=dist(ndist)
      kright=ndist-1
      hinvr=1/(dist(ndist)-dist(kright))
      write(6,*)
      write(6,*) 'Set up smooth transition domains...'
      write(6,'(1x,a,f10.3 )')
     $     'Anisotropy fraction kept at shorter range', fact
      write(6,'(1x,a,f10.3,a,f10.3)')
     $     'Left  domain', dist(1),  '  -> ', dist(kleft+1)
      write(6,'(1x,a,f10.3,a,f10.3)')
     $     'Right domain', dist(kright), '  -> ', dist(ndist)


c In order to avoid the ondulation of the spline in the last
c interval (where the interaction energy is small), a single
c additional point obtained from the long-range extrapolation at
c r=dist(ndist)+1 (i.e. 16 au) is added to cc.
c
      np1=ndist+1
      if (np1.gt.ndx) stop 'increase ndx'
      dist(np1)=dist(ndist)+1
      do itt=1,ntt
         cc(np1,itt)=c(itt)*(1/dist(np1))**b(itt)
      enddo

c Initialize spline coefficients for original data

      write(6,*)
      write(6,*) 'Set up spline coefficients...'
      do itt=1,ntt
         call cubspl (np1,dist,
     $        cc(1,itt),c1(1,itt),c2(1,itt),c3(1,itt),0,0)
c scale c2 and c3 to obtain directly taylor coefficients
c optionally one may put a unit conversion there
c be clever enough in this case to take into account also the unit 
c conversions in the short and long range extraoplations...
c be also consistent with the R12 correction too...
         do id=1,np1
!            cc(id,itt) = cc(id,itt)
!            c1(id,itt) = c1(id,itt)
            c2(id,itt) = c2(id,itt) / 2.d0
            c3(id,itt) = c3(id,itt) / 6.d0
         enddo
      enddo

!	  print*,'check5'
c-------------------------------------------------------------------
c  R12-CCSD(T) corrections
c-------------------------------------------------------------------

c  Read fitted data in cal.DDDD.dat files

c     r12_flag has been moved to argument list

      read(1,*) ndist_r12
      write(6,*)
      write(6,*) '********************** R12-CCSD(T) corrections'
      write(6,*)
      write(6,*) 'Number of distances', ndist_r12
      if (ndist_r12.gt.ndx_r12) stop 'increase dx_r12'

      do id=1,ndist_r12

         read(1,*) fdist, cfit
         write(6,'(f8.3,1x,2a)') fdist, dir(:le), cfit

         name=dir(:le)//cfit
         open(unit=2, file=name, status='old')

         read(2,*)
         read(2,*)
         read(2,*) Dist_r12(id),
     $        Ngeom_r12(id), Nfunc, Sinv(id), RMS(id), ERR(id)
         if (dist_r12(id).ne.fdist)
     $        stop 'inconsistent distance in cal.DDDD.dat'
         if (id.eq.1) then
            Nfunc0=Nfunc
            nt_r12=Nfunc
            write(6,*) 'Nfunc = ', Nfunc
            if (Nfunc.gt.ntx_r12) stop 'increase ntx_r12'
         else
            if (Nfunc0.ne.Nfunc) stop 'inconsistent Nfunc'
            if (Dist_r12(id).le.Dist_r12(id-1))
     $           stop 'non monotonic Dist'
         endif

c     R12 calculations have been performed with 2 different H basis set:
c     1. For R=3, 3.5, 4.0 and 4.5 au, 'old' H basis set using the PMG
c     grid extended with 20 random geometries (giving a total of 103
c     geometries)
c     2. For R=5-12au, a 'new' H basis set using 40 random geometries.
c     Below the shift between old and new R12 calculations is computed
c     at 5au (id=5) and is used below to correct 'old' R12 values at
c     R=3-4.5 au.

         if (id.eq.5) then
            do it=1,nt_r12
               read(2,*) (ll_r12(i,it),i=1,4), cc_r12(id,it), 
     $              cc_r12_old(id,it)
               cc_r12_shift(it)=cc_r12(id,it)-cc_r12_old(id,it)

            p1=ll_r12(1,it)
            q1=ll_r12(2,it)
            if (q1.eq.0) then
             dq1=1
            else
             dq1=0
            endif
            xn = two / (one+dq1) / (2*p1+1)
            z = sqrt(xn)
            cc_r12_shift(it) = cc_r12_shift(it)/z
            cc_r12(id,it) = cc_r12(id,it)/z
            enddo
         else
            do it=1,nt_r12
               read(2,*) (ll_r12(i,it),i=1,4), cc_r12(id,it)
            p1=ll_r12(1,it)
            q1=ll_r12(2,it)
            if (q1.eq.0) then
             dq1=1
            else
             dq1=0
            endif
            xn = two / (one+dq1) / (2*p1+1)
            z = sqrt(xn)
            cc_r12(id,it) = cc_r12(id,it)/z
            enddo
         endif
c
         if (id.eq.1) then
            if (nt_r12.gt.ntx) stop 'increase ntx'
            do it=1,nt_r12
               do i=1,4
                  l0(i,it)=ll_r12(i,it)
               enddo
            enddo
         else
            k=0
            do it=1,nt_r12
               do i=1,4
                  k=k+abs(l0(i,it)-ll_r12(i,it))
               enddo
            enddo
            if (k.ne.0) stop 'inconsistent R12 expansion terms'
         endif

         close(unit=2)

      enddo

c     The R12 shift is added to distances R=3.0, 3.5, 4.0, 4.5 au.

      do id=1,4
         do it=1,nt_r12
            cc_r12(id,it)=cc_r12(id,it)+cc_r12_shift(it)
         enddo
      enddo
c
      write(6,*)
      write(6,*) ndist_r12, ' R12 Distances'
      write(6,*)
      write(6,*)
     $     '   Dist  Ngeom Nfunc  Sinv       RMS         ERR'
      do id=1,ndist_r12
         write(6,'(f8.3,2i6,f8.3,2f12.6)') Dist_r12(id),
     $        Ngeom_r12(id), nt_r12, Sinv(id), RMS(id), ERR(id)
      enddo

      write(6,*)
      write(6,*) 'Original terms from cal.DDD.dat files', nt_r12
      write(6,*)
      write(6,'(3(i4,2x,4i3,3x))')
     $     (it, (ll_r12(i,it), i=1,4), it=1,nt_r12)

c set up smooth transition domains using step function f
      dleft_r12=dist_r12(1)
      kleft_r12=1
      hinvl_r12=1/(dist_r12(kleft_r12+1)-dist_r12(1))
      dright_r12=dist_r12(ndist_r12)
      kright_r12=ndist_r12-1
      hinvr_r12=1/(dist_r12(ndist_r12)-dist_r12(kright_r12))
      write(6,*)
      write(6,*) 'Set up smooth transition domains for R12...'
      write(6,'(1x,a,f10.3,a,f10.3)') 'Left  domain',
     $     dist_r12(1),  '  -> ', dist_r12(kleft_r12+1)
      write(6,'(1x,a,f10.3,a,f10.3)') 'Right domain',
     $     dist_r12(kright_r12), '  -> ', dist_r12(ndist_r12)


c Retrieve and set up long range terms as C/R**beta
c All terms are extrapolated from values at 11 and 12 au

      write(6,*)
      write(6,*) 'Set up long range extrapolation for R12 correction'
      write(6,*) ' #term          term            C                beta'
      
      do itt=1,nt_r12
         b_r12(itt)=log(cc_r12(ndist_r12-1,itt)
     $        /cc_r12(ndist_r12,itt))
     $        /log(dist_r12(ndist_r12)/dist_r12(ndist_r12-1))
         c_r12(itt)=cc_r12(ndist_r12,itt)
     $        *dist_r12(ndist_r12)**b_r12(itt)
         write(6,'(i6,4i6,1pe16.6,0pf12.6)') itt, 
     $        (ll_r12(i,itt),i=1,4), c_r12(itt), b_r12(itt)
      enddo


c In order to avoid ondulation of the spline in the last interval
c (where the interaction energy is small), a single additional point
c obtained from the long-range extrapolation at r=dist(ndist)+1
c (ie 13 au) is added to cc_r12.
c
      np1_r12=ndist_r12+1
      dist_r12(np1_r12)=dist_r12(ndist_r12)+1
      do itt=1,nt_r12
         cc_r12(np1_r12,itt)=c_r12(itt)*
     $        (1/dist_r12(np1_r12))**b_r12(itt)
      enddo


c Initialize spline coefficients for R12 corrections

      write(6,*)
      write(6,*) 'Set up spline coefficients...'
      do itt=1,nt_r12
         call cubspl (np1_r12,dist_r12,cc_r12(1,itt),
     $        c1_r12(1,itt),c2_r12(1,itt),c3_r12(1,itt),0,0)
c scale c2 and c3 to obtain directly taylor coefficients
         do id=1,np1_r12
!            cc_r12(id,itt) = cc_r12(id,itt)
!            c1_r12(id,itt) = c1_r12(id,itt)
            c2_r12(id,itt) = c2_r12(id,itt) / 2.d0
            c3_r12(id,itt) = c3_r12(id,itt) / 6.d0
         enddo
      enddo


c Set up pointers to connect the R12 corrections to the reference expansion

      if (ntx_r12.gt.ntx) stop 'ntx should be larger than ntx_r12'
      do it=1,nttx
         chk(it)=.false.
      enddo
      do itt=1,ntx_r12
         indll(itt)=0
      enddo
      do itt=1,nt_r12
         do it=1,ntt
            if (       ll(1,it).eq.ll_r12(1,itt)
     $           .and. ll(2,it).eq.ll_r12(2,itt)
     $           .and. ll(3,it).eq.ll_r12(3,itt)
     $           .and. ll(4,it).eq.ll_r12(4,itt) ) then
               if (chk(it)) then
                  write(6,*) 'll    =', (ll(i,it),i=1,4)
                  write(6,*) 'll_r12=', (ll_r12(i,itt),i=1,4)
                  write(6,*) 'itt=', itt, '     it=', it
                  stop 'error -- chk(it) already set'
               else
                  indll_r12(itt)=it
                  chk(it)=.true.
               endif
            endif
         enddo
      enddo

      write(6,*)
      write(6,*) 'R12 correction pointers to final expansion'
      write(6,*)
      write(6,'(5(i4,a,i3,3x))')
     $     (itt, ' ->', indll_r12(itt), itt=1,nt_r12)
      do itt=1,nt_r12
         if (indll_r12(itt).eq.0) stop 'some invalid pointer(s)'
      enddo

      if (r12_flag) then
         write(6,*) 
         write(6,*) '----------------------'
         write(6,*) 'R12 correction applied'
         write(6,*) '----------------------'
         write(6,*)
      else
         write(6,*)
         write(6,*) '----------------------'
         write(6,*) 'R12 correction ignored'
         write(6,*) '----------------------'
         write(6,*)
      endif

c-------------------------------------------------------------------

c return proper data
c MLD change
      if (ntt.gt.ndim) stop 'increase NDIM'
      nlam=ntt
      ind = 1
      do itt=1,ntt
         do i=1,4
            lam(ind)=ll(i,itt)
            ind = ind + 1
         enddo
      enddo
      close(unit=1)
      RD = RM
      P(1) = EPSIL
      write(6,*) 'p(1)'
      NDIM = NTT

c MLD end change

      write(6,*)
      write(6,*) 'Initialization done.'
      write(6,*)
      return      

!	  print*,'check6'
c-------------------------------------------------------------------
c Interpolation or extrapolation of the PES expansion (ICNTRL=0)
c-------------------------------------------------------------------

 1000 if (ntt.gt.ndim) stop 'invalid NDIM'
      nlam=ntt
      r=rd

c------------------------------------------ Reference PES

c short range extrapolation
c fact=0   -->   squeeze completely anisotropy for RD <= dleft
c fact=1   -->   keep anisotropy for RD <= dleft in a proportional way:
c the extrapolation of the isotropic term is used for all terms with a 
c scaling factor computed at r=dist(2) (ie 3.25 au)

      if (r.lt.dleft) then
         v0000=a(n0000)*exp(-alpha(n0000)*r)
         ccinv0000=1/cc(2,n0000)
         do itt=1,ntt
            if (itt.eq.n0000) then
               p(itt)=v0000
            else
               p(itt)=fact*v0000*cc(2,itt)*ccinv0000
            endif
         enddo

c long range extrapolation
      else if (r.gt.dright) then
         Rinv=1/r
         do itt=1,ntt
            if (c(itt).ne.0) then
               p(itt)=c(itt)*Rinv**b(itt)
            else
               p(itt)=0
            endif
         enddo

c spline domain
      else
         call splget (ndist,dist,r,k)
         h = r - dist(k)
         do itt=1,ntt
c remember y2 and y3 have been scaled for optimisation
            p(itt) = cc(k,itt)+h*
     $           (c1(k,itt)+h*(c2(k,itt)+h*c3(k,itt)))
         enddo

c branch smooth step function for first interval
         if (k.le.kleft) then
            fu=f((r-dist(1))*hinvl)
            v0000=a(n0000)*exp(-alpha(n0000)*r)         
            ccinv0000=1/cc(2,n0000)
            do itt=1,ntt
               if (itt.eq.n0000) then
                  p(itt)=fu*p(itt)+(1-fu)*v0000
               else
                  p(itt)=fu*p(itt)
     $                 +(1-fu)*fact*v0000*cc(2,itt)*ccinv0000
               endif
            enddo
            
c branch smooth step function for last interval
         else if (k.ge.kright) then
            fu=f((r-dist(kright))*hinvr)
            do itt=1,ntt
               p(itt)=(1-fu)*p(itt)
            enddo

            Rinv=1/r
            do itt=1,ntt
               if (c(itt).ne.0) then
                  p(itt)=p(itt)+fu*c(itt)*Rinv**b(itt)
               endif
            enddo

         endif
      endif
c------------------------------------------ R12 corrections

      if (.not.r12_flag) then
         goto 2000
      endif


c To be consistent with the reference PES short-range extrapolation, an 
c additional R12 point pm_r12(itt) at r=dist(2) is obtained from the spline

      call splget (ndist_r12,dist_r12,dist(2),k2)
      h = dist(2) - dist_r12(k2)
      do itt=1,nt_r12
c     remember y2 and y3 have been scaled for optimisation
         pm_r12(itt) = cc_r12(k2,itt)+h*
     $        (c1_r12(k2,itt)+h*(c2_r12(k2,itt)+h*c3_r12(k2,itt)))
      enddo

c short range extrapolation
c fact=0   -->   squeeze completely anisotropy for RD <= dleft
c fact=1   -->   keep anisotropy for RD <= dleft in a proportional way
c the extrapolation of the isotropic term (reference PES) is used for all 
c terms with a scaling factor computed at r=dist(2) (ie 3.25 au)

      if (r.lt.dleft_r12) then
         do itt=1,nt_r12
            it=indll_r12(itt)
            p(it)=p(it)+ fact*v0000*pm_r12(itt)*ccinv0000
         enddo

c long range extrapolation
      else if (r.gt.dright_r12) then
         Rinv=1/r
         do itt=1,nt_r12
            it=indll_r12(itt)
            p(it)=p(it)+c_r12(itt)*Rinv**b_r12(itt)
         enddo

c spline domain
      else
         call splget (ndist_r12,dist_r12,r,k2)
         h = r - dist_r12(k2)
         do itt=1,nt_r12
c remember y2 and y3 have been scaled for optimisation
            p_r12(itt) = cc_r12(k2,itt)+h*
     $           (c1_r12(k2,itt)+h*(c2_r12(k2,itt)+h*c3_r12(k2,itt)))
         enddo

c branch smooth step function for first interval
         if (k2.le.kleft_r12) then
            fu=f((r-dist_r12(1))*hinvl_r12)
            do itt=1,nt_r12
               p_r12(itt)=fu*p_r12(itt)+(1-fu)*fact*v0000*pm_r12(itt)
     $              *ccinv0000
            enddo
            
c branch smooth step function for last interval
         else if (k2.ge.kright_r12) then
            fu=f((r-dist_r12(kright_r12))*hinvr_r12)
            Rinv=1/r
            do itt=1,nt_r12
               p_r12(itt)=(1-fu)*p_r12(itt)
     $              +fu*c_r12(itt)*Rinv**b_r12(itt)    
            enddo

         endif

c apply R12 corrections
         do itt=1,nt_r12
            it=indll_r12(itt)
            p(it)=p(it)+p_r12(itt)
         enddo

      endif
	

 2000 end





      subroutine cubspl (n,tau,c1,c2,c3,c4,ibcbeg,ibcend)
c
c******  piecewise cubic spline interpolants computation; adapted from
c  'a practical guide to splines' , carl de boor , applied mathematical
c  sciences, springer-verlag, vol.27, p57-59 (1978).
c
c     ************************* input **************************
c
c     n = number of data points. assumed to be .ge. 2.
c     (tau(i), c1(i), i=1,...,n) = abscissae and ordinates of the
c        data points. tau is assumed to be strictly monotonous.
c     ibcbeg, ibcend = boundary condition indicators, and
c     c2(1) , c2(n)  = boundary condition information. specifically,
c        ibcbeg = 0  means no boundary condition at tau(1) is given.
c           in this case, the not-a-knot condition is used, i.e. the
c           jump in the third derivative across tau(2) is forced to
c           zero, thus the first and the second cubic polynomial pieces
c           are made to coincide.
c        ibcbeg = 1  means that the slope at tau(1) is made to equal
c           c2(1), supplied by input.
c        ibcbeg = 2  means that the second derivative at tau(1) is
c           made to equal c2(1), supplied by input.
c        ibcend = 0, 1, or 2 has analogous meaning concerning the
c           boundary condition at tau(n), with the additional infor-
c           mation taken from c2(n).
c
c     ********************** output ****************************
c
c     n, tau, c1, c2, ibcbeg, ibcend  are not altered by cubspl.
c     cj(i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
c        of the cubic interpolating spline with interior knots (or
c        joints) tau(2), ..., tau(n-1). precisely, in the interval
c        (tau(i), tau(i+1)), the spline f is given by
c           f(x) = c1(i)+h*(c2(i)+h*(c3(i)+h*c4(i)/3.)/2.)
c        where h = x - tau(i).
c     in other words, for i=1,...,n, c2(i) and c3(i) are respectively
c        equal to the values of the first and second derivatives of
c        the interpolating spline, and c4(i) is equal to the third deri-
c        vative of the interpolating spline in the interval (tau(i),
c        tau(i+1)). c4(n) is meaningless and is set to 0. for clarity.
c
c     **********************************************************
c
c      implicit double precision (a-h,o-z)
      implicit none
      integer n, ibcbeg, ibcend, i, l, m, j
      double precision tau, c1, c2, c3, c4, taum1, g, dtau, divdf1,
     :                 divdf3
      dimension tau(n),c1(n),c2(n),c3(n),c4(n)
c***** a tridiagonal linear system for the unknown slopes s(i) of
c  f at tau(i), i=1,...,n, is generated and then solved by gauss elim-
c  ination, with s(i) ending up in c2(i), all i.
c     c3(.) and c4(.) are used initially for temporary storage.
c
check -- n.ge.2
      if (n.lt.2) then
         write (6,111)
 111     format (/,' cubspl -- less than two pivots',/)
         stop
      endif
check -- tau strictly monotonous
      taum1=tau(2)
      if (tau(2)-tau(1)) 101,102,103
  101 if (n.eq.2) goto 200
      do 1 i=3,n
      if ((tau(i)-taum1).ge.0.d0) goto 102
    1 taum1=tau(i)
      goto 200
  102 write (6,222)
  222 format (/,' cubspl -- non monotonous abscissae',/)
      stop
  103 if (n.eq.2) goto 200
      do 3 i=3,n
      if ((tau(i)-taum1).le.0.d0) goto 102
    3 taum1=tau(i)
c
  200 l = n-1
compute first differences of tau sequence and store in c3(.). also,
compute first divided difference of data and store in c4(.).
      do 10 m=2,n
         c3(m) = tau(m) - tau(m-1)
   10    c4(m) = (c1(m) - c1(m-1))/c3(m)
construct first equation from the boundary condition, of the form
c             c4(1)*s(1) + c3(1)*s(2) = c2(1)
      if (ibcbeg-1)                     11,15,16
   11 if (n .gt. 2)                     goto 12
c     no condition at left end and n = 2.
      c4(1) = 1.d0
      c3(1) = 1.d0
      c2(1) = 2.d0*c4(2)
                                        goto 25
c     not-a-knot condition at left end and n .gt. 2.
   12 c4(1) = c3(3)
      c3(1) = c3(2) + c3(3)
      c2(1) = ((c3(2)+2.d0*c3(1))*c4(2)*c3(3)+c3(2)**2*c4(3))/c3(1)
                                        goto 19
c     slope prescribed at left end.
   15 c4(1) = 1.d0
      c3(1) = 0.d0
                                        goto 18
c     second derivative prescribed at left end.
   16 c4(1) = 2.d0
      c3(1) = 1.d0
      c2(1) = 3.d0*c4(2) - c3(2)/2.d0*c2(1)
   18 if(n .eq. 2)                      goto 25
c  if there are interior knots, generate the corresp. equations and car-
c  ry out the forward pass of gauss elimination, after which the m-th
c  equation reads    c4(m)*s(m) + c3(m)*s(m+1) = c2(m).
   19 do 20 m=2,l
         g = -c3(m+1)/c4(m-1)
         c2(m) = g*c2(m-1) + 3.d0*(c3(m)*c4(m+1)+c3(m+1)*c4(m))
   20    c4(m) = g*c3(m-1) + 2.d0*(c3(m) + c3(m+1))
construct last equation from the second boundary condition, of the form
c           (-g*c4(n-1))*s(n-1) + c4(n)*s(n) = c2(n)
c     if slope is prescribed at right end, one can go directly to back-
c     substitution, since c array happens to be set up just right for it
c     at this point.
      if (ibcend-1)                     21,30,24
   21 if (n .eq. 3 .and. ibcbeg .eq. 0) goto 22
c     not-a-knot and n .ge. 3, and either n.gt.3 or also not-a-knot at
c     left endpoint.
      g = c3(l) + c3(n)
      c2(n) = ((c3(n)+2.d0*g)*c4(n)*c3(l)
     *            + c3(n)**2*(c1(l)-c1(n-2))/c3(l))/g
      g = -g/c4(l)
      c4(n) = c3(l)
                                        goto 29
c     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
c     knot at left endpoint).
   22 c2(n) = 2.d0*c4(n)
      c4(n) = 1.d0
                                        goto 28
c     second derivative prescribed at right endpoint.
   24 c2(n) = 3.d0*c4(n) + c3(n)/2.d0*c2(n)
      c4(n) = 2.d0
                                        goto 28
   25 if (ibcend-1)                     26,30,24
   26 if (ibcbeg .gt. 0)                goto 22
c     not-a-knot at right endpoint and at left endpoint and n = 2.
      c2(n) = c4(n)
                                        goto 30
   28 g = -1.d0/c4(l)
complete forward pass of gauss elimination.
   29 c4(n) = g*c3(l) + c4(n)
      c2(n) = (g*c2(l) + c2(n))/c4(n)
carry out back substitution
   30 do 40 j=l,1,-1
   40    c2(j) = (c2(j) - c3(j)*c2(j+1))/c4(j)
c****** generate cubic coefficients in each interval, i.e., the deriv.s
c  at its left endpoint, from value and slope at its endpoints.
      do 50 i=2,n
         dtau = c3(i)
         divdf1 = (c1(i) - c1(i-1))/dtau
         divdf3 = c2(i-1) + c2(i) - 2.d0*divdf1
         c3(i-1) = 2.d0*(divdf1 - c2(i-1) - divdf3)/dtau
   50    c4(i-1) = (divdf3/dtau)*(6.d0/dtau)
c****** compute in addition c3(n). set c4(n) to 0.
         c3(n) = c3(l) + c4(l)*dtau
         c4(n) = 0.d0
                                        return
      end



      subroutine splget (n,x,t,k)
      integer n, k
      double precision x, t
      dimension x(n)
c
c***********************************************************************
c     *** the subroutine splget modifies the index k so that the
c     argument t lies within the interval | x(k) ... x(k+1) |.
c     in case of extrapolation, k is forced to the value 1 or n-1.
c
c     n      number of data points (n is assumed .ge. 2).
c     (x(i), i=1,...,n) abcissae of the points
c            (x is assumed to be strictly increasing).
c     t      argument for which the spline function is to be determined.
c     k      initial guess for k.
c
c                                       p. valiron  8-june-84
c***********************************************************************
c
      if(k.lt.1) k=1
      if(k.gt.n-1) k=n-1
      if(t.le.x(k+1)) go to 11
   10 if(k.eq.n-1) goto 20
      k=k+1
      if(t.gt.x(k+1)) go to 10
      go to 20
   11 if(k.eq.1) goto 20
      if(t.ge.x(k)) go to 20
      k=k-1
      go to 11
   20 return
      end 
	
		subroutine xyz(a1,b1,g1,a2,b2,g2,theta,phi,thetap,phip)
		implicit none
		real*8 o(3),h1(3),h2(3)
		real*8 safe_o(3),safe_h1(3),safe_h2(3)
		real*8 h3(3),h4(3)
		real*8 com_h2(3),r_com_h2,r_com_h,temp_com_h2(3)
		real*8 tempo(3),temph1(3),temph2(3)
		real*8 temph3(3),temph4(3)
		real*8 r_oh1,r_oh2,r_h1h2
		real*8 rr_oh1,rr_oh2,rr_h1h2
		real*8 temp1(3,3),temp2(3),temp3(3)
		real*8 theta,phi,temp00,temp01,temp02,temp03,temp04
		real*8 thetap,phip,temp10,temp11,temp12,temp13,temp14
		real*8,parameter :: r_w_h=6.00d0
		real*8,parameter :: r_h2=0.740d0
		real*8,parameter :: pi=4.0d0*datan(1.0d0)
		real*8 a1,b1,g1,rot_mat1(3,3),inv_rot_mat1(3,3)
		real*8 temp_rot_mat1(3,3)
		real*8 a2,b2,g2,rot_mat2(3,3)
		
		! coordinates of H2O
		o(1)=0.0d0
		o(2)=0.0d0
		o(3)=0.06530d0
		h1(1)=0.774880d0
		h1(2)=0.0d0
		h1(3)=-0.53470d0
		h2(1)=-0.774880d0
		h2(2)=0.0d0
		h2(3)=-0.53470d0
		
		safe_o=o
		safe_h1=h1
		safe_h2=h2
		
		! coordinates of hydrogens in H2
		h3(1)=0.0d0
		h3(2)=0.0d0
		h3(3)=r_h2*0.50d0
		
		h4(1)=0.0d0
		h4(2)=0.0d0
		h4(3)=-r_h2*0.50d0
		
		! find the COM of H2 molecule
		com_h2(1)=(h3(1)+h4(1))*0.50d0
		com_h2(2)=(h3(2)+h4(2))*0.50d0
		com_h2(3)=(h3(3)+h4(3))*0.50d0

		!find inter-atomic distance of H2O
		rr_oh1=dsqrt((o(1)-h1(1))**2+(o(2)-h1(2))**2+(o(3)-h1(3))**2)
		rr_oh2=dsqrt((o(1)-h2(1))**2+(o(2)-h2(2))**2+(o(3)-h2(3))**2)
		rr_h1h2=dsqrt((h2(1)-h1(1))**2+(h2(2)-h1(2))**2
     & +(h2(3)-h1(3))**2)
		
		!Euler rotation matrix for H2O
		call rot_mat(a1,b1,g1,rot_mat1)
		
		!Euler rotation matrix for H2
		call rot_mat(a2,b2,g2,rot_mat2)
		
		!SF rotation matrix for H2
		call sf_rot_mat(b2,a2,rot_mat2)
		
		tempo=o
		temph1=h1
		temph2=h2
		o=matmul(rot_mat1,tempo)
		h1=matmul(rot_mat1,temph1)
		h2=matmul(rot_mat1,temph2)
		
		!find inter-atomic distance of H2O after rotation
		r_oh1=dsqrt((o(1)-h1(1))**2+(o(2)-h1(2))**2+(o(3)-h1(3))**2)
		r_oh2=dsqrt((o(1)-h2(1))**2+(o(2)-h2(2))**2+(o(3)-h2(3))**2)
		r_h1h2=dsqrt((h2(1)-h1(1))**2+(h2(2)-h1(2))**2
     & +(h2(3)-h1(3))**2)
		
		if(abs(rr_oh1-r_oh1).gt.1d-9) then
		print*,"Interatomic distance wrong after rotation: 
     & O_H1",rr_oh1,r_oh1
		stop
		endif
		if(abs(rr_oh2-r_oh2).gt.1d-9) then
		print*,"Interatomic distance wrong after rotation: 
     & O_H2",rr_oh2,r_oh2
		stop
		endif
		if(abs(rr_h1h2-r_h1h2).gt.1d-9) then
		print*,"Interatomic distance wrong after rotation: 
     & H1_H2",rr_h1h2,r_h1h2
		stop
		endif

		call inv_rot_mat(a1,b1,g1,inv_rot_mat1)
	!	call inverse(rot_mat1,inv_rot_mat1,3)
	!	call rot_mat(-a1,-b1,-g1,inv_rot_mat1)
		tempo=o
		temph1=h1
		temph2=h2
		o=matmul(inv_rot_mat1,tempo)
		h1=matmul(inv_rot_mat1,temph1)
		h2=matmul(inv_rot_mat1,temph2)
		
		!find inter-atomic distance of H2O after inverse-rotation
		r_oh1=dsqrt((o(1)-h1(1))**2+(o(2)-h1(2))**2+(o(3)-h1(3))**2)
		r_oh2=dsqrt((o(1)-h2(1))**2+(o(2)-h2(2))**2+(o(3)-h2(3))**2)
		r_h1h2=dsqrt((h2(1)-h1(1))**2+(h2(2)-h1(2))**2
     & +(h2(3)-h1(3))**2)
		
		
		if(abs(rr_oh1-r_oh1).gt.1d-9) then
		print*,"Interatomic distance wrong after inverse-rotation: 
     & O_H1",rr_oh1,r_oh1
		stop
		endif
		if(abs(rr_oh2-r_oh2).gt.1d-9) then
		print*,"Interatomic distance wrong after inverse-rotation: 
     & O_H2",rr_oh2,r_oh2
		stop
		endif
		if(abs(rr_h1h2-r_h1h2).gt.1d-9) then
		print*,"Interatomic distance wrong after inverse-rotation: 
     & H1_H2",rr_h1h2,r_h1h2
		stop
		endif
		
		!check if H2O is shifted back to its reference position
		if(abs(safe_o(1)-o(1)).gt.1d-9) then
		print*,"Oxygen is not in its initial position along 
     & X-coordinate",safe_o(1),o(1)
		stop
		endif
		if(abs(safe_o(2)-o(2)).gt.1d-9) then
		print*,"Oxygen is not in its initial position along 
     & Y-coordinate",safe_o(2),o(2)
		stop
		endif
		if(abs(safe_o(3)-o(3)).gt.1d-9) then
		print*,"Oxygen is not in its initial position along 
     & Z-coordinate",safe_o(3),o(3)
		stop
		endif
		if(abs(safe_h1(1)-h1(1)).gt.1d-9) then
		print*,"Hydrogen1 is not in its initial position along 
     & X-coordinate",safe_h1(1),h1(1)
		stop
		endif
		if(abs(safe_h1(2)-h1(2)).gt.1d-9) then
		print*,"Hydrogen1 is not in its initial position along 
     & Y-coordinate",safe_h1(2),h1(2)
		stop
		endif
		if(abs(safe_h1(3)-h1(3)).gt.1d-9) then
		print*,"Hydrogen1 is not in its initial position along 
     & Z-coordinate",safe_h1(3),h1(3)
		stop
		endif
		if(abs(safe_h2(1)-h2(1)).gt.1d-9) then
		print*,"Hydrogen2 is not in its initial position along 
     & X-coordinate",safe_h2(1),h2(1)
		stop
		endif
		if(abs(safe_h2(2)-h2(2)).gt.1d-9) then
		print*,"Hydrogen2 is not in its initial position along 
     & Y-coordinate",safe_h2(2),h2(2)
		stop
		endif
		if(abs(safe_h2(3)-h2(3)).gt.1d-9) then
		print*,"Hydrogen2 is not in its initial position along 
     & Z-coordinate",safe_h2(3),h2(3)
		stop
		endif
		
		! rotation on H2 molecule
		temph3=h3
		temph4=h4
		
		h3=matmul(rot_mat2,temph3)
		h4=matmul(rot_mat2,temph4)
	!	h3(1)=dsqrt(temph3(1)**2+temph3(2)**2+temph3(3)**2)*dsin(b2)*dcos(a2)
	!	h3(2)=dsqrt(temph3(1)**2+temph3(2)**2+temph3(3)**2)*dsin(b2)*dsin(a2)
	!	h3(3)=dsqrt(temph3(1)**2+temph3(2)**2+temph3(3)**2)*dcos(b2)
	!	h4(1)=dsqrt(temph4(1)**2+temph4(2)**2+temph4(3)**2)*dsin(b2)*dcos(a2)
	!	h4(2)=dsqrt(temph4(1)**2+temph4(2)**2+temph4(3)**2)*dsin(b2)*dsin(a2)
	!	h4(3)=dsqrt(temph4(1)**2+temph4(2)**2+temph4(3)**2)*dcos(b2)
		
		! shifting of hydrogens of H2
		h3(3)=h3(3)+r_w_h
		h4(3)=h4(3)+r_w_h	
		com_h2(3)=com_h2(3)+r_w_h
		
		temph3=h3
		temph4=h4
		temp_com_h2=com_h2
		h3=matmul(inv_rot_mat1,temph3)
		h4=matmul(inv_rot_mat1,temph4)	
		com_h2=matmul(inv_rot_mat1,temp_com_h2)	
		
		r_com_h2=dsqrt(com_h2(1)**2+com_h2(2)**2+com_h2(3)**2)
		temp01=com_h2(3)/r_com_h2
		theta=dacos(temp01)
		
	!	if(com_h2(1).ge.0.00d0 .and. com_h2(2).ge.0.00d0) then
	!	temp03=com_h2(2)/com_h2(1)
	!	if(com_h2(2).eq.0.00d0 .and. com_h2(1).eq.0.00d0) temp03=0.d0
	!	phi=datan(temp03)
	!	else if(com_h2(1).lt.0.00d0 .and. com_h2(2).ge.0.00d0) then
	!	temp03=com_h2(2)/abs(com_h2(1))
	!	phi=-datan(temp03)+pi
	!	else if(com_h2(1).lt.0.00d0 .and. com_h2(2).lt.0.00d0) then
	!	temp03=abs(com_h2(2)/com_h2(1))
	!	phi=datan(temp03)+pi
	!	else if(com_h2(1).ge.0.00d0 .and. com_h2(2).lt.0.00d0) then
	!	temp03=abs(com_h2(2))/com_h2(1)
	!	phi=-datan(temp03)+pi*.00d0
	!	endif
		phi=datan2(com_h2(2),com_h2(1))
		if(phi.lt.0.00d0) phi=phi+2.0d0*pi
		
		if(theta.ne.theta) then
		print*,'Something is wrong in the transformation 
     & of angle theta',theta
		stop
		endif
		if(phi.ne.phi) then
		print*,'Something is wrong in the transformation 
     & of angle phi',phi
		stop
		endif
		
		! shifting COM of H2 and H2 molecule to origin
		h3=h3-com_h2
		h4=h4-com_h2
		
		r_com_h=dsqrt(h3(1)**2+h3(2)**2+h3(3)**2)
		temp11=h3(3)/r_com_h
		thetap=dacos(temp11)
		
	!	if(h3(1).ge.0.00d0 .and. h3(2).ge.0.00d0) then
	!	temp13=h3(2)/h3(1)
	!	if(h3(2).eq.0.00d0 .and. h3(1).eq.0.00d0) temp13=0.d0
	!	phip=datan(temp13)
	!	else if(h3(1).lt.0.00d0 .and. h3(2).ge.0.00d0) then
	!	temp13=h3(2)/abs(h3(1))
	!	phip=-datan(temp13)+pi
	!	else if(h3(1).lt.0.00d0 .and. h3(2).lt.0.00d0) then
	!	temp13=abs(h3(2)/h3(1))
	!	phip=datan(temp13)+pi
	!	else if(h3(1).ge.0.00d0 .and. h3(2).lt.0.00d0) then
	!	temp13=abs(h3(2))/h3(1)
	!	phip=-datan(temp13)+pi*2.00d0
	!	endif
		phip=datan2(h3(2),h3(1))
		if(phip.lt.0.00d0) phip=phip+2.0d0*pi
		
		if(thetap.ne.thetap) then
		print*,'Something is wrong in the transformation 
     & of angle theta_prime',thetap
		stop
		endif
		if(phip.ne.phip) then
		print*,'Something is wrong in the transformation 
     & of angle phi_prime',phip
		stop
		endif
		
!		print*,'theta_num=',theta*180d0/pi,'theta_ana=',dacos(cos(b1))*180d0/pi
!		print*,'phi_num=',phi*180d0/pi,'phi_ana=',datan2(dsin(b1)*dsin(g1),dsin(b1)*dcos(g1))*180d0/pi,dsin(b1)*dcos(g1)
!		print*,'theta"_num=',thetap*180d0/pi,'theta"_ana=',dacos(dsin(b1)*dcos(a2)*dsin(b2)+dcos(b1)*dcos(b2))*180d0/pi
!		print*,'phi"_num=',phip*180d0/pi,'phi"_ana=',datan2((dcos(b1)*dsin(g1)*dcos(a2)*dsin(b2)+dcos(g1)*dsin(a2)*dsin(b2)+dsin(b1)*dsin(g1)*dcos(b2)),(dcos(b1)*dcos(g1)*dcos(a2)*dsin(b2)+dsin(g1)*dsin(a2)*dsin(b2)-dcos(g1)*dsin(b1)*dcos(b2)))*180d0/pi
		end subroutine

		subroutine rot_mat(a,b,g,mat)
		implicit none
		real*8 a,b,g,mat(3,3)	!angles are supplied in radian
		
		mat(1,1)=dcos(a)*dcos(b)*dcos(g)-dsin(a)*dsin(g)
		mat(1,2)=-dcos(a)*dcos(b)*dsin(g)-dsin(a)*dcos(g)
		mat(1,3)=dcos(a)*dsin(b)
		
		mat(2,1)=dsin(a)*dcos(b)*dcos(g)+dcos(a)*dsin(g)
		mat(2,2)=-dsin(a)*dcos(b)*dsin(g)+dcos(a)*dcos(g)
		mat(2,3)=dsin(a)*dsin(b)
		
		mat(3,1)=-dsin(b)*dcos(g)
		mat(3,2)=dsin(b)*dsin(g)
		mat(3,3)=dcos(b)
		endsubroutine
		
		subroutine inv_rot_mat(a,b,g,mat)
		implicit none
		real*8 a,b,g,mat(3,3)	!angles are supplied in radian
		
		mat(1,1)=dcos(a)*dcos(b)*dcos(g)-dsin(a)*dsin(g)
		mat(2,1)=-dcos(a)*dcos(b)*dsin(g)-dsin(a)*dcos(g)
		mat(3,1)=dcos(a)*dsin(b)
		
		mat(1,2)=dsin(a)*dcos(b)*dcos(g)+dcos(a)*dsin(g)
		mat(2,2)=-dsin(a)*dcos(b)*dsin(g)+dcos(a)*dcos(g)
		mat(3,2)=dsin(a)*dsin(b)
		
		mat(1,3)=-dsin(b)*dcos(g)
		mat(2,3)=dsin(b)*dsin(g)
		mat(3,3)=dcos(b)
		endsubroutine

		subroutine sf_rot_mat(sf_t,sf_p,mat)
		implicit none
		real*8 sf_t,sf_p,mat(3,3)	!angles are supplied in radian
		
		mat(1,1)=dcos(sf_p)*dcos(sf_t)
		mat(1,2)=-dsin(sf_p)
		mat(1,3)=dcos(sf_p)*dsin(sf_t)
		
		mat(2,1)=dsin(sf_p)*dcos(sf_t)
		mat(2,2)=dcos(sf_p)
		mat(2,3)=dsin(sf_p)*dsin(sf_t)
		
		mat(3,1)=-dsin(sf_t)
		mat(3,2)=0.00d0
		mat(3,3)=dcos(sf_t)
		endsubroutine
		
		subroutine mat_vec(row,column,mat,vec,res)
		implicit none
		integer i,j
		integer row,column
		real*8 mat(row,column),vec(column),res(row)
		
		res=0.00d0
		do i=1,row
		res(i)=0.00d0
		do j=1,column
		res(i)=res(i)+mat(i,j)*vec(j)
		enddo
		enddo
		
		end subroutine
		
		subroutine mat_print(m,n,a)
		implicit none
		integer m,n,i,j
		real*8 a(m,n)
		
		do i=1,m
		do j=1,n
		if(abs(a(i,j)).le.1d-9 .and. a(i,j).ge.0.00d0) a(i,j)=0.00d0
		enddo
		enddo
		
		do i=1,m
		do j=1,n
		write(*,'(f12.5,2x)',advance='no')a(i,j)
		enddo
		write(*,*)' '
		enddo
		write(*,*)' '
		endsubroutine
		
		subroutine inverse(a,c,n)
		!============================================================
		! Inverse matrix
		! Method: Based on Doolittle LU factorization for Ax=b
		! Alex G. December 2009
		!-----------------------------------------------------------
		! input ...
		! a(n,n) - array of coefficients for matrix A
		! n      - dimension
		! output ...
		! c(n,n) - inverse matrix of A
		! comments ...
		! the original matrix a(n,n) will be destroyed 
		! during the calculation
		!===========================================================
		implicit none 
		integer n
		double precision a(n,n), c(n,n)
		double precision L(n,n), U(n,n), b(n), d(n), x(n)
		double precision coeff
		integer i, j, k

		! step 0: initialization for matrices L and U and b
		! Fortran 90/95 aloows such operations on matrices
		L=0.0
		U=0.0
		b=0.0

		! step 1: forward elimination
		do k=1, n-1
		   do i=k+1,n
			  coeff=a(i,k)/a(k,k)
			  L(i,k) = coeff
			  do j=k+1,n
				 a(i,j) = a(i,j)-coeff*a(k,j)
			  end do
		   end do
		end do

		! Step 2: prepare L and U matrices 
		! L matrix is a matrix of the elimination coefficient
		! + the diagonal elements are 1.0
		do i=1,n
		  L(i,i) = 1.0
		end do
		! U matrix is the upper triangular part of A
		do j=1,n
		  do i=1,j
			U(i,j) = a(i,j)
		  end do
		end do

		! Step 3: compute columns of the inverse matrix C
		do k=1,n
		  b(k)=1.0
		  d(1) = b(1)
		! Step 3a: Solve Ld=b using the forward substitution
		  do i=2,n
			d(i)=b(i)
			do j=1,i-1
			  d(i) = d(i) - L(i,j)*d(j)
			end do
		  end do
		! Step 3b: Solve Ux=d using the back substitution
		  x(n)=d(n)/U(n,n)
		  do i = n-1,1,-1
			x(i) = d(i)
			do j=n,i+1,-1
			  x(i)=x(i)-U(i,j)*x(j)
			end do
			x(i) = x(i)/u(i,i)
		  end do
		! Step 3c: fill the solutions x(n) into column k of C
		  do i=1,n
			c(i,k) = x(i)
		  end do
		  b(k)=0.0
		end do
		end subroutine inverse 


