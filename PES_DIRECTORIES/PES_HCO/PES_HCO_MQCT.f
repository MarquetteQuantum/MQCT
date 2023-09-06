      SUBROUTINE USER_DEFINED_PES(V,R,rvib,rvib2,alpha,beta,gamma,
     &										   aalpha,bbeta,ggamma)
!     INPUT:  R - distance betweenn COMs, rvib - vibrational coordinate
! 		      alpha,beta,gamma - Euler's angles of the first molecule
!   		  aalpha,bbeta, ggamma - Euler's angles of the second molecule
!     OUTPUT: V - value of the potential
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 V,R,alpha,beta,gamma,aalpha,bbeta,ggamma,t,rvib,rvib2
      REAL*8 R1,COG, v3d
	  logical ifirst/.true./

	  COG=dcos(beta)
!	  rvib = 2.1322d0
	  Call Vint3D(V,rvib,R,COG)
      return
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

! *********************************************************************
!                H-CO 3D potential
! *********************************************************************
! This program evaluates the 3D HCO potential V_HCO = V_int(3D) + V_CO
! V_int(3D) is a fit to our UCCSD(T)/CBS ab initio data, as described
! in the paper.
! V_CO is obtained by cubic spline interpolation of the CO(X) RKR1
! potential, A. Le Floch, Mol. Phys. 72, 133-144 (1991)
! Jacobi coordinates used:
!    r_CO  (C-O distance)
!    R     (distance from CO center of mass to H nucleus)
!    theta (angle between vector R and CO bond axis,
!           0 for linear C-O-H, 180 for linear H-C-O)
!
! For Vint3D:
!
! INPUT:
!      zsmallr: r_CO (a_0)
!      bigR   : R    (a_0)
!      xang   : cos(theta*pi/180)
! OUTPUT:
!      pot    : interaction energy between H and CO (hartree)
!
! V_HCO minimum is at R = 3.021 a_0, r_CO = 2.221 a_0, theta = 144.8 deg
! with V_HCO = V_int + V_CO = -3.06972964770039075d-02 hartree,
! V_HCO = 0 at R infinitely large and CO at its equilibrum bond length
!
! NOTE:
! The fit of V_int is valid for R >= 2 a_0
!             and 1.65 <= r_CO <= 2.85 a_0
!
! If R < 2 a_0 the potential is calculated at R = 2 a_0
!
! *********************************************************************
      subroutine Vint3D(pot,zsmallr,bigR, xang)
      implicit double precision(a-h, o-z)
      parameter(nsmallr_ab = 9, numpl=5,
     & shiftE = -3.06972964770039075d-02)
      double precision Vint2D_MAT(nsmallr_ab), w_MAT(nsmallr_ab),
     & wVint2D_MAT(nsmallr_ab), wPL_MAT(nsmallr_ab,numpl+1),
     & PL_new(6),
     & smallr_ab(nsmallr_ab), smallr_abSC(nsmallr_ab),
     & U_MAT(9,6), V_MAT(9,6), sig(6), tem_MAT(6), S_IMAT(6,6),
     & UT(6,9), UTP(6), SIUTP(6),COEF(6),V_MATeco(6,6)
      logical ifirst/.true./
c
      if (ifirst) then
         write(*,'(a)')
     +'=============================================================',
     +'UCCSD(T) H-CO interaction potential, April 2013.',
     +'Lei Song, Ad van der Avoird, Gerrit C. Groenenboom,',
     +'Theoretical Chemistry, IMM, Radboud University Nijmegen',
     +'http://www.theochem.ru.nl/, e-mail: gerritg@theochem.ru.nl.',
     +'============================================================='
c
c       Save the point
c
        Rsave = bigR
        zsave = xang
        r1save= zsmallr
c
c       ... and first check if this point gives the right value:
c
        pi=dacos(-1d0)
        bigR = 3.021d0
        xang = cos(144.8d0*pi/180d0)
        zsmallr = 2.221d0
        pot1 = -3.50062597810989751d-02
        go to 20
      endif
10    continue
      if (ifirst) then
        if (dabs((pot/pot1)-1d0).gt.1d-10) then
           write(*,*)'** ERROR in vint3D: integrity check not passed'
           write(*,*)'** v(',zsmallr,',',bigR,',',xang,')=',pot
           write(*,*)'** this should be ',pot1
           stop 230
        endif
c       OK, check passed, do the required point
        bigR = Rsave
        xang = zsave
        zsmallr = r1save
        ifirst=.false.
      endif
c
20    continue
      pot=0d0
c
      if (dabs(xang).gt.1d0) then
         write(*,*)'** ERROR in Vint3D: illegal cos(theta) = z = ',xang
         stop 231
      endif
c
      smallr_ab(1)  = 1.65d+00
      smallr_ab(2)  = 1.85d+00
      smallr_ab(3)  = 2.05d+00
      smallr_ab(4)  = 2.1322d+00
      smallr_ab(5)  = 2.20d+00
      smallr_ab(6)  = 2.25d+00
      smallr_ab(7)  = 2.45d+00
      smallr_ab(8)  = 2.65d+00
      smallr_ab(9)  = 2.85d+00
! change angle

!      xang=-xang1

! normalized zsmallr into -1 to 1 range
      zsmallrSC = ((zsmallr-1.64d+00)+(zsmallr-2.86d+00))/
     & (2.86d+00-1.64d+00)
      do i = 1, nsmallr_ab
         smallr_abSC(i) = ((smallr_ab(i)-1.64d+00)+
     &   (smallr_ab(i)-2.86d+00))/(2.86d+00-1.64d+00)
      enddo

! calculate Vint_2D

      do id_smallr = 1 , nsmallr_ab
         call Vint_2D(Vint2D_MAT(id_smallr),id_smallr,bigR,xang)
      enddo

! weighted matrix of Vint_2D, and weighted Vint_2D

      do id_smallr = 1,  nsmallr_ab
          if (Vint2D_MAT(id_smallr) .GT. 0.1d+00) then
            w_MAT(id_smallr) = exp((-Vint2D_MAT(id_smallr)+0.1d+00)
     &                         *4.00d+00)
          else
            w_MAT(id_smallr) = 1.0d+00
          endif
          wVint2D_MAT(id_smallr) = w_MAT(id_smallr)*
     &                             Vint2D_MAT(id_smallr)
      enddo

! weighted Legendre Polynomial
      do ids = 1, nsmallr_ab
         do id_npl = 1, numpl+1
            wPL_MAT(ids,id_npl) = pl(smallr_abSC(ids),id_npl-1)
     &                                  *w_MAT(ids)
         enddo
      enddo

      do id_npl =1, numpl+1
      PL_new(id_npl)  = pl(zsmallrSC, id_npl-1)
      enddo

! Vint_3D potential

      call svd(9,9,6,wPL_MAT,sig,.true.,U_MAT,.true.,V_MAT,ierr,tem_MAT)
! change sig(6) to diag_matrix
      do i = 1, 6
         do j= 1,6
            if (j .eq. i) then
               S_IMAT(i,j) = 1d+00/sig(i)
            else
               S_IMAT(i,j) = 0d+00
            endif
         enddo
      enddo
! simpify the V_MAT ,remove the term with value = 0

      V_MATeco = V_MAT(1:6,1:6)
      UT     = TRANSPOSE(U_MAT)
      UTP    = MATMUL(UT, wVint2D_MAT)
      SIUTP  = MATMUL(S_IMAT,UTP)
      COEF   = MATMUL(V_MATeco,SIUTP)

      Vint_3D = DOT_PRODUCT(PL_new,COEF)

! the CO POTENTIAL

!      call VCO_1D(VCO1D, zsmallr)

! add Vint_3d and VCO

!      pot  = VCO1D + Vint_3D-shiftE
      pot  = Vint_3D
!
      if (ifirst) go to 10
!
      END

! *********************************************************************
      subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
c
      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      double precision a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
      double precision c,f,g,h,s,x,y,z,scale,anorm
      double precision dsqrt,dmax1,dabs,dsign
      logical matu,matv
c
c     DOWNLOAD FROM: http://www.netlib.org/fmm/svd.f (double)
c     this subroutine is a translation of the algol procedure svd,
c     num. math. 14, 403-420(1970) by golub and reinsch.
c     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
c
c     this subroutine determines the singular value decomposition
c          t
c     a=usv  of a real m by n rectangular matrix.  householder
c     bidiagonalization and a variant of the qr algorithm are used.
c
c     INPUT:
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.  note that nm must be at least
c          as large as the maximum of m and n.
c
c        m is the number of rows of a (and u).
c
c        n is the number of columns of a (and u) and the order of v.
c
c        a contains the rectangular input matrix to be decomposed.
c
c        matu should be set to .true. if the u matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c        matv should be set to .true. if the v matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c     OUTPUT:
c
c        a is unaltered (unless overwritten by u or v).
c
c        w contains the n (non-negative) singular values of a (the
c          diagonal elements of s).  they are unordered.  if an
c          error exit is made, the singular values should be correct
c          for indices ierr+1,ierr+2,...,n.
c
c        u contains the matrix u (orthogonal column vectors) of the
c          decomposition if matu has been set to .true.  otherwise
c          u is used as a temporary array.  u may coincide with a.
c          if an error exit is made, the columns of u corresponding
c          to indices of correct singular values should be correct.
c
c        v contains the matrix v (orthogonal) of the decomposition if
c          matv has been set to .true.  otherwise v is not referenced.
c          v may also coincide with a if u is not needed.  if an error
c          exit is made, the columns of v corresponding to indices of
c          correct singular values should be correct.
c
c        ierr is set to
c          zero       for normal return,
c          k          if the k-th singular value has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.
c
c     this is a modified version of a routine from the eispack
c     collection by the nats project
c
c     modified to eliminate machep
c
      ierr = 0
c
      do 100 i = 1, m
c
         do 100 j = 1, n
            u(i,j) = a(i,j)
  100 continue
c     .......... householder reduction to bidiagonal form ..........
      g = 0.0d0
      scale = 0.0d0
      anorm = 0.0d0
c
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m) go to 210
c
         do 120 k = i, m
  120    scale = scale + dabs(u(k,i))
c
         if (scale .eq. 0.0d0) go to 210
c
         do 130 k = i, m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i)**2
  130    continue
c
         f = u(i,i)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,i) = f - g
         if (i .eq. n) go to 190
c
         do 150 j = l, n
            s = 0.0d0
c
            do 140 k = i, m
  140       s = s + u(k,i) * u(k,j)
c
            f = s / h
c
            do 150 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  150    continue
c
  190    do 200 k = i, m
  200    u(k,i) = scale * u(k,i)
c
  210    w(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m .or. i .eq. n) go to 290
c
         do 220 k = l, n
  220    scale = scale + dabs(u(i,k))
c
         if (scale .eq. 0.0d0) go to 290
c
         do 230 k = l, n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k)**2
  230    continue
c
         f = u(i,l)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,l) = f - g
c
         do 240 k = l, n
  240    rv1(k) = u(i,k) / h
c
         if (i .eq. m) go to 270
c
         do 260 j = l, m
            s = 0.0d0
c
            do 250 k = l, n
  250       s = s + u(j,k) * u(i,k)
c
            do 260 k = l, n
               u(j,k) = u(j,k) + s * rv1(k)
  260    continue
c
  270    do 280 k = l, n
  280    u(i,k) = scale * u(i,k)
c
  290    anorm = dmax1(anorm,dabs(w(i))+dabs(rv1(i)))
  300 continue
c     .......... accumulation of right-hand transformations ..........
      if (.not. matv) go to 410
c     .......... for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (g .eq. 0.0d0) go to 360
c
         do 320 j = l, n
c     .......... double division avoids possible underflow ..........
  320    v(j,i) = (u(i,j) / u(i,l)) / g
c
         do 350 j = l, n
            s = 0.0d0
c
            do 340 k = l, n
  340       s = s + u(i,k) * v(k,j)
c
            do 350 k = l, n
               v(k,j) = v(k,j) + s * v(k,i)
  350    continue
c
  360    do 380 j = l, n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
  380    continue
c
  390    v(i,i) = 1.0d0
         g = rv1(i)
         l = i
  400 continue
c     .......... accumulation of left-hand transformations ..........
  410 if (.not. matu) go to 510
c     ..........for i=min(m,n) step -1 until 1 do -- ..........
      mn = n
      if (m .lt. n) mn = m
c
      do 500 ii = 1, mn
         i = mn + 1 - ii
         l = i + 1
         g = w(i)
         if (i .eq. n) go to 430
c
         do 420 j = l, n
  420    u(i,j) = 0.0d0
c
  430    if (g .eq. 0.0d0) go to 475
         if (i .eq. mn) go to 460
c
         do 450 j = l, n
            s = 0.0d0
c
            do 440 k = l, m
  440       s = s + u(k,i) * u(k,j)
c     .......... double division avoids possible underflow ..........
            f = (s / u(i,i)) / g
c
            do 450 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  450    continue
c
  460    do 470 j = i, m
  470    u(j,i) = u(j,i) / g
c
         go to 490
c
  475    do 480 j = i, m
  480    u(j,i) = 0.0d0
c
  490    u(i,i) = u(i,i) + 1.0d0
  500 continue
c     .......... diagonalization of the bidiagonal form ..........
c     .......... for k=n step -1 until 1 do -- ..........
  510 do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
c     .......... test for splitting.
c                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            if (dabs(rv1(l)) + anorm .eq. anorm) go to 565
c     .......... rv1(1) is always zero, so there is no exit
c                through the bottom of the loop ..........
            if (dabs(w(l1)) + anorm .eq. anorm) go to 540
  530    continue
c     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0d0
         s = 1.0d0
c
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            if (dabs(f) + anorm .eq. anorm) go to 565
            g = w(i)
            h = dsqrt(f*f+g*g)
            w(i) = h
            c = g / h
            s = -f / h
            if (.not. matu) go to 560
c
            do 550 j = 1, m
               y = u(j,l1)
               z = u(j,i)
               u(j,l1) = y * c + z * s
               u(j,i) = -y * s + z * c
  550       continue
c
  560    continue
c     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
c     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0d0 * h * y)
         g = dsqrt(f*f+1.0d0)
         f = ((x - z) * (x + z) + h * (y / (f + dsign(g,f)) - h)) / x
c     .......... next qr transformation ..........
         c = 1.0d0
         s = 1.0d0
c
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = dsqrt(f*f+h*h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            if (.not. matv) go to 575
c
            do 570 j = 1, n
               x = v(j,i1)
               z = v(j,i)
               v(j,i1) = x * c + z * s
               v(j,i) = -x * s + z * c
  570       continue
c
  575       z = dsqrt(f*f+h*h)
            w(i1) = z
c     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0d0) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600
c
            do 590 j = 1, m
               y = u(j,i1)
               z = u(j,i)
               u(j,i1) = y * c + z * s
               u(j,i) = -y * s + z * c
  590       continue
c
  600    continue
c
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
         go to 520
c     .......... convergence ..........
  650    if (z .ge. 0.0d0) go to 700
c     .......... w(k) is made non-negative ..........
         w(k) = -z
         if (.not. matv) go to 700
c
         do 690 j = 1, n
  690    v(j,k) = -v(j,k)
c
  700 continue
c
      go to 1001
c     .......... set error -- no convergence to a
c                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end
! *********************************************************************
!                 VCO RKR 1st potential
! *********************************************************************
! Cubic spline interpolation of the CO RKR1 potential
! A. Le Floch, Mol. Phys. 72, 133-144 (1991), obtained with
! the RKR1 2.0 code of R. J. Le Roy, University of Waterloo
! Chem. Phys. Res. Rep. CP-657R (2004)
!
! NOTE: 1.58 <= r_CO <= 3.49 a_0
!
      SUBROUTINE VCO_1D(VCO1d, zsmallr)
! INPUT:
!      zsmallr: CO distance
! OUTPUT:
!      VCO1D: CO potential, the minimum is zero energy
!
      implicit double precision(a-h,o-z)
      parameter (n = 158, ZeroP = 4.9280667951504186d-03)
      double precision coR(n),VCO(n),b(n),c(n),d(n)
      logical first/.true./
c
      if (first) then
         write(*,'(a)')
     +'===============================================================',
     +'Cubic spline interpolation of the CO RKR1 potential',
     +'A. Le Floch, Mol. Phys. 72, 133-144 (1991), obtained with',
     +'the RKR1 2.0 code of R. J. Le Roy, University of Waterloo',
     +'Chem. Phys. Res. Rep. CP-657R (2004)',
     +' ',
     +'NOTE: 1.58 <= r_CO <= 3.49 a_0',
     +'==============================================================='
         first=.false.
      endif
c
      call load_COrkr(coR, VCO)
      call spline_gerrit(n,coR,VCO,b,c,d)
      VCO1d  = seval(n, zsmallr, coR, VCO, b, c, d)+ ZeroP
      END
! *********************************************************************
      subroutine spline_gerrit (n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
c
c  DOWNLOAD FROM: http://www.netlib.org/fmm/spline.f
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      double precision t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.d0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0d0
      c(n) = 0.0d0
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.0d0*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0d0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.0d0*c(i)
   40 continue
      c(n) = 3.0d0*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.0d0
      d(1) = 0.0d0
      b(2) = b(1)
      c(2) = 0.0d0
      d(2) = 0.0d0
      return
      end
c ---------------------------------------------------------------------
      double precision function seval(n, u, x, y, b, c, d)
      integer n
      double precision  u, x(n), y(n), b(n), c(n), d(n)
c
c  DOWNLOAD FROM: http://www.netlib.org/fmm/seval.f
c  this subroutine evaluates the cubic spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    n = the number of data points
c    u = the abscissa at which the spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
      integer i, j, k
      double precision dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
c
c  binary search
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
c
c  evaluate spline
c
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end
! *********************************************************************
! LOAD the CO RKR1 potential
!
      SUBROUTINE load_COrkr(coR, VCO)
      implicit double precision(a-h,o-z)
      common/COrkr/qqtab(316)
      double precision coR(158), VCO(158)

      do id = 1, 158
         coR(id) = qqtab(id)
         VCO(id) = qqtab(158+id)
      enddo

      END
! ---------------------------------------------------------------------
! DATA BLOCK
      BLOCK DATA DATA_COPES
      implicit double precision(a-h,o-z)
      common/COrkr/qqtab(316)
      DATA qqtab/
! CO distance, unit bohr
     &    1.56318145797966d+00,  1.58094488363852d+00,
     &    1.59870830929738d+00,  1.61647173495624d+00,
     &    1.63423516061510d+00,  1.65216025012155d+00,
     &    1.65410619101545d+00,  1.65608647072397d+00,
     &    1.65810186174322d+00,  1.66015316878529d+00,
     &    1.66224123057826d+00,  1.66436692177746d+00,
     &    1.66653115500035d+00,  1.66873488299881d+00,
     &    1.67097910098300d+00,  1.67326484911334d+00,
     &    1.67559321517793d+00,  1.67796533747495d+00,
     &    1.68038240792144d+00,  1.68284567541245d+00,
     &    1.68535644945669d+00,  1.68791610411861d+00,
     &    1.69052608230002d+00,  1.69318790039798d+00,
     &    1.69590315338122d+00,  1.69867352033226d+00,
     &    1.70150077050837d+00,  1.70438676998230d+00,
     &    1.70733348893180d+00,  1.71034300965673d+00,
     &    1.71341753541391d+00,  1.71655940017381d+00,
     &    1.71977107941812d+00,  1.72305520211659d+00,
     &    1.72641456404346d+00,  1.72985214261959d+00,
     &    1.73337111349890d+00,  1.73697486915400d+00,
     &    1.74066703976260d+00,  1.74445151675027d+00,
     &    1.74833247941261d+00,  1.75231442512145d+00,
     &    1.75640220372027d+00,  1.76060105683845d+00,
     &    1.76491666300948d+00,  1.76935518967116d+00,
     &    1.77392335337247d+00,  1.77862848982173d+00,
     &    1.78347863581083d+00,  1.78848262556572d+00,
     &    1.79365020474617d+00,  1.79899216620396d+00,
     &    1.80452051278804d+00,  1.81024865407269d+00,
     &    1.81619164604900d+00,  1.82236648580929d+00,
     &    1.82879247744459d+00,  1.83549169134464d+00,
     &    1.84248954774634d+00,  1.84981556817535d+00,
     &    1.85750435777050d+00,  1.86559691144181d+00,
     &    1.87414238451677d+00,  1.88320054691150d+00,
     &    1.89284527339419d+00,  1.90316965975594d+00,
     &    1.91429379764115d+00,  1.92637711830542d+00,
     &    1.93963908490212d+00,  1.95439637494429d+00,
     &    1.96759715091227d+00,  1.97478516881344d+00,
     &    1.98245019510032d+00,  1.99068046329780d+00,
     &    1.99959496735479d+00,  2.00936091825219d+00,
     &    2.02022645676947d+00,  2.03258897761561d+00,
     &    2.04716229970483d+00,  2.06550706219985d+00,
     &    2.09299352140611d+00,  2.13222199469650d+00,
     &    2.17349905623490d+00,  2.20508937810319d+00,
     &    2.22754694865605d+00,  2.24624206096582d+00,
     &    2.26273539114351d+00,  2.27774079523438d+00,
     &    2.29165570526266d+00,  2.30472829950874d+00,
     &    2.31712582664160d+00,  2.32896731848708d+00,
     &    2.34034104663734d+00,  2.36194113040569d+00,
     &    2.38725024220528d+00,  2.41112314073306d+00,
     &    2.43387712070690d+00,  2.45573227057232d+00,
     &    2.47684865517457d+00,  2.49734700827020d+00,
     &    2.51732107652513d+00,  2.53684539460178d+00,
     &    2.55598040160374d+00,  2.57477593163107d+00,
     &    2.59327366826273d+00,  2.61150891553392d+00,
     &    2.62951190444638d+00,  2.64730877566722d+00,
     &    2.66492233136810d+00,  2.68237261919547d+00,
     &    2.69967739201813d+00,  2.71685247429849d+00,
     &    2.73391205727900d+00,  2.75086893920456d+00,
     &    2.76773472261230d+00,  2.78451997773011d+00,
     &    2.80123437886159d+00,  2.81788681904725d+00,
     &    2.83448550711325d+00,  2.85103805033147d+00,
     &    2.86755152524376d+00,  2.88403253868625d+00,
     &    2.90048728065071d+00,  2.91692157030934d+00,
     &    2.93334089628342d+00,  2.94975045204304d+00,
     &    2.96615516716977d+00,  2.98255973509022d+00,
     &    2.99896863778753d+00,  3.01538616791675d+00,
     &    3.03181644868256d+00,  3.04826345178415d+00,
     &    3.06473101368548d+00,  3.08122285043286d+00,
     &    3.09774257121026d+00,  3.11429369079588d+00,
     &    3.13087964106352d+00,  3.14750378165170d+00,
     &    3.16416940990996d+00,  3.18087977021740d+00,
     &    3.19763806275806d+00,  3.21444745182804d+00,
     &    3.23131107374174d+00,  3.24823204439710d+00,
     &    3.26521346655475d+00,  3.28225843688081d+00,
     &    3.29937005279881d+00,  3.31655141919326d+00,
     &    3.33380565500402d+00,  3.35113589974905d+00,
     &    3.36854532001028d+00,  3.38603711591700d+00,
     &    3.40361452765913d+00,  3.42128084206281d+00,
     &    3.43903939925982d+00,  3.45689359948264d+00,
     &    3.47484691001665d+00,  3.49290287234249d+00,
! CO potential data, unit hartree
     &    4.31861952297040d-01,  3.93693428484224d-01,
     &    3.58279628075282d-01,  3.25421735389813d-01,
     &    2.94935283801373d-01,  2.66401295619077d-01,
     &    2.63431461362142d-01,  2.60434136543447d-01,
     &    2.57409304425866d-01,  2.54356946016885d-01,
     &    2.51277040176472d-01,  2.48169563721602d-01,
     &    2.45034491527570d-01,  2.41871796626246d-01,
     &    2.38681450301365d-01,  2.35463422181008d-01,
     &    2.32217680327351d-01,  2.28944191323807d-01,
     &    2.25642920359641d-01,  2.22313831312157d-01,
     &    2.18956886826535d-01,  2.15572048393390d-01,
     &    2.12159276424133d-01,  2.08718530324182d-01,
     &    2.05249768564092d-01,  2.01752948748661d-01,
     &    1.98228027684047d-01,  1.94674961442944d-01,
     &    1.91093705427862d-01,  1.87484214432539d-01,
     &    1.83846442701518d-01,  1.80180343987912d-01,
     &    1.76485871609386d-01,  1.72762978502377d-01,
     &    1.69011617274563d-01,  1.65231740255598d-01,
     &    1.61423299546137d-01,  1.57586247065141d-01,
     &    1.53720534595492d-01,  1.49826113827910d-01,
     &    1.45902936403189d-01,  1.41950953952747d-01,
     &    1.37970118137507d-01,  1.33960380685106d-01,
     &    1.29921693425429d-01,  1.25854008324493d-01,
     &    1.21757277516662d-01,  1.17631453335218d-01,
     &    1.13476488341280d-01,  1.09292335351086d-01,
     &    1.05078947461650d-01,  1.00836278074793d-01,
     &    9.65642809195738d-02,  9.22629100731290d-02,
     &    8.79321199799366d-02,  8.35718654695314d-02,
     &    7.91821017726896d-02,  7.47627845361129d-02,
     &    7.03138698356419d-02,  6.58353141880334d-02,
     &    6.13270745613401d-02,  5.67891083839373d-02,
     &    5.22213735522430d-02,  4.76238284371853d-02,
     &    4.29964318894760d-02,  3.83391432437527d-02,
     &    3.36519223216600d-02,  2.89347294339464d-02,
     &    2.41875253816568d-02,  1.94102714565114d-02,
     &    1.55668067090264d-02,  1.36378471377684d-02,
     &    1.17040662320129d-02,  9.76546160525316d-03,
     &    7.82203088344660d-03,  5.87377170499502d-03,
     &    3.92068172072316d-03,  1.96275859385680d-03,
     &    2.25514051876985d-17, -1.96759637288797d-03,
     &   -3.94003282451607d-03, -4.92806679515042d-03,
     &   -3.94003282451607d-03, -1.96759637288797d-03,
     &    2.25514051876985d-17,  1.96275859385680d-03,
     &    3.92068172072316d-03,  5.87377170499502d-03,
     &    7.82203088344660d-03,  9.76546160525316d-03,
     &    1.17040662320129d-02,  1.36378471377684d-02,
     &    1.55668067090264d-02,  1.94102714565114d-02,
     &    2.41875253816568d-02,  2.89347294339464d-02,
     &    3.36519223216600d-02,  3.83391432437527d-02,
     &    4.29964318894760d-02,  4.76238284371853d-02,
     &    5.22213735522430d-02,  5.67891083839373d-02,
     &    6.13270745613401d-02,  6.58353141880334d-02,
     &    7.03138698356419d-02,  7.47627845361129d-02,
     &    7.91821017726896d-02,  8.35718654695314d-02,
     &    8.79321199799366d-02,  9.22629100731290d-02,
     &    9.65642809195738d-02,  1.00836278074793d-01,
     &    1.05078947461650d-01,  1.09292335351086d-01,
     &    1.13476488341280d-01,  1.17631453335218d-01,
     &    1.21757277516662d-01,  1.25854008324493d-01,
     &    1.29921693425429d-01,  1.33960380685106d-01,
     &    1.37970118137507d-01,  1.41950953952747d-01,
     &    1.45902936403189d-01,  1.49826113827910d-01,
     &    1.53720534595492d-01,  1.57586247065141d-01,
     &    1.61423299546137d-01,  1.65231740255598d-01,
     &    1.69011617274563d-01,  1.72762978502377d-01,
     &    1.76485871609386d-01,  1.80180343987912d-01,
     &    1.83846442701518d-01,  1.87484214432539d-01,
     &    1.91093705427862d-01,  1.94674961442944d-01,
     &    1.98228027684047d-01,  2.01752948748661d-01,
     &    2.05249768564092d-01,  2.08718530324182d-01,
     &    2.12159276424133d-01,  2.15572048393390d-01,
     &    2.18956886826535d-01,  2.22313831312157d-01,
     &    2.25642920359641d-01,  2.28944191323807d-01,
     &    2.32217680327351d-01,  2.35463422181008d-01,
     &    2.38681450301365d-01,  2.41871796626246d-01,
     &    2.45034491527570d-01,  2.48169563721602d-01,
     &    2.51277040176472d-01,  2.54356946016885d-01,
     &    2.57409304425866d-01,  2.60434136543447d-01,
     &    2.63431461362142d-01,  2.66401295619077d-01/
      END
! *********************************************************************
!                 VCO RKR 1st potential ending
! *********************************************************************
! *********************************************************************
!                     Reproduce the 2D Vint
! *********************************************************************
      SUBROUTINE Vint_2D(Vint2d, id_smallr, bigR1, xang)

! INPUT:
!      id_smallr : 1, 2, 3, 4, 5, 6, 7, 8, 9 corresponding to
!      1.65, 1.85, 2.05, 2.1322, 2.20, 2.25, 2.45, 2.65, 2.85 bohr
!      bigR      : unit bohr
!      xang      : cosine of angle
! OUTPUT:
!      Vint_2D   : unit hartree, get Vint_2D in arbitary (bigR, xang)
!                  based on ab initio smallr_i grid

      implicit double precision(a-h,o-z)

      parameter (n_bigR_ii   = 32,
     &           n_ln_ii     = 13, !the order is l_aim = 0,1,2,3...12
     &           n_smallr_ii = 9)

      double precision  Rs(n_bigR_ii),
     &        CMAT_lr(n_ln_ii,n_smallr_ii),
     &        CMAT_RKHS(n_bigR_ii,n_ln_ii,n_smallr_ii),
     &        pes_Rn(n_ln_ii),
     &        Qvet(n_bigR_ii)

      integer mRKHS(n_ln_ii), lt_n(n_ln_ii), l_aim(n_ln_ii)

      common/maxqn/n_bigR_i, n_ln_i, n_smallr_i
      n_bigR_i   = n_bigR_ii
      n_ln_i     = n_ln_ii
      n_smallr_i = n_smallr_ii

! when bigR less than 2bohr, replace it into 2 bohr
      bigR = bigR1
      if (bigR1 .lt. 2.00d0) then
          bigR = 2.d0
!          write(*,*) 'Warning: bigR < 2, set to 2'
      endif

! load data for expansion coefficents

      call load_coeff(Rs,l_aim, lt_n, mRKHS, nRKHS, CMAT_lr, CMAT_RKHS)

      ! reproduce the Rn term
      Vint2d  = 0d+00
      do id_ln = 1, n_ln_ii
         ! long_range part
         re_lt       = CMAT_lr(id_ln,id_smallr)*bigR**(-lt_n(id_ln))
         call tt_df(re_tt,bigR,lt_n(id_ln))
         pes_Rnlr     = re_lt*re_tt
         ! short_range part
         pes_Rnsr     = 0d+00
         ! calculate the Qnm matrix
         call Qnm(Qvet, Rs, n_bigR_i, bigR, nRKHS, mRKHS(id_ln))
         do id_Rs  = 1, n_bigR_ii
            pes_Rnsr  = pes_Rnsr +
     &                 Qvet(id_Rs)*CMAT_RKHS(id_Rs,id_ln,id_smallr)
            pes_Rn(id_ln) = pes_Rnlr + pes_Rnsr
         enddo
         Vint2d = Vint2d + pl(xang,id_ln-1)*pes_Rn(id_ln)
      enddo
      end

! *********************************************************************
      function pl(x,n)
! Legendre Polynomial
! download from
! http://ww2.odu.edu/~agodunov/computing/programs/f90/Legendre.f90
! INPUT:
!      n: the number of Legendre Polynomial terms
!      x: the coordinate (-1<x<1)
! OUTPUT:
!     pl: Pl_0(x),Pl_1(x)... Pl_n-1(x)
!         p(0)= Pl_0(x), p(1)=Pl_1(x)...
! calculates Legendre polynomials Pn(x)
! using the recurrence relation
! if n > 100 the function retuns 0.0

      double precision pl
      double precision x
      double precision pln(0:n)
      integer n, k

      pln(0) = 1.0d0
      pln(1) = x

      if (n <= 1) then
        pl = pln(n)
      else
        do k=1,n-1
          pln(k+1) = ((2.0d0*k+1.0d0)*x*pln(k) - float(k)*pln(k-1))/
     &(float(k+1))
        enddo
       pl = pln(n)
      endif
      return
      END
!
! *********************************************************************
      SUBROUTINE tt_df(y, x, ip_lr)
! INPUT:
!      x
!      p_lr
! OUTPUT:
!      y
!     Tang-Toennies damping function

      IMPLICIT double precision (a-h, o-z)
      integer i, ip_lr, ifac

      a_tt    = 2.00d0      ! Calculated for H--CO, bbeta
      sum_val = 0
      do i = 0, ip_lr
        call fact(ifac,i)
        sum_val = sum_val + (a_tt*x)**i/ifac
      enddo

      y = 1.0d0-sum_val*exp(-a_tt*x)
      end
! ---------------------------------------------------------------------
      SUBROUTINE fact(fac,n)
!     Factorial

      IMPLICIT NONE

      integer fac, n, i

      fac = 1
      do i = 2, n
	fac = fac * i
      enddo

      end

! *********************************************************************
      SUBROUTINE Qnm(Qvec, R, n_R, X, nRKHS, mRKHS)

! set up the reproducing kenel qmn(R,X) vector in RKHS method
! [Eq.(17) in JCP vol 104, p 2584 (1995)]
! INPUT:
!      R:     (vector) the points used to produce the Hilbert space
!             the distance-like independent internal coordinate
!      n_R:   the number of R
!      X:     (scalar) the point which will be evalued in R coordinate
!      nRKHS: the order of smoothness
!      mRKHS: the reciprocal power of the weighting factor
! OUTPUT:
!      Qvet : (vector) the r.k. qmn (R,X)

      implicit double precision(a-h,o-z)
      parameter (lp = 10)
      double precision Qvec(n_R),
     &       R(n_R),
     &       P(0:lP)
      save

      if (lP .lt. nRKHS-1) then
         write(*,*) 'lP too small in subroutine qnm'
         stop
      endif

      call setup(P, nRKHS, mRKHS)
!     Return polynomial coefficients in P(0), P(1),...
!     P are coefficients of (2F_1coef)*nRKHS/binom(nRKHS+mRKHS,mRKHS)
!     without {(x<)/(x>)] and (x>)^[-(mRKHS+1)] terms

      do i = 1,n_R
           Qvec(i) = eval(X, R(i), P, nRKHS, mRKHS)
!          Return qnm(R,X)(vector)
      enddo

      end
! ---------------------------------------------------------------------
      double precision FUNCTION EVAL(X, R, P, n, m)
!     Return EVAL=qnm(R,X)(scalar) in x< / x>
!     R(scalar), X(scalar), P(vector)
      implicit double precision(a-h,o-z)

      double precision P(0:n)

      xlarge = max(X, R)
      xsmall = min(X, R)
      xsl    = xsmall/xlarge

      xx  = xsl
      Pol = P(0)
      do i = 1, n-1
          Pol = Pol + P(i)*xx
          xx  = xx*xsl
      enddo

      eval = Pol/xlarge**(m+1)

      end
! ---------------------------------------------------------------------
      SUBROUTINE SETUP(P, n, m)
! INPUT:
!      n : nRKHS
!      m : mRKHS
! OUTPUT:
!      P : (vector), P(0), P(1),...
!          coefficients of (2F_1coef)*nRKHS/binom(nRKHS+mRKHS,mRKHS)
      implicit double precision(a-h,o-z)
      parameter (kmax = 50)
      common /bin/  binom(0:kmax, 0:kmax), jmax
      double precision P(0:n)

!     Fill out BINOM by Pascal's triangle relation (if necessary)
      IF (jmax .EQ. 0) CALL Pascal

      if (n+m .ge. jmax) then
          write(*,'(a, i5)')
     &        'Increase kmax in subroutine Pascal to', n+m
          stop
      endif

      Call gauss_2f1(-n+1, m+1, n+m+1, P, n)
!     Get expansion coefficients of hypergeometric function {_2}F_1
      obin = dble(n)/binom(n+m, n)
!     n*(m!*n!)/(m+n)!
      do i = 0,n
        P(i) = P(i)*obin
      enddo

      end
! ---------------------------------------------------------------------
      SUBROUTINE GAUSS_2f1(i, j, k, P, n)
!     Hypergeometric function {_2}F_1,
!     see Abramowitz & Stegun 15.1.1,  15.4.1.

!     {_2}F_1 is a polynomial of degree n,
!     if either i or j is 0, -1, -2, ...

!     GAUSS_2F1 returns the expansion coefficients of
!     this polynomial in P(0..n).

      implicit double precision(a-h,o-z)
      double precision P(0:n)
      parameter (ihuge = 2**30)

      nmax = ihuge
      if (i .le. 0) then
        nmax = -i
      endif
      if (j .le. 0) then
         nmax = min(nmax,-j)
      endif

      if (nmax .eq. ihuge) then
        write(*,*) ' Either i or j must be 0,-1,-2,...'
      endif

      lfac=1
      do  l = 0, nmax
         Pi   = pochhammer(dble(i), l)
         Pj   = pochhammer(dble(j), l)
         Pk   = pochhammer(dble(k), l)
         p(l) = Pi*Pj/(dble(lfac)*Pk)
         lfac = lfac*(l+1)
      enddo

      end
! ---------------------------------------------------------------------
      double precision FUNCTION POCHHAMMER(z, n)
!     Pochhammer's Symbol, see Abramowitz & Stegun 6.1.22
!     POCHHAMMER(Z,N) returns P=(Z)_N.
!
!     Definition:  (Z)_0 = 1
!                  (Z)_N = Z*(Z+1)*...*(Z+N-1)

      implicit double precision(a-h,o-z)
      double precision P(0:n)

      P(0) = 1.d0
      do i = 1, n
         P(i) = P(i-1)*(z + dble(i-1))
      enddo
      pochhammer = P(n)

      end
! ---------------------------------------------------------------------
      SUBROUTINE PASCAL

!     Fill common bin with binomial coefficients (Pascal's triangle)

      implicit double precision(a-h,o-z)
      parameter (kmax = 50)
      common /bin/  binom(0:kmax, 0:kmax), jmax

      if ( kmax .eq. jmax ) return  ! bin has been set

!     Set jmax as sign that binom has been set:
      jmax = kmax

      binom(0,0) = 1.d0
      do  i=1,jmax
         binom(i,0) = 1.d0
         binom(i,i) = 1.d0
         do  j=1,i-1
            binom(i,j) = binom(i-1,j-1) + binom(i-1,j)
         enddo
      enddo

      end

! *********************************************************************
! Load RKHS parameters for expanxion coefficients of
! potential obtained from bigR fitting

      SUBROUTINE load_coeff(Rs,l_aim, lt_n, mRKHS, nRKHS,
     $ CMAT_lr, CMAT_RKHS)
      implicit double precision(a-h,o-z)
      integer q_tab(52)
      common/maxqn/n_bigR_i, n_ln_i, n_smallr_i
      common/coefs1/coe_lr(117)     ! long-range coefficients
      common/coefs2/coe_RKHS(3744)  ! RKHS coefficients
      common/qnrs/q_tab
      common/grid/Rss(32)

      integer mRKHS(13), l_aim(13), lt_n(13)
      double precision CMAT_lr(n_ln_i, n_smallr_i),
     &       CMAT_RKHS(n_bigR_i,n_ln_i, n_smallr_i),
     &       Rs(n_bigR_i)

      nRKHS = 2
      do iid=1,n_bigR_i
         Rs(iid)=Rss(iid)
      enddo
      do iln = 1 , n_ln_i
         mRKHS(iln)  = q_tab(4*iln)
         l_aim(iln)  = q_tab(4*iln-2)
         lt_n(iln)   = q_tab(4*iln-1)
         do ismallr = 1 , n_smallr_i
            CMAT_lr(iln,ismallr)  = coe_lr(iln+(ismallr-1)*13)
            do ibigR = 1 , n_bigR_i
               CMAT_RKHS(ibigR,iln,ismallr) =
     &                    coe_RKHS(ibigR+(iln-1)*32+(ismallr-1)*32*13)
            enddo
         enddo
      enddo

      end
! ---------------------------------------------------------------------
! DATA BLOCK: /bin/, initio jmax
      BLOCK DATA BLPAS
      implicit double precision(a-h, o-z)
      PARAMETER (kmax=50)
      common /bin/ binom(0:kmax, 0:kmax), jmax
      data jmax/0/
      END
! ---------------------------------------------------------------------
! DATA BLOCK:
      BLOCK DATA DATA_POT
      implicit double precision(a-h,r-z)

      integer q_tab(52)

      common/maxqn/n_bigR_i, n_ln_i, n_smallr_i
      common/coefs1/coe_lr(117)     ! long-range coefficients
      common/coefs2/coe_RKHS(3744)  ! RKHS coefficients
      common/qnrs/q_tab
      common/grid/Rss(32)
! ---------------------------------------------------------------------
! the ab initio grid, bigR data, unit bohr
      DATA Rss/
     &   2.0000d+00,
     &   2.2000d+00,
     &   2.4000d+00,
     &   2.6000d+00,
     &   2.8000d+00,
     &   3.0000d+00,
     &   3.2000d+00,
     &   3.4000d+00,
     &   3.6000d+00,
     &   3.8000d+00,
     &   4.0000d+00,
     &   4.2000d+00,
     &   4.4000d+00,
     &   4.6000d+00,
     &   4.8000d+00,
     &   5.0000d+00,
     &   5.2000d+00,
     &   5.4000d+00,
     &   5.6000d+00,
     &   5.8000d+00,
     &   6.0000d+00,
     &   6.4000d+00,
     &   6.8000d+00,
     &   7.2000d+00,
     &   7.6000d+00,
     &   8.0000d+00,
     &   9.3199d+00,
     &   10.8577d+00,
     &   12.6491d+00,
     &   14.7361d+00,
     &   17.1675d+00,
     &   20.0000d+00/
! ---------------------------------------------------------------------
! Quantum number and related parameters (id_ln,l_aim(l),lt_n(n),mRKHS)
      DATA q_tab/
     &    1,  0,  6,  7,
     &    2,  1,  7,  8,
     &    3,  2,  6,  7,
     &    4,  3,  7,  8,
     &    5,  4,  8,  9,
     &    6,  5,  9, 10,
     &    7,  6, 10,  9,
     &    8,  7, 11, 10,
     &    9,  8, 12, 11,
     &   10,  9, 13, 12,
     &   11, 10, 14, 13,
     &   12, 11, 15, 14,
     &   13, 12, 16, 15/
! ---------------------------------------------------------------------
! The expansion coefficients for the analytic longe-range fits
      DATA coe_lr/
! id_smallr = 1
! smallr = 1.65 bohr, order is n6l0 n7l1 n6l2 n7l3 n8l4 ... n16l12
     &    -1.86215448965293d+01,     8.12002668755660d+01,
     &    -6.28303246398490d-01,     1.27329153579483d+00,
     &    -4.35815840478602d+00,     6.11108032937551d+01,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,
! id_smallr = 2, smallr = 1.85 bohr
     &    -1.95843781466434d+01,     7.53313175179620d+01,
     &    -7.55234216532640d-01,    -6.97182052154721d-01,
     &    -3.17868174474587d-01,    -2.50350649754062d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,
! id_smallr = 3, smallr = 2.05 bohr
     &    -2.09037403790486d+01,     7.03976596342499d+01,
     &    -1.06059832127932d+00,    -2.72797375121022d+00,
     &    -3.89981985119670d+00,    -2.90031221740528d+02,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,
! id_smallr = 4, smallr = 2.1322 bohr
     &    -2.12201496820663d+01,     7.03027053745700d+01,
     &    -1.42621172122518d+00,    -4.75132290767043d+00,
     &     2.20728522320447d+01,    -2.30874534893996d+02,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,
! id_smallr = 5, smallr = 2.20 bohr
     &    -2.18317783255118d+01,     6.01838031242214d+01,
     &    -1.69050164717435d+00,    -3.94277706975864d+00,
     &     4.50934487392890d+01,     1.26790952683031d+03,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,
! id_smallr = 6, smallr = 2.25 bohr
     &    -2.20936866931153d+01,     6.50057356832982d+01,
     &    -1.63455122117228d+00,    -4.60369621308745d+00,
     &     2.02893918711756d+01,    -1.88785229959318d+02,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,
! id_smallr = 7, smallr = 2.45 bohr
     &    -2.38077043779279d+01,     6.35249845552372d+01,
     &    -2.58206289561288d+00,    -4.47302434991739d+00,
     &    -3.20559368988716d+00,    -2.47746224682310d+02,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,
! id_smallr = 8, smallr = 2.65 bohr
     &    -2.53251228458261d+01,     6.42843942034100d+01,
     &    -3.50240139546976d+00,    -2.06187578173496d+00,
     &     5.56070227345916d-01,     4.43605084624622d+01,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,
! id_smallr = 9, smallr = 2.85 bohr
     &    -2.55662461933946d+01,     6.80023603358404d+01,
     &    -4.61577025194795d+00,    -9.03109787449004d-01,
     &     9.14573428963305d+00,     8.50542704064894d+01,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00,     0.00000000000000d+00,
     &     0.00000000000000d+00/
! ---------------------------------------------------------------------
      DATA coe_RKHS(1:1000)/
! id_smallr =1
! smallr = 1.65 bohr, order is n6l0 n7l1 n6l2 n7l3 n8l4 ... n16l12
     &     5.48931736051164d+03,    -4.17543803842683d+03,
     &     3.15215829590817d+03,     1.02193295152996d+03,
     &     2.00921841634506d+03,     1.16101844482344d+03,
     &     4.83768353786980d+03,     7.90598149260684d+03,
     &     8.72151612627623d+03,     9.81094289430708d+03,
     &     1.01400569714898d+04,     9.90052246233053d+03,
     &     8.69900335604214d+03,     8.25197448870962d+03,
     &     3.06778981625441d+03,     4.39528306137449d+03,
     &    -1.89485990578225d+03,    -6.76970690649785d+03,
     &    -8.55707321595167d+03,    -2.21821404110379d+03,
     &    -3.09702374361717d+04,    -3.43289797186881d+04,
     &    -3.64642721470450d+04,    -3.05560959007333d+04,
     &    -2.12125605461722d+04,    -1.44594051756717d+04,
     &     2.89772699623958d+04,     4.60696031777950d+04,
     &     1.55216016756357d+04,     1.48393353439313d+03,
     &    -3.03591253946435d+03,    -3.81245546445855d+03,
     &    -1.10208324496196d+04,     7.52456581417157d+03,
     &    -6.14364559427528d+03,     1.00602260860078d+04,
     &     1.97194047577366d+04,     3.55906793713611d+04,
     &     1.56204448511048d+04,    -4.85175250415953d+03,
     &    -4.94878099184585d+03,    -1.31077246301801d+04,
     &    -2.06777140170489d+04,    -3.42322436727361d+04,
     &    -4.11630676484042d+04,    -5.85691673205617d+04,
     &    -7.47959553735090d+04,    -7.94423707306840d+04,
     &    -9.91271252560566d+04,    -5.71696876934581d+04,
     &    -1.57869465222174d+05,     9.99086228975717d+04,
     &    -1.77849734460893d+05,     8.55328153455472d+04,
     &     1.30130476347595d+05,     2.61802276499590d+05,
     &     1.39545568508006d+05,     4.87703610028195d+05,
     &     2.60828591749968d+05,    -3.40490885775940d+05,
     &    -2.48975367739521d+05,    -8.32776516906848d+04,
     &     1.42220240535092d+04,     6.49121824702933d+03,
     &     1.05178729339123d+04,    -1.19974156228032d+04,
     &     2.91155150730545d+03,    -2.63049965922817d+03,
     &    -3.42147953926019d+03,    -9.47071580200798d+03,
     &     2.58303344244484d+03,     1.05891821680745d+04,
     &     7.80195229643723d+03,     8.75689087799203d+03,
     &     8.86674328705007d+03,     9.21625681395010d+03,
     &     1.04690334935088d+04,     6.80031682621972d+03,
     &     1.40959958342214d+04,     3.23123155302915d+03,
     &     7.99755634951078d+03,     7.78427125144900d+03,
     &     1.67925830693973d+03,    -2.35713497791938d+04,
     &     9.98684805854456d+03,    -3.83768993209611d+04,
     &    -3.06413019181720d+04,    -3.95314179397660d+04,
     &    -2.50280881939015d+04,    -3.31667577024155d+04,
     &     1.37936958752203d+04,     5.08448757678873d+04,
     &     1.87542710488061d+04,     1.52571262816081d+03,
     &    -7.78645287237081d+02,    -7.63587934932224d+03,
     &    -8.15716348782588d+03,     2.69704695062012d+03,
     &    -8.86288558856988d+03,     4.64091735759222d+03,
     &     1.69453116386197d+04,     3.76218274239139d+04,
     &    -6.81143620952712d+03,    -3.89389585936256d+04,
     &    -1.75883594864467d+04,    -1.38700223809767d+04,
     &    -8.55974897225366d+03,     1.66333480353327d+02,
     &    -5.87433283265676d+03,    -5.58855891057998d+03,
     &    -4.50091094286210d+02,    -3.09681042760395d+04,
     &     1.09370439136072d+04,    -8.06671431211520d+04,
     &     3.00169159623988d+04,    -4.15167449370785d+04,
     &    -2.22367859761419d+04,     4.07323399600286d+03,
     &     3.44968627463326d+04,     8.79239809759889d+04,
     &     4.85795190251338d+04,     1.73841001447709d+05,
     &     7.89683343728781d+04,    -1.21212723999124d+05,
     &    -8.25193393966622d+04,    -9.30396843929712d+03,
     &    -3.50807714370255d+04,     4.12363719734584d+04,
     &     2.64115626518565d+04,    -2.38435067251875d+04,
     &     1.12966160258493d+04,     5.71052248022450d+03,
     &    -8.36003841771784d+03,    -8.26844189961896d+04,
     &     7.65578189684018d+04,     2.04111981454787d+05,
     &     1.07407242956413d+05,     8.61369957955674d+04,
     &     4.57346660301085d+04,     5.33776998844021d+04,
     &    -4.25848769088246d+04,     7.46301000574125d+04,
     &    -1.37539895490310d+05,     1.89994801543483d+05,
     &    -7.42091726059165d+04,    -7.05189145944927d+04,
     &     5.69453959036696d+04,     5.88604919380630d+05,
     &    -6.20165878790539d+05,     2.48053625808502d+05,
     &    -1.86261059583681d+05,    -2.83453694342704d+05,
     &    -1.46968682020546d+05,    -6.24862901466053d+05,
     &    -4.57026983824238d+05,     4.99865823162220d+05,
     &     3.22098877323305d+05,     1.69438576909205d+05,
     &    -1.06473484204171d+05,     3.93661398271473d+05,
     &    -7.70931489461168d+03,     2.84428529262724d+03,
     &    -2.33663708932545d+04,    -4.96317432337230d+03,
     &     7.25181939567487d+04,     2.88436245099372d+05,
     &    -2.68521587492078d+05,    -7.80185685093580d+05,
     &    -3.82544634595887d+05,    -3.07496767624776d+05,
     &    -1.04515222544528d+05,    -1.69885550340671d+04,
     &    -1.89672031714589d+04,     4.43858819172949d+05,
     &    -1.39080220011918d+05,     8.49491174561211d+05,
     &    -1.28927129711125d+06,     3.20111994198274d+06,
     &    -2.45889202645911d+06,    -7.78847350200724d+05,
     &     2.05894155274159d+06,    -1.35055208019653d+06,
     &     5.92332394363318d+05,    -4.81906930668558d+05,
     &     3.43839353814544d+05,     8.26431394821775d+05,
     &     1.22151440414879d+06,    -7.49964432771027d+05,
     &     4.49427689512918d+05,    -5.09851920150993d+06,
     &     1.15381116668184d+07,    -1.39232561773813d+07,
     &     1.13907474258771d+04,    -1.67768186852881d+04,
     &     5.15366447374541d+03,     6.65659802187780d+02,
     &    -1.14548759338712d+04,    -8.15718167025155d+04,
     &     1.03932178983036d+05,     2.08490106466502d+05,
     &     4.50683298945152d+04,     1.44971265754257d+03,
     &    -6.86322629464708d+03,    -7.19153842576606d+04,
     &    -7.29328301710832d+03,    -1.06228837167355d+05,
     &     7.31592506911718d+04,    -2.44383973867478d+05,
     &     1.62059814831483d+05,    -1.00652942588887d+05,
     &    -2.10365694034773d+04,    -2.10340038472740d+05,
     &     2.82069089437395d+05,    -1.13110513735792d+05,
     &     9.41929813243884d+04,    -1.38127888943775d+05,
     &     1.92468252105072d+05,    -1.29165631982796d+05,
     &     5.20861813882086d+04,    -1.60926909822293d+04,
     &     3.34121136217180d+04,     8.48229254747855d+04,
     &    -1.40471975320751d+05,    -5.86844705745131d+04,
     &    -6.21864824559988d+03,    -4.68036756186719d+03,
     &    -2.21111777801620d+03,     4.25935924472745d+03,
     &     7.16307117473997d+04,     2.77784536570096d+05,
     &    -2.83950549278522d+05,    -6.71392805353761d+05,
     &    -1.03563710754417d+05,    -9.08782390021902d+03,
     &     1.44235885117081d+05,    -1.14008089369433d+05,
     &     5.17413079583357d+05,    -3.72356525458564d+05,
     &     6.38788182715885d+05,    -6.12842078931008d+05,
     &     2.10670550748174d+06,    -4.31646022906792d+06,
     &     4.31331596377265d+06,     2.24845200347468d+05,
     &    -2.36958958492427d+06,     1.20192110218022d+06,
     &    -5.67299463926153d+05,     9.24065675772925d+04,
     &     9.61491299235529d+04,    -1.44993134137724d+05,
     &    -6.93369000837597d+04,     2.56380078337772d+04,
     &    -1.70125910956194d+05,     2.74456248594340d+05,
     &    -4.96886109754463d+05,     3.09074167120576d+06,
     &     1.56042153175801d+04,    -2.77402179465397d+04,
     &     2.15460006947751d+04,     1.69136604843507d+04,
     &    -9.37490344316795d+04,    -7.88725471705708d+05,
     &     7.31788792817884d+05,     2.03594403953832d+06,
     &     2.00143797969998d+05,    -2.52376788143137d+05,
     &    -2.68830881178340d+05,    -4.89393381415663d+05,
     &     1.15283499858587d+04,    -1.40330444550705d+05,
     &    -1.95022846488769d+06,     4.41991728469832d+06,
     &    -9.87826327750388d+06,     1.33526742512662d+07,
     &    -1.01393286178746d+07,     2.19430232299834d+06,
     &     5.22434723427516d+05,    -6.47280454861538d+05,
     &     4.39973384784958d+05,     5.75190756605746d+05,
     &     5.83091157064112d+05,    -1.62289201907249d+06,
     &     1.62765178177744d+06,    -1.31902094814025d+06,
     &     7.45715526000190d+06,    -4.47302754741117d+07,
     &     1.38280218749993d+08,    -1.66432908180505d+08,
     &     6.87045298190631d+03,    -4.62075563741187d+04,
     &    -1.26938619873567d+04,    -8.35552613925572d+04,
     &     3.07030958674713d+05,     2.11308261146132d+06,
     &    -1.81618641417549d+06,    -6.01976566562472d+06,
     &    -3.58015441068382d+05,     1.48392446879004d+06,
     &     2.77215648499047d+05,     4.77677238499489d+06,
     &    -7.19129393126712d+06,     1.33666755262331d+07,
     &    -1.87286307637704d+07,     2.83813840679797d+07,
     &    -6.24995785040826d+07,     1.42066202722493d+08,
     &    -1.52973263592058d+08,     1.78173196767188d+07,
     &     6.79783352170583d+07,    -3.60588335817472d+07,
     &     1.42692168167400d+07,    -2.25032995684289d+07,
     &     3.66718882285421d+07,    -2.43936698256192d+07,
     &     4.70150008811313d+06,     1.48773424267511d+07,
     &    -1.64019211973635d+08,     9.78454573562111d+08,
     &    -2.77633426085565d+09,     1.41749605248058d+09,
     &     5.48037889914584d+04,    -1.51407760518661d+05,
     &     1.07440708607986d+05,     1.08613176005324d+05,
     &    -2.33528702376700d+05,    -5.01629386439778d+06,
     &     4.50875608124226d+06,     1.66044083227464d+07,
     &     1.09474058322559d+06,    -3.64180131305125d+06,
     &    -1.04261220465737d+07,     7.86603752862179d+06,
     &    -1.44158031910054d+07,    -1.92166769229953d+07,
     &     9.30576699223970d+07,    -2.13637759354028d+08,
     &     4.79168721649359d+08,    -8.78918322384193d+08,
     &     8.30990615931845d+08,    -1.16209935141023d+08,
     &    -3.06558196286374d+08,     2.33057645929647d+08,
     &    -3.17508297032972d+08,     7.22544564321232d+08,
     &    -1.09887424014644d+09,     7.30274750594245d+08,
     &    -2.46295629166137d+08,     2.52907638051176d+08,
     &    -2.18256180798257d+08,    -3.26541220934402d+09,
     &     3.01345336292347d+10,    -8.42387369244792d+10,
     &     3.54227886494833d+04,    -2.02985878296531d+05,
     &     1.09052472510328d+05,    -3.13096507181262d+05,
     &     9.97515269279548d+05,     1.09298663644472d+07,
     &    -9.89818368621505d+06,    -4.16788000050709d+07,
     &    -5.05301303444534d+06,     1.75345010936569d+07,
     &     1.00720758491408d+06,     2.67474051885747d+07,
     &     3.07691668324441d+06,     2.12858699342467d+07,
     &    -1.23547364178636d+07,    -1.72350445118057d+08,
     &     9.98713902285902d+08,    -2.77507193109638d+09,
     &     3.39797630877783d+09,    -9.89871384146691d+08,
     &    -8.39013322815515d+08,     3.13021908317587d+08,
     &     1.31880335197338d+08,     2.05756461735044d+08,
     &    -3.83605366593175d+08,    -2.78064740352485d+07,
     &     3.02240740136327d+08,    -2.92444739372760d+09,
     &     4.10772769086199d+10,    -3.41664711965578d+11,
     &     1.64264529816935d+12,    -3.46919747678655d+12,
     &     5.96118764305352d+04,    -1.87781551955929d+05,
     &     1.97243600624297d+05,    -9.68493180209492d+04,
     &    -5.71045486745387d+05,    -1.85151031481794d+07,
     &     1.74835117983421d+07,     7.36139559914519d+07,
     &     2.23655081435217d+07,    -2.88931453091609d+07,
     &    -4.59258578801645d+06,    -2.16665342128242d+08,
     &     5.21872925153146d+08,    -5.48241597110049d+08,
     &    -7.66413760036389d+08,     4.79844828693621d+09,
     &    -1.66504513072353d+10,     3.69152457283682d+10,
     &    -3.99486949042720d+10,     9.70931305824867d+09,
     &     1.08808084115969d+10,    -4.37706369821226d+09,
     &    -5.61032516729369d+09,     1.56907041783528d+10,
     &    -2.48714626829258d+10,     1.94961712122761d+10,
     &    -1.07203132217864d+10,     3.12081027628633d+10,
     &    -5.04980273306594d+11,     6.01287228050618d+12,
     &    -3.43571201045907d+13,     7.08723201767870d+13, ! id_smallr = 2, smallr = 1.85 bohr
     &     5.89886682815496d+03,    -5.15397004139595d+03,
     &     2.98422644219024d+03,     1.90213596659279d+03,
     &     1.90012134168367d+03,     9.30264299418906d+02,
     &    -1.12917151100738d+03,     3.29671132698473d+03,
     &     1.36814497701109d+04,     1.37572545084827d+04,
     &     1.28217378251229d+04,     1.33831813277906d+04,
     &     1.01619483832925d+04,     1.09675843428351d+04,
     &     4.28731865341115d+03,     3.30638242878471d+03,
     &     2.99348460308781d+03,    -1.40624001516600d+04,
     &     8.51394644089388d+03,    -3.47794868374132d+04,
     &    -7.74704909974438d+03,    -4.69082892673807d+04,
     &    -3.71393459297390d+04,    -3.56429064672700d+04,
     &    -2.30919699655383d+04,    -1.68150484739942d+04,
     &     3.18125078378120d+04,     5.06183232792370d+04,
     &     1.66235341543603d+04,     6.88006344432119d+02,
     &     7.44640400869021d+01,    -9.35160377885024d+03,
     &    -1.24714854569634d+04,     3.02744678764073d+03,
     &    -1.06872405941078d+04,     1.52040536716257d+04,
     &     1.59427181562678d+04,     4.14931668549463d+04,
     &     7.35589997732270d+04,     5.08366790488928d+04,
     &    -4.34347186768207d+04,    -4.24201293316948d+04,
     &    -3.88495535253623d+04,    -5.26959267793558d+04,
     &    -6.60429830142279d+04,    -6.16387414011749d+04,
     &    -1.01108173490219d+05,    -7.96122235230217d+04,
     &    -1.07557646453037d+05,    -6.14691328042281d+04,
     &    -2.18489488370184d+05,     2.26813040564395d+05,
     &    -2.88934277255655d+05,     1.36950184914566d+05,
     &     1.37211181198556d+05,     2.86580696374779d+05,
     &     1.77751639000850d+05,     5.13440033325586d+05,
     &     2.74822532385361d+05,    -3.77920672124497d+05,
     &    -2.67178894897806d+05,    -7.10488358342224d+04,
     &    -2.33038380832311d+03,     1.78036599525552d+04,
     &     1.64100511821015d+04,    -2.02829272223377d+04,
     &     4.64479061435556d+03,     9.91982050094934d+02,
     &    -1.81671936590517d+03,    -7.60990024749507d+03,
     &    -2.08960134632880d+04,    -6.58913834304958d+03,
     &     3.07957921953239d+04,     2.03857950000696d+04,
     &     1.40177481282455d+04,     1.42336082568915d+04,
     &     1.31461388576859d+04,     1.25789593251766d+04,
     &     7.00140212252084d+03,     1.88495154788235d+04,
     &    -8.96872299314697d+03,     2.54054816499963d+04,
     &    -2.04682586274526d+04,    -1.40082238884124d+03,
     &    -7.04401579238760d+03,    -4.01603364379001d+04,
     &    -4.02156397493465d+04,    -4.44082144141766d+04,
     &    -3.01545783706270d+04,    -3.85809309528951d+04,
     &     1.77986832553758d+04,     5.96039270082918d+04,
     &     2.03769432968681d+04,     2.19485138233158d+03,
     &    -5.99855164536808d+02,    -1.00249200498664d+04,
     &    -5.90139954984246d+03,    -4.69991080613766d+03,
     &    -1.72925533887879d+04,     1.02773105307772d+04,
     &     1.92432749968431d+03,     3.06456642861916d+04,
     &     8.08625297260154d+04,     3.67025815704315d+04,
     &    -1.23516327444757d+05,    -5.32723770463839d+04,
     &    -1.27076164458941d+04,    -1.26469067357512d+04,
     &    -7.03672212191499d+03,    -1.30735937375710d+04,
     &    -2.97309218176781d+03,    -3.97285964186089d+04,
     &     9.15196437109342d+03,    -1.23321900913262d+05,
     &     1.45962069803053d+05,    -2.52291659260439d+05,
     &     1.21631451278848d+05,    -4.03867293837506d+04,
     &     7.91337181181183d+04,     1.01783361169241d+05,
     &     8.36784810485193d+04,     2.21988422581488d+05,
     &     1.00403382951704d+05,    -1.68913737646237d+05,
     &    -9.84903005351894d+04,    -3.01302459420521d+04,
     &    -2.87791233258126d+04,     6.77221367312097d+04,
     &     4.31133359765135d+04,    -3.21180660845004d+04,
     &    -5.80138100257817d+03,     3.90692532505034d+04,
     &     1.68722319216422d+04,    -3.86642488429455d+02,
     &    -2.03026580222609d+05,    -4.82455358050411d+04,
     &     6.08317716286350d+05,     2.36814060225466d+05,
     &     3.75500236344384d+04,     1.02573303579441d+04,
     &     2.86592914310625d+04,    -1.23543410324433d+05,
     &     1.98802569191338d+05,    -2.17457579849380d+05,
     &     2.63589077250618d+05,    -2.56688228869253d+05,
     &     2.45729010981872d+05,     4.52158080386007d+05,
     &    -4.98345616810165d+05,     2.49245599282706d+05,
     &    -2.92261520981632d+05,    -3.21091585711200d+05,
     &    -2.30813981807371d+05,    -9.49522618931016d+05,
     &    -6.16992004337282d+05,     7.00188613818910d+05,
     &     4.77754759327415d+05,     2.85485138907451d+05,
     &    -5.30376831745501d+05,     1.05725185090996d+06,
     &    -3.58040101856795d+04,     7.46866882388452d+04,
     &    -9.44690422954201d+04,     6.44501333715127d+04,
     &    -5.94231988279565d+04,     7.00687885831826d+04,
     &     6.40173083119075d+05,     6.72834677063777d+04,
     &    -2.56395564305750d+06,    -9.00696896244148d+05,
     &     3.30314411453406d+04,     7.61533287419704d+04,
     &     5.34783287986760d+05,    -1.69733339651386d+05,
     &     1.63842919665709d+06,    -1.54861889990687d+06,
     &     1.68169736381048d+06,     1.56473470007949d+06,
     &    -4.74090666247152d+06,     5.64873201723657d+06,
     &    -3.18320328825531d+06,     2.84847269767259d+05,
     &    -3.76150312216017d+05,     3.73957436846703d+05,
     &    -4.84690331244004d+05,     2.31825859591237d+06,
     &     1.95906800719627d+06,    -7.73604131604642d+05,
     &    -9.48099572712465d+05,    -5.44878006770972d+06,
     &     3.48594174459958d+07,    -7.59960234031616d+07,
     &     1.28416974723167d+04,    -6.62298304333324d+03,
     &    -1.11349260269912d+04,     1.87103548462631d+04,
     &    -8.51946201539171d+03,    -4.85710888246740d+03,
     &    -1.91856127390304d+05,     4.28419313619127d+04,
     &     7.23121192254250d+05,     9.72102729329752d+04,
     &    -1.36835021215310d+05,    -1.17770813842127d+05,
     &    -1.52226311036634d+05,    -9.58654317367947d+04,
     &    -7.08120571238228d+04,    -4.16430933122365d+04,
     &    -2.31964929728517d+05,     1.75956014772789d+05,
     &     4.16283870643610d+04,    -4.97929701804196d+05,
     &     4.58671654973243d+05,    -1.25810360715548d+05,
     &     7.06631271963049d+04,    -1.91885160647372d+04,
     &     2.71711593647345d+04,    -4.78483264617137d+04,
     &     9.72675948386272d+03,     5.08183311138758d+04,
     &    -2.90036134637124d+04,     2.16877093842136d+05,
     &    -3.70001946494937d+05,     2.74427996694503d+05,
     &    -2.21962143352518d+04,     5.83795187803970d+03,
     &    -1.62461583347192d+04,     5.85968825616889d+04,
     &    -1.99503595113277d+04,     9.85445612443625d+04,
     &     6.42690143252701d+05,     3.13382108886963d+03,
     &    -2.64097498866459d+06,    -2.89334428501115d+05,
     &     5.32395382366136d+05,     4.03468787318636d+05,
     &     3.33015097952563d+05,     7.96374592898144d+05,
     &    -1.02131917840975d+06,     2.15559226337454d+06,
     &    -1.23532953730095d+06,    -5.04418799212151d+05,
     &     3.45354719346847d+06,    -4.27230875612092d+06,
     &     2.52435111074508d+06,    -7.36654798782196d+05,
     &     2.15975355670695d+05,    -2.65674353468634d+05,
     &     7.12072450262751d+04,    -1.20374386300091d+05,
     &     3.25585960309774d+04,    -2.61210234501129d+05,
     &     3.68538265956578d+05,     1.37525084432332d+06,
     &    -1.53841993319502d+07,     4.22200675492475d+07,
     &     3.80681589283219d+04,    -8.06323535189497d+04,
     &     3.33403374918570d+04,     1.48280079304685d+05,
     &    -3.43939546416862d+04,    -1.10853690982455d+05,
     &    -2.08626666142765d+06,    -5.65382736635246d+05,
     &     8.77071474183005d+06,     7.56295082538464d+05,
     &    -2.15308927410658d+06,    -1.31428558976232d+06,
     &    -1.83724247815279d+06,     1.00233737391787d+06,
     &    -4.04101879993223d+06,     2.48579720912842d+06,
     &     2.25729725655038d+06,    -8.59057559074672d+06,
     &    -1.59609102312268d+06,     1.92489544643915d+07,
     &    -1.84505547294738d+07,     6.48492568235183d+06,
     &    -2.46920054032651d+06,     2.30745800857709d+06,
     &    -2.95482769002035d+06,     2.68527552805965d+06,
     &    -5.76875545279626d+05,    -1.21684644001486d+06,
     &     1.27043984628558d+07,    -2.64465586364208d+07,
     &    -1.23971078332190d+08,     5.79766798095404d+08,
     &     2.96670921005742d+04,    -3.69761365454926d+04,
     &    -1.57705427222910d+05,     2.47517578815307d+05,
     &    -4.02312839923039d+05,     3.27600515792782d+05,
     &     5.49019832837318d+06,     3.19636290201507d+06,
     &    -2.73095785887496d+07,    -2.13719825527216d+06,
     &     7.93890842693930d+06,     5.51416182404585d+06,
     &     2.86246470060342d+06,     3.67918504572831d+06,
     &     9.35188076105544d+06,    -2.76738122279730d+07,
     &     1.39614780522085d+07,     6.50147230733437d+07,
     &    -1.34569854048089d+08,     1.29837504852561d+08,
     &    -6.28303249886812d+07,     1.54943143977402d+07,
     &    -4.19106087634234d+06,    -5.45213811699838d+06,
     &     1.39975698875350d+07,    -1.59256423144573d+07,
     &     1.56085390240676d+07,    -4.24205970947047d+07,
     &     2.06599949927650d+05,     1.01968814106340d+09,
     &    -4.86269174868558d+09,     6.29070654188238d+09,
     &     1.16991366035116d+05,     2.32893770690720d+03,
     &    -4.75026558103835d+05,     7.73392218095447d+05,
     &     6.21153070513756d+04,     5.67086275510586d+05,
     &    -1.31678131852059d+07,    -1.13759207277494d+07,
     &     7.76735480418248d+07,     9.69564362469366d+06,
     &    -3.15811910415294d+07,    -1.20206219166426d+07,
     &    -2.19825562817433d+07,    -7.40094925721388d+06,
     &     2.45426636493230d+07,    -7.99865460548530d+06,
     &    -7.69168663130648d+07,     1.09612243281014d+07,
     &     4.61735520621215d+08,    -9.73774191208228d+08,
     &     7.38812097087622d+08,    -2.44723925388528d+08,
     &     9.49483780825990d+07,    -4.74758708297458d+07,
     &     4.76372135291504d+07,    -3.50068586680652d+07,
     &    -1.23229950333307d+07,     1.04758816992195d+08,
     &    -4.58254460382236d+08,     7.04145980427565d+09,
     &    -5.71310096728234d+10,     1.97359195680233d+11,
     &    -2.47848304039871d+05,     5.12643411463882d+05,
     &    -7.49384300588570d+05,     1.01082858546369d+06,
     &    -6.18655756706256d+05,    -5.25216816054747d+05,
     &     2.55069311328006d+07,     3.51876624281710d+07,
     &    -2.07919800302913d+08,    -2.01560829315055d+07,
     &     8.15775138086878d+07,     4.50788744926742d+07,
     &    -1.04423321230887d+07,     2.84205586530613d+08,
     &    -5.15475883650578d+08,     2.60854926378140d+08,
     &     1.00897032419610d+09,    -2.67934071923971d+09,
     &     3.21981132410267d+09,    -2.31818122160897d+09,
     &     9.85742627551260d+08,    -2.98949881333317d+08,
     &     6.71431553729399d+07,     5.74422030692832d+08,
     &    -1.53685833190701d+09,     1.52622069143891d+09,
     &    -9.76098293113335d+08,    -5.00812811653056d+08,
     &     2.11678734530344d+10,    -1.03110720727033d+11,
     &    -6.75651312193866d+11,     4.90746029712830d+12,
     &    -6.30839098723299d+05,     9.11338353842426d+05,
     &     5.76812584011409d+04,     7.38685197248577d+05,
     &    -1.84659524469858d+06,     9.99359791905132d+05,
     &    -3.74488555158426d+07,    -7.67805152471465d+07,
     &     4.04615670150462d+08,     7.67521702793507d+07,
     &    -2.28052919494732d+08,     3.32457632153163d+07,
     &    -2.84810961506784d+08,     6.07935113797101d+08,
     &    -2.74321563052089d+09,     6.80687412739306d+09,
     &    -1.12862860150104d+10,     1.72364874783800d+10,
     &    -2.89630715895638d+10,     3.89934637808907d+10,
     &    -2.71573109345099d+10,     9.92050728308061d+09,
     &    -4.92497892581958d+09,     7.21904255201417d+08,
     &     6.08545485048174d+09,    -8.75258324113822d+09,
     &     2.12763064276512d+09,     7.27570497132766d+10,
     &    -5.61627230017486d+11,     3.68524658446988d+11,
     &     2.78969466799225d+13,    -1.75459036695504d+14, ! id_smallr = 3, smallr = 2.05 bohr
     &     6.64662083248580d+03,    -6.32570710383072d+03,
     &     4.09626708245348d+03,    -3.40618783534509d+03,
     &     5.39829096972405d+03,     4.35924955146297d+03,
     &    -2.52409319401719d+03,    -4.43502713069963d+03,
     &    -2.27634941228706d+03,     2.11664568832531d+04,
     &     2.82259013825671d+04,     1.45614852381448d+04,
     &     1.78481957767045d+04,     9.91488204522983d+03,
     &     7.95013063085017d+03,     7.54356406200534d+03,
     &    -9.31439253032707d+03,     7.60873165678197d+03,
     &    -2.02734448014904d+04,    -6.33288124026146d+03,
     &    -2.69763349002886d+04,    -4.42571952282710d+04,
     &    -4.78449263488553d+04,    -3.13305048345460d+04,
     &    -3.43307062843437d+04,    -1.53111621427194d+04,
     &     3.36843718184483d+04,     5.47075522729851d+04,
     &     2.21106053987144d+04,    -9.23527746903472d+03,
     &     1.61641336315648d+04,    -2.55618975177552d+04,
     &    -1.31542926081654d+04,     2.13266106378908d+03,
     &    -1.17187247705316d+04,    -2.38179044773684d+04,
     &     3.45467516977259d+04,     6.64687288728088d+04,
     &     5.83951454548565d+04,     1.31827057967456d+05,
     &     1.32000030720433d+05,    -9.82711455853871d+04,
     &    -1.91388636292184d+05,    -5.78541510175138d+04,
     &    -1.33181792669640d+05,    -6.34036798799159d+04,
     &    -1.17601195833292d+05,    -9.09286680433361d+04,
     &    -1.38758145164485d+05,    -5.42995130363723d+04,
     &    -1.66035015454579d+05,     8.80480418111082d+04,
     &    -1.75274957557513d+05,     8.90549720379643d+04,
     &     2.38004938700198d+05,     2.07360008341938d+05,
     &     3.09324773644746d+05,     5.11230245916883d+05,
     &     2.92293062599357d+05,    -4.33770298949970d+05,
     &    -2.01396632144685d+05,    -2.84556151434543d+05,
     &     2.75784913451646d+05,    -1.16867857870410d+05,
     &     2.27338601394261d+04,    -2.33090438038363d+04,
     &     5.14650526422688d+03,    -1.72477678960503d+04,
     &     1.15336534243694d+04,     9.28912559077849d+03,
     &    -2.15082413828561d+04,    -3.58020209248958d+04,
     &    -3.38966188655450d+04,     5.57655511235721d+04,
     &     7.03484582621504d+04,     1.18588159796485d+04,
     &     2.65342302219879d+04,     1.20085726543505d+04,
     &     1.31591681191585d+04,     5.61057436138434d+03,
     &     1.89400336081343d+04,    -1.99539670368033d+04,
     &     3.28463240078914d+04,    -4.58330846962224d+04,
     &     4.78968578785853d+03,    -4.89036065059028d+04,
     &    -5.28998127930849d+04,    -4.10708411045887d+04,
     &    -4.74077770509471d+04,    -4.04251452969036d+04,
     &     2.22551177185981d+04,     6.84113004952470d+04,
     &     2.67845154851591d+04,    -1.11695651342892d+04,
     &     2.36292633644179d+04,    -3.01814784999526d+04,
     &    -8.62718337115751d+03,     7.80643392067432d+03,
     &    -1.97747665484209d+04,    -5.23317128255873d+04,
     &     2.04039559550408d+04,     5.95329152195388d+04,
     &     1.59047527899940d+04,     1.53035486254007d+05,
     &     1.57945915822101d+05,    -2.14782773265796d+05,
     &    -2.63690823915187d+05,     4.66312126444690d+04,
     &    -6.79046046530205d+04,     3.11130088740961d+04,
     &    -4.80139991104008d+04,     1.06114133067145d+04,
     &    -9.81923141225766d+04,     7.57473512068687d+04,
     &    -2.03049935056574d+05,     1.38636818283080d+05,
     &    -1.19480877925925d+05,     4.24947336025149d+04,
     &     7.84964468954186d+04,     1.26666624339550d+05,
     &     1.25512868281230d+05,     2.82082928501985d+05,
     &     1.27646652725034d+05,    -2.28269018301962d+05,
     &    -1.33309641531744d+05,     7.28602890847621d+04,
     &    -3.43331691354098d+05,     4.11763716311542d+05,
     &     5.07397121026375d+04,    -2.29281083428616d+04,
     &     5.46414109807556d+04,    -1.31483771064122d+05,
     &     9.06805339788299d+04,     1.87508811162010d+05,
     &    -7.34311784672670d+04,    -3.57626598391648d+05,
     &    -4.90594400280177d+05,     1.06218556152347d+06,
     &     1.04412918586809d+06,    -3.01407405624459d+05,
     &     6.82816499573532d+04,    -2.29902576622864d+05,
     &     1.28972281766274d+05,    -1.58681770802591d+05,
     &     2.27760556387806d+05,    -3.61013486175608d+04,
     &    -3.20614556667966d+04,     4.78772376086284d+05,
     &    -2.72242668864446d+05,     2.05191767780176d+05,
     &    -3.05814226234727d+05,    -7.81514246370376d+05,
     &     1.79604140335948d+05,    -1.59488001872054d+06,
     &    -7.27347804589873d+05,     6.57625564273690d+05,
     &     1.23197647310090d+06,    -5.17739977423198d+05,
     &     4.86466170059082d+06,    -1.01253670026297d+07,
     &    -8.97543151579852d+04,     1.02098872061247d+05,
     &     1.11694267003188d+05,    -2.83005900035513d+05,
     &     6.20254554429825d+04,     3.04034351084194d+05,
     &    -3.79094315015388d+05,     1.18941462808670d+06/
      DATA coe_RKHS(1001:2000)/
     &     1.38578405743177d+06,    -4.94945121915224d+06,
     &    -3.80965452794128d+06,     1.38279832808158d+06,
     &     1.34831178001408d+06,     5.73900885237471d+05,
     &     1.13790355781823d+06,     3.58499440534695d+05,
     &     1.12721395056476d+06,    -1.67080659186872d+06,
     &     5.59379282222132d+06,    -9.18415184596017d+06,
     &     7.27550440022228d+06,    -3.88904525204601d+06,
     &    -6.26851116728886d+05,     6.45930104114913d+06,
     &    -7.78397892002385d+06,     7.12941012062866d+06,
     &     3.58566054126392d+06,    -1.64828757996207d+06,
     &     2.13580314777597d+06,    -1.46971039185316d+07,
     &     9.06720679119594d+07,    -2.26279723661473d+08,
     &     3.42671673107114d+04,    -5.63578924365859d+04,
     &     6.73194882769993d+04,    -8.04008482455340d+04,
     &     3.99171833397639d+04,     8.89512902350759d+04,
     &    -8.59569275541851d+04,    -2.87247107791777d+05,
     &    -2.65418167516358d+05,     1.51020063481164d+06,
     &     6.47559042040027d+05,    -6.23381063229856d+05,
     &    -4.31996025760194d+05,    -1.62165517122520d+05,
     &    -3.63784391369448d+05,     7.49939973582905d+04,
     &    -3.65720936534850d+05,     4.46322850813470d+05,
     &    -1.15711437137396d+06,     1.79108707555930d+06,
     &    -1.31964289770063d+06,     3.80998152765824d+05,
     &     5.16969432417838d+05,    -1.27925088455366d+06,
     &     1.38228407028821d+06,    -6.95819664618859d+05,
     &     5.30577614643836d+04,     2.74401905118977d+05,
     &    -1.33038128221551d+06,     6.60180481354495d+06,
     &    -1.90138452767019d+07,     2.17336997790897d+07,
     &     8.70323314887144d+03,    -8.14653993151970d+04,
     &     2.50596426018349d+04,    -8.18213865246754d+04,
     &     1.65857502028059d+05,     2.53172434532840d+05,
     &    -1.98146979138827d+05,     1.21112590675145d+06,
     &     1.13755764099249d+06,    -6.07629980702312d+06,
     &    -2.31476719325839d+06,     2.37731522859296d+06,
     &     1.62727492853550d+06,     8.46039712001126d+05,
     &     7.64360553455307d+05,     4.92797494420023d+05,
     &     3.19753910700644d+05,    -1.07060432215893d+06,
     &     4.18604763506105d+06,    -4.74661554962902d+06,
     &     1.68885531865363d+06,     2.15645626026425d+06,
     &    -5.81084309613149d+06,     8.12502978118535d+06,
     &    -7.72021036893004d+06,     3.28969102551780d+06,
     &    -3.95671695696958d+05,    -7.77374149877851d+05,
     &     2.04418668935701d+06,     3.65693623143645d+06,
     &    -1.79377287957714d+07,     1.67173942010067d+07,
     &     9.41391816468075d+04,    -9.74001960914966d+04,
     &    -3.13947203985934d+04,    -2.47914784524819d+05,
     &     5.86965095442035d+05,     6.22302678147615d+05,
     &    -1.12214022905301d+06,    -3.64729776196481d+06,
     &    -6.13561573517122d+06,     2.25225991047256d+07,
     &     7.89014352832024d+06,    -8.48583361953188d+06,
     &    -6.36688569277422d+06,    -3.87019566820889d+06,
     &     5.23657144750294d+05,    -2.67399187744596d+06,
     &    -2.87895340858630d+06,     4.51195039582062d+06,
     &     1.10002545395922d+06,    -2.15684732583862d+07,
     &     2.60386573682412d+07,    -1.80572426669233d+07,
     &     1.76770899049777d+07,    -1.40450596108176d+07,
     &     3.00246596364408d+06,     7.44044378351138d+06,
     &    -1.33155327401404d+07,     4.23281172934760d+07,
     &    -1.44720716629440d+08,     1.05212362659095d+09,
     &    -4.91677978123471d+09,     7.26575571570076d+09,
     &    -2.38968155052875d+05,     6.63658320238169d+05,
     &    -2.70573703146772d+05,    -1.24198495002447d+06,
     &     1.03861056724297d+06,     8.01534728774312d+05,
     &    -1.67770340190780d+06,     1.02620028647873d+07,
     &     2.68662642849984d+07,    -7.60693166251142d+07,
     &    -3.03246096577224d+07,     4.15129917492965d+07,
     &     8.81636039977316d+06,     1.90661738976289d+07,
     &     5.43366092208674d+06,    -3.48493482517563d+07,
     &     9.31530339748544d+07,    -5.39826779544285d+07,
     &    -1.24658607335071d+08,     2.22787777537322d+08,
     &    -7.19434452418032d+07,    -1.22782148219670d+08,
     &     2.46207549400814d+08,    -3.50898494726531d+08,
     &     3.48681022734644d+08,    -1.69919300777617d+08,
     &    -6.69625342677651d+07,     8.79521195952594d+08,
     &    -4.92556322359158d+09,     2.08409024007907d+10,
     &    -7.93864900410233d+10,     1.74783989179979d+11,
     &    -1.89721450633417d+05,     7.80773563407222d+05,
     &     4.45555192904183d+05,    -3.15267309040585d+06,
     &     2.35046897121873d+06,     5.27134631874357d+06,
     &    -3.50247202747070d+06,    -1.57208348505621d+07,
     &    -1.02650540329749d+08,     2.37636377294842d+08,
     &     1.10383538898282d+08,    -1.74922473521223d+08,
     &    -9.20531580654369d+05,    -4.09022798125410d+07,
     &    -7.87207835056341d+07,     2.05296542249458d+08,
     &    -3.95867057075698d+08,     1.28029295627528d+08,
     &     5.17995610827469d+08,    -1.53784677562199d+08,
     &    -8.81341070419707d+08,     1.66308805229037d+09,
     &    -2.95362363296892d+09,     4.82129510415425d+09,
     &    -5.01261986461645d+09,     2.39582550952367d+09,
     &    -6.96077379523630d+08,     1.28658030719220d+09,
     &     3.37572494348541d+09,    -1.80512848617842d+11,
     &     1.36024756634634d+12,    -2.63143008457756d+12,
     &     3.44116861896802d+05,    -1.98250813107429d+06,
     &     2.65965057429537d+06,    -2.32028117389527d+06,
     &     2.63678238266184d+05,     1.15421365606916d+07,
     &    -1.98954023306505d+07,     1.78833849509024d+07,
     &     3.09254673895922d+08,    -6.51884141973217d+08,
     &    -3.44415386081899d+08,     5.54231424995403d+08,
     &    -2.47764750503129d+07,     1.84410753075414d+08,
     &    -1.83153422378839d+06,     1.75590874461472d+08,
     &    -5.27470465651660d+08,     3.26387125900698d+07,
     &     3.76632073414267d+09,    -8.28117691155975d+09,
     &     6.01122414144854d+09,    -2.67582315310332d+09,
     &     1.18905244016409d+10,    -3.35262309973909d+10,
     &     4.13465027540315d+10,    -1.96165698135227d+10,
     &    -1.55384534229574d+09,     2.31145984382949d+09,
     &    -1.10122060090003d+11,     3.09524872701535d+12,
     &    -1.80908480656646d+13,     6.48609060142481d+12,
     &     2.25768262771872d+06,    -8.21273992988308d+06,
     &     5.36622009253429d+06,     5.87607320179183d+06,
     &    -7.45293230155509d+06,     7.90362947963592d+06,
     &    -4.06291211250502d+06,     1.98942628528257d+07,
     &    -6.61436176499867d+08,     1.29452246119808d+09,
     &     9.47190526908077d+08,    -1.63007056792381d+09,
     &     9.24849622436427d+08,    -1.81715068594841d+09,
     &     1.72096171884883d+09,    -3.11509557818521d+09,
     &     6.09286308264579d+09,     3.80361987785378d+09,
     &    -3.97432904513647d+10,     5.27088537124489d+10,
     &    -8.35595277813118d+09,    -3.09733052589303d+10,
     &    -1.38370127417078d+10,     1.44911301022777d+11,
     &    -2.11180862676060d+11,     1.20017104872318d+11,
     &    -1.30753203342447d+11,     7.02525780284299d+11,
     &     2.27000677957914d+12,    -6.46813267448248d+13,
     &     3.03318775784464d+14,     2.27160578034688d+12, ! id_smallr = 4, smallr = 2.1322 bohr
     &     9.30291290023762d+03,    -1.68707934107300d+04,
     &     2.37417431850213d+04,    -2.29840700312851d+04,
     &     1.08473933467024d+04,     6.97008505320651d+03,
     &    -9.88098996165775d+02,    -5.62547862122911d+03,
     &    -7.54706070056746d+03,     7.54511039810655d+03,
     &     4.17894296949179d+04,     2.10500620013404d+04,
     &     1.76491438038639d+04,     1.30996877745128d+04,
     &     1.01575557291427d+04,     1.19671836864320d+03,
     &     3.86776559314529d+03,    -9.71626590177418d+03,
     &    -8.41369058355553d+03,    -5.38493524296729d+03,
     &    -3.44810479932883d+04,    -4.37304918146010d+04,
     &    -4.76346627437031d+04,    -3.81756174686813d+04,
     &    -3.18045955624292d+04,    -1.76003222579292d+04,
     &     3.39646897698536d+04,     5.77585584159300d+04,
     &     1.66102921127371d+04,     6.61978681682646d+03,
     &    -1.71318136688128d+04,     9.54211300587337d+03,
     &    -1.02609271379744d+04,    -1.31153899444754d+04,
     &     1.82783502928212d+04,    -5.23461891896562d+04,
     &    -1.39695699096334d+03,     1.16036194079585d+05,
     &     5.40412348194845d+04,     1.28043149742137d+05,
     &     1.86537434917144d+05,     5.72706712842228d+04,
     &    -3.29160006163420d+05,    -1.27216763738072d+05,
     &    -1.22646932982334d+05,    -1.06115579089227d+05,
     &    -1.25308218728203d+05,    -9.20292478735126d+04,
     &    -1.06383397249265d+05,    -1.25434133696157d+05,
     &    -8.38132508831122d+04,    -4.73356655634123d+04,
     &    -5.36794978002169d+04,     5.29481205700296d+04,
     &     2.26770485526983d+05,     3.12049673881423d+05,
     &     2.15467415424744d+05,     5.77608241100703d+05,
     &     2.95901301822163d+05,    -4.53382441957139d+05,
     &    -2.35933505331527d+05,    -2.66639362945767d+05,
     &     2.86047032782104d+05,    -6.27838197118589d+04,
     &     2.03594166654363d+04,    -2.88348136283803d+03,
     &    -3.54232159434473d+04,     2.95821723792001d+04,
     &    -4.20803719149624d+04,     5.10097849876044d+04,
     &    -2.67566591169507d+04,    -3.11021373652771d+04,
     &    -5.70896721951448d+04,    -6.92589096324756d+02,
     &     1.31410830901166d+05,     2.90631439933774d+04,
     &     2.52135924658174d+04,     1.70378748424333d+04,
     &     1.43257325228574d+04,     9.57754340834935d+03,
     &     4.41060595405581d+02,     9.58953187326263d+03,
     &    -8.84611731574196d+03,    -9.26570310525177d+02,
     &    -2.57069651673345d+04,    -4.35742104809553d+04,
     &    -5.75021992764376d+04,    -5.29418933814586d+04,
     &    -4.00898714078886d+04,    -4.71415428946747d+04,
     &     2.35665357087974d+04,     7.55798961143935d+04,
     &     1.99371709344841d+04,     1.25670629400131d+04,
     &    -1.53501452838937d+04,    -4.85719329185508d+03,
     &    -2.11798335130932d+04,     6.09053777744420d+04,
     &    -1.18979809723524d+05,     7.62267571704875d+04,
     &    -1.53673966048491d+05,     1.83774576521849d+05,
     &    -2.27882955632471d+04,     1.22190200240889d+05,
     &     2.43179166240442d+05,     5.26740792076120d+04,
     &    -5.55851005765813d+05,    -2.39472677019525d+04,
     &    -2.65065453056867d+04,    -8.88199581881798d+03,
     &    -2.21058561284014d+04,    -2.22530229584819d+04,
     &    -4.31008315691789d+04,     3.66671162131884d+04,
     &    -1.85843022963521d+05,     1.14943678663031d+05,
     &    -1.06201167925748d+05,     4.25740116618361d+04,
     &     8.95990184767700d+04,     1.55548003422931d+05,
     &     1.27891179658177d+05,     3.13685480861916d+05,
     &     1.37913136532929d+05,    -2.41937205475236d+05,
     &    -1.39663575680448d+05,    -9.55669111935759d+04,
     &     4.57441241288375d+05,    -7.43788917687190d+05,
     &     7.77690225140225d+04,    -1.30364445159137d+05,
     &     2.77950760960584d+05,    -2.66489845337639d+05,
     &    -9.53305775658538d+04,     4.64510972993346d+05,
     &    -9.05461749264735d+04,    -1.69029256883071d+05,
     &    -8.14920458115779d+05,    -2.82840115491422d+04,
     &     2.43017017551793d+06,    -1.75715459171100d+05,
     &    -1.03135573901269d+05,    -1.52974774750526d+05,
     &    -1.74635502365039d+04,    -1.00350405154954d+05,
     &     1.38532649775909d+05,    -1.23591321387724d+05,
     &     5.42803355292536d+05,    -6.06247076213001d+05,
     &     6.68396888056152d+05,    -3.14658497578969d+05,
     &    -1.18015713476827d+05,    -6.06955267151825d+05,
     &    -4.23982029348827d+05,    -1.44216192141375d+06,
     &    -9.95506172736682d+05,     1.01452961695206d+06,
     &     7.82352981537541d+05,     6.02427536277511d+05,
     &    -3.14388748842122d+06,     6.85055373622612d+06,
     &    -2.15309931933656d+04,    -3.49707550108509d+05,
     &     1.09242299692768d+06,    -1.07358021673676d+06,
     &    -9.55954950209211d+04,     9.44139954870627d+05,
     &    -6.89409041383544d+05,     4.12726490771434d+05,
     &     2.39162094655096d+06,    -7.51334071754981d+05,
     &    -1.03580830810029d+07,     1.83900129463001d+06,
     &     1.22935579409801d+06,     1.84766786456173d+06,
     &     4.43816437162523d+05,     2.08231859438139d+06,
     &    -7.81233100955325d+05,     6.20606870438124d+05,
     &     1.15148431285078d+06,    -1.09950902028930d+06,
     &    -2.34066095790628d+05,    -2.31200798468006d+05,
     &    -1.79832731713498d+05,    -6.98792072972715d+05,
     &     1.70381571710989d+06,     2.88946986040786d+06,
     &     4.48682638412953d+06,    -8.03432971536915d+04,
     &    -1.37216575864775d+07,     1.01136579746543d+08,
     &    -3.98300759274408d+08,     5.69044840419978d+08,
     &     4.13472649724447d+04,    -4.78578636635290d+04,
     &    -5.56667161375484d+03,     1.07668414427839d+05,
     &    -2.84717703593399d+05,     4.28914785975915d+05,
     &    -2.28775137985450d+05,    -9.07618876530943d+04,
     &    -5.94300018255979d+05,     5.83988803249101d+05,
     &     2.45613072176292d+06,    -9.20293220287996d+05,
     &    -4.71677701540620d+05,    -6.23680397977306d+05,
     &    -3.33742657623344d+04,    -6.22109843853254d+05,
     &     6.42067358628366d+05,    -1.00886148708751d+06,
     &     5.65521842458561d+05,     2.23642009858758d+05,
     &    -4.96869756376775d+05,     3.78337815227441d+05,
     &    -1.89601664964814d+05,     1.62525380311480d+05,
     &    -1.07225641681148d+04,    -1.28951510363557d+05,
     &     1.28745153855679d+04,    -3.20529514718610d+04,
     &     4.80308446149429d+05,    -4.01793244535609d+05,
     &    -1.37991832168804d+06,     3.59378808287094d+06,
     &    -9.48618751135024d+04,     4.88314293862175d+05,
     &    -1.34167586097214d+06,     1.65927679701004d+06,
     &    -1.47656380465543d+06,     1.65342583870197d+06,
     &    -9.28184752857631d+05,     8.20826038941556d+05,
     &     2.18113400380602d+06,    -2.45161466410185d+06,
     &    -9.86762468645535d+06,     3.83491468359703d+06,
     &     1.99052249194472d+06,     2.36181239809365d+06,
     &     7.49106234277107d+04,     2.40680620024362d+06,
     &    -2.35061199967506d+06,     2.86142375186811d+06,
     &     7.12931020789355d+05,    -4.51124460532876d+06,
     &     4.48855534701402d+06,    -2.76622217255638d+06,
     &     1.49430124267922d+06,    -5.79984478191557d+05,
     &    -1.49545509923683d+06,     1.60617543903276d+06,
     &    -8.16333425465225d+05,     6.89350450206680d+05,
     &     3.27701873918806d+06,    -6.01852138775644d+07,
     &     2.77199082106380d+08,    -4.39804165019974d+08,
     &     2.31782550630047d+04,     4.29666059276274d+05,
     &    -1.34061031774737d+06,     1.62929140252772d+06,
     &    -2.07555426837561d+06,     3.82510905334715d+06,
     &    -2.62009374667171d+06,    -1.55929421602120d+06,
     &    -1.04355908843156d+07,     9.20738034348208d+06,
     &     3.63843961463283d+07,    -1.37945408328676d+07,
     &    -8.89997558118545d+06,    -7.34188039476992d+06,
     &    -2.81477446217466d+04,    -7.60344686702689d+06,
     &     6.93597134343321d+06,    -1.27737500226754d+07,
     &     1.70477944405578d+07,    -2.08485933956844d+07,
     &     1.17280219818938d+07,    -7.03149693100445d+05,
     &     4.01181988750192d+06,    -1.38854136725153d+07,
     &     2.39260859154215d+07,    -1.60012258839184d+07,
     &     1.92017005521169d+06,     3.39634263285099d+07,
     &    -1.69769877393401d+08,     6.77257965323039d+08,
     &    -2.12625892451557d+09,     3.10128992235120d+09,
     &     1.32641310747878d+05,    -1.98676811833206d+06,
     &     6.84616983511904d+06,    -9.24120811866869d+06,
     &     2.66201648087142d+06,     4.86584749402878d+06,
     &    -5.01018351476964d+06,     3.96490388387544d+06,
     &     3.82571604831476d+07,    -2.35519700466045d+07,
     &    -1.32758985729969d+08,     5.19755441758269d+07,
     &     3.51692435218462d+07,     2.11509740390095d+07,
     &    -2.31773912620034d+05,     1.97662698441947d+07,
     &    -1.99216575522884d+06,    -3.73775073156785d+06,
     &    -3.71818729134628d+07,     1.36488848861317d+08,
     &    -1.44615707769147d+08,     1.05371439122428d+08,
     &    -1.43260854441343d+08,     1.74459056696873d+08,
     &    -1.45924655070200d+08,     5.57692869994057d+07,
     &     4.98121411646122d+07,    -4.33201060898652d+08,
     &     1.72711993758630d+09,    -5.40259934214818d+09,
     &     3.11441896717059d+10,    -9.97247989532397d+10,
     &     7.80431652999817d+05,    -5.48790873670375d+06,
     &     1.72445420313564d+07,    -2.21522255623060d+07,
     &     6.50991885034209d+06,     1.42766257420568d+07,
     &    -9.61739258744781d+06,     7.51676330616430d+06,
     &    -1.28589563402776d+08,     4.93870734656328d+07,
     &     4.43737393778269d+08,    -1.75248608481338d+08,
     &    -1.26018377496863d+08,    -7.18213202453148d+07,
     &     1.03090415307197d+07,    -2.28903046191373d+07,
     &    -1.34433889487662d+08,     2.61089545473369d+08,
     &    -1.64852754337466d+08,    -1.21202621239336d+08,
     &     3.87060253014961d+08,    -9.01243103996030d+08,
     &     1.69709603070000d+09,    -1.80807538767104d+09,
     &     9.53901082746309d+08,    -1.84018439028282d+08,
     &     5.35328082183621d+08,    -5.24646997855084d+09,
     &     5.50866250222750d+10,    -5.87336534755395d+11,
     &     3.42409631271169d+12,    -7.13690887269403d+12,
     &    -9.15153644254639d+05,     7.02076470692982d+06,
     &    -2.50940699206119d+07,     4.08182946459471d+07,
     &    -4.51281625420640d+07,     5.31468948430189d+07,
     &    -3.30876867104832d+07,    -5.15228957971898d+07,
     &     3.43330024912763d+08,    -5.68822366289855d+07,
     &    -1.31165914411000d+09,     5.04394152431357d+08,
     &     3.93038444942547d+08,     2.19564232262092d+08,
     &     6.76884176160782d+07,    -4.91577148066683d+08,
     &     1.99424642584139d+09,    -3.66960602474209d+09,
     &     4.26372863256495d+09,    -3.97572779416603d+09,
     &     1.57503668723481d+09,     3.91956153055589d+09,
     &    -1.03267527382572d+10,     1.09002549727691d+10,
     &    -3.57757243482766d+09,    -5.57018175316439d+08,
     &    -1.52459238921908d+10,     1.81112007742583d+11,
     &    -1.65567743269742d+12,     1.08163868742774d+13,
     &    -3.22129695445357d+13,    -2.03302153308514d+12,
     &    -3.74639124355539d+06,     3.83050100222511d+07,
     &    -1.39649128004096d+08,     2.15147162409359d+08,
     &    -1.73855579868416d+08,     1.17028961055323d+08,
     &    -9.70292648275673d+07,     1.91055716345504d+08,
     &    -6.63737074260371d+08,    -7.96889100155252d+07,
     &     3.01191756637674d+09,    -1.18421482600117d+09,
     &    -7.56153137689906d+08,    -5.88397835819699d+08,
     &    -5.15728862146983d+08,     3.18500780937457d+09,
     &    -1.14246019380329d+10,     2.40202493632043d+10,
     &    -3.33933458658759d+10,     3.31267396967859d+10,
     &    -1.71544934312941d+10,    -9.17714443448830d+09,
     &     3.36579509548115d+10,    -3.27778304465617d+10,
     &    -1.48144261702267d+09,     1.91778741858747d+10,
     &     4.29565254435815d+10,    -9.11658157068787d+11,
     &     1.06713778412268d+13,    -9.24331008451165d+13,
     &     3.93782980634984d+14,    -4.35165744905883d+14, ! id_smallr = 5, smallr = 2.20 bohr
     &     9.17242467918183d+03,    -1.26411945987086d+04,
     &     1.08463489448761d+04,    -6.27018447813184d+03,
     &    -2.17040622419556d+03,     8.47532578320291d+03,
     &     6.29423648867025d+03,    -8.85086468341266d+03,
     &    -8.02186713204110d+03,    -5.92222942777384d+03,
     &     4.14724681187346d+04,     3.69910845385652d+04,
     &     1.78567246761802d+04,     1.48519703646931d+04,
     &     9.91325158454362d+03,     5.41107992463380d+03,
     &    -2.83874395231217d+03,    -2.37617159914392d+03,
     &    -1.00918748356505d+04,    -1.52133235480273d+04,
     &    -2.48418726703417d+04,    -4.95330510930411d+04,
     &    -4.60164278255304d+04,    -4.35619324020375d+04,
     &    -2.80366553827516d+04,    -2.19787164639954d+04,
     &     3.54643611281099d+04,     5.90613901226298d+04,
     &     1.85456604861770d+04,    -1.32175187233111d+03,
     &     1.42801329035675d+04,    -3.10525713294892d+04,
     &    -1.15323239015661d+04,    -4.66632676867119d+03,
     &    -5.82054433188367d+03,    -1.85559882085074d+04,
     &    -4.52253536486745d+04,     9.80389670842606d+04,
     &     1.22663800765429d+05,     9.75045875212051d+04,
     &     2.03355533723879d+05,     2.01879777581620d+05,
     &    -3.11984358190269d+05,    -2.97604917997135d+05,
     &    -1.24234662327585d+05,    -1.21857327157151d+05,
     &    -1.25189134733690d+05,    -1.17841105155783d+05,
     &    -8.93589408878967d+04,    -1.41881243760288d+05,
     &    -1.16580107797339d+04,    -1.86313690089269d+05,
     &     5.33510690271718d+04,     3.37094798563802d+04,
     &     2.31034592837973d+05,     3.49326526878919d+05,
     &     2.07297108774440d+05,     6.12179360716936d+05,
     &     2.99529126869479d+05,    -4.25618475271928d+05,
     &    -2.87759426027507d+05,    -2.41892104797751d+05,
     &     1.10015199786832d+06,    -1.83943827052499d+06,
     &     2.32769244513184d+04,    -1.30230875667336d+04,
     &    -7.99170073799028d+03,    -2.35619242202981d+02,
     &    -3.31487642314422d+04,     3.20572672059548d+04,
     &     1.26396544419808d+04,    -4.32084483155370d+04,
     &    -5.84330691865856d+04,    -5.40802752089900d+04,
     &     1.29658387674769d+05,     9.40540091699616d+04,
     &     2.00534444847565d+04,     1.73907062019566d+04,
     &     1.74134485190659d+04,     6.77040769735249d+03,
     &     7.22274105143554d+03,    -1.57142713677743d+03,
     &     1.17607154352345d+03,    -1.74357243144735d+04,
     &    -1.48583801882590d+04,    -5.21778551387908d+04,
     &    -5.62374506988947d+04,    -6.20761331573555d+04,
     &    -3.66954536012030d+04,    -5.27020103605329d+04,
     &     2.61019875827574d+04,     7.78374217511275d+04,
     &     2.61184620215446d+04,    -5.96074154850863d+03,
     &     3.23015817221496d+04,    -5.59977277354638d+04,
     &    -2.00724236001847d+04,     3.77410487636665d+04,
     &    -4.62950420237300d+04,     2.11639431944099d+03,
     &    -1.34818535784739d+05,     9.86057006204442d+04,
     &     1.10550102222965d+05,     2.33331901617839d+04,
     &     2.74640303252947d+05,     2.73189183632919d+05,
     &    -5.30991850231234d+05,    -3.41737173160796d+05,
     &     2.60707856165125d+04,     8.69507064859980d+02,
     &    -3.90605472925713d+04,    -2.84933016146093d+03,
     &    -5.61458415742818d+04,     1.29978931258344d+04,
     &    -1.89646924946526d+05,     2.53902505012862d+05,
     &    -2.59275788573175d+05,     1.04646776942282d+05,
     &     7.26343749612872d+04,     2.00412545038257d+05,
     &     9.80608473298074d+04,     3.61264177792625d+05,
     &     1.44717756495353d+05,    -2.64205787901140d+05,
     &    -1.53479002984930d+05,    -8.64410108268964d+04,
     &     3.37557011862711d+05,    -5.01664450850111d+05,
     &     8.31764905235303d+04,    -1.06633367759903d+05,
     &     1.68279959898972d+05,    -3.51793042952898d+04,
     &    -2.82844864100825d+05,     2.97028788719765d+05,
     &     3.34763812038309d+05,    -2.29173486136594d+05,
     &    -7.55385630882274d+05,    -9.52551048980182d+05,
     &     2.34561553659933d+06,     1.27002277791487d+06,
     &    -4.92637703526506d+05,    -2.45309996227908d+05,
     &    -4.28391295217129d+04,    -1.39436467443997d+04,
     &    -1.30721577164775d+05,     1.71041210096792d+05,
     &     5.21691878491443d+05,    -9.16017584966114d+05,
     &     8.79478284688830d+05,    -3.10886403940294d+05,
     &    -1.08164511839617d+05,    -9.62659720682966d+05,
     &    -3.50268824929575d+04,    -1.84086590363379d+06,
     &    -1.17675325071525d+06,     1.24316984959677d+06,
     &     6.11844133526727d+05,     2.52064303877337d+06,
     &    -1.55729668629771d+07,     2.96127254015804d+07,
     &    -6.02336155987831d+04,    -1.36036115552655d+05,
     &     3.68248207798201d+05,     5.29932363524687d+04,
     &    -7.82363184462642d+05,     4.29163056788472d+05,
     &     5.02560612653864d+05,    -1.00833451650623d+06,
     &     2.36516657146105d+06,     2.43351243732399d+06,
     &    -1.06193845700710d+07,    -4.63279905098213d+06,
     &     3.48692370038165d+06,     1.53366318222775d+06,
     &     1.70666308687465d+06,    -7.52687810057620d+05,
     &     3.14195922644517d+06,    -3.16596101869206d+06,
     &     1.78066671696522d+06,    -2.19083437955835d+06,
     &     4.03773802779365d+05,    -1.31706093958787d+06,
     &    -3.75250971043689d+06,     3.94877561860552d+06,
     &    -3.87602952221606d+06,     3.95905134281776d+06,
     &     3.59917372342901d+06,    -8.65119651140544d+06,
     &     5.75973473594887d+06,     1.24028852134717d+06,
     &    -3.56427909389637d+08,     1.00826746012387d+09,
     &     4.92811008606204d+04,    -4.92607670220989d+04,
     &    -3.19088019963024d+04,     1.48557036543581d+05,
     &    -2.48171773412126d+05,     1.75019129012662d+05,
     &     1.77938713243337d+05,    -2.42933184722636d+05,
     &    -4.82815247205190d+05,    -3.75073409527164d+05,
     &     2.93881719879948d+06,     4.00704021392660d+05,
     &    -1.30175131500927d+06,    -4.15066901234994d+05,
     &    -6.98951473941550d+05,     1.61456070234374d+05,
     &    -6.30529670363409d+05,     4.64410940920203d+05,
     &    -8.42523051928481d+05,     1.46602483201702d+06,
     &    -1.19755204912054d+06,     3.86782303457663d+05,
     &     3.11492748628126d+05,    -4.75728866776577d+05,
     &     2.44724323035028d+05,    -6.24690303618450d+04,
     &    -1.40897004809141d+05,     3.88813058264749d+05,
     &    -1.02950703322120d+06,     3.35609976699166d+06,
     &    -8.44655950762188d+06,     1.36947021229480d+07,
     &    -7.57077447668953d+04,     2.76313913006469d+05,
     &    -6.41815323082184d+05,     5.87255233546438d+05,
     &    -5.68646233885007d+05,     6.70542738055656d+05,
     &     4.64999829672516d+05,    -4.50866600565859d+05,
     &     2.18701935989035d+06,     1.48490134855013d+06,
     &    -1.27251584743216d+07,    -1.00857872874597d+06,
     &     5.51234320946460d+06,     1.36480642244738d+06,
     &     2.92266812375583d+06,     1.30530289798387d+05,
     &     1.74658089490664d+05,     4.22838059951528d+05,
     &     5.76269112753749d+06,    -1.37586266971833d+07,
     &     1.12826052615897d+07,    -3.46847244302267d+06,
     &    -4.74206218974850d+05,     9.60039939654587d+05,
     &    -1.18208387060691d+06,     6.87911849662210d+05,
     &    -1.02738585144536d+05,    -9.00939772994119d+05,
     &    -5.08987125165571d+05,     3.54726974995734d+07,
     &    -2.24103555016444d+08,     4.82966199567799d+08,
     &     8.04142396481334d+04,     1.20114374002230d+05,
     &    -3.53497032994515d+05,     1.21563981391822d+05,
     &    -7.22907548476936d+05,     1.57696253416502d+06,
     &     1.51642902701512d+06,    -4.05551205243260d+06,
     &    -6.81754218385143d+06,    -1.04521210544041d+07,
     &     5.37824530261914d+07,    -2.47575453287578d+05,
     &    -1.97009630236119d+07,    -6.94078613650425d+06,
     &    -4.75792517258865d+06,    -1.18598435488540d+07,
     &     2.05116728221667d+07,    -2.54167177819975d+07,
     &    -1.02729160551946d+07,     5.71928744950402d+07,
     &    -4.95067873005274d+07,     1.72505788122825d+07,
     &    -1.13126316945729d+07,     1.59129718012006d+07,
     &    -7.33986222140573d+06,    -2.51552471153505d+06,
     &     8.24090692348324d+06,    -2.80423828529786d+07,
     &     2.42975978382250d+08,    -1.79970911435223d+09,
     &     8.12518814908539d+09,    -1.56896786768236d+10,
     &    -3.10573688827010d+04,    -9.14359249339569d+05,
     &     3.37891555606293d+06,    -3.65919450988108d+06,
     &    -9.88582636489250d+05,     1.54394072471276d+06,
     &     4.04508828481583d+06,    -4.97636164038288d+06,
     &     2.68082639144597d+07,     5.31781181221126d+07,
     &    -2.03860487883405d+08,     1.01170621743523d+07,
     &     6.76880899370552d+07,     3.83619165816167d+07,
     &    -1.79792806985483d+07,     8.90934010850842d+07,
     &    -1.47998521868465d+08,     1.93947325215299d+08,
     &    -9.98625589905230d+07,    -4.48866035974501d+07,
     &     6.28665660676505d+07,    -4.18099592449293d+07,
     &     1.17046542145090d+08,    -1.90713356742753d+08,
     &     1.10451675242705d+08,     9.18879435323640d+06,
     &    -1.15603594286114d+08,     4.89537340401369d+08,
     &    -2.01473854403526d+09,     3.67241279294392d+09,
     &     1.49103155305112d+10,    -6.96716508328596d+10,
     &     5.19913471401159d+05,    -2.76273078738843d+06,
     &     7.35678519145546d+06,    -5.58751811544328d+06,
     &    -4.57020861213818d+06,     5.41429981158780d+06,
     &     1.63311121145237d+07,    -9.71876567359197d+06,
     &    -5.78422902322119d+07,    -2.39579671856842d+08,
     &     7.12001974206359d+08,    -6.12220958328496d+07,
     &    -2.24819604011597d+08,    -1.27160561853599d+08,
     &     5.45919801738316d+07,    -2.52695912047123d+08/
      DATA coe_RKHS(2001:3000)/
     &     4.95850244444153d+08,    -9.66574804809376d+08,
     &     1.55139925441219d+09,    -1.93420053533636d+09,
     &     1.28441100324358d+09,    -5.47108174606552d+07,
     &    -1.13310016930181d+09,     2.20998150797873d+09,
     &    -2.20907793365365d+09,     1.03152465000363d+09,
     &    -6.26470660619169d+05,    -3.46694794066358d+09,
     &     2.88793422360804d+10,     3.31276163496615d+10,
     &    -1.47669030512660d+12,     5.46451627821445d+12,
     &    -5.93278731194634d+05,     3.77576808659335d+06,
     &    -1.38821016970298d+07,     2.24123543586404d+07,
     &    -2.21037462578316d+07,     9.67604694737298d+06,
     &     4.10180037502739d+07,    -1.11744161266905d+08,
     &     1.39445418807501d+08,     8.16304621891909d+08,
     &    -2.13306364528401d+09,     1.92423628509697d+08,
     &     7.69941018812756d+08,     1.01246484621224d+08,
     &     3.90490661013374d+08,    -2.05776723929008d+08,
     &    -3.56065557822399d+08,     3.20105778689202d+09,
     &    -9.17140978171102d+09,     1.29974373017626d+10,
     &    -6.91266890968739d+09,    -3.65848832494708d+09,
     &     1.32361009901736d+10,    -2.27499870681827d+10,
     &     2.51478031253842d+10,    -1.37035023892805d+10,
     &     4.29672637450223d+09,    -6.97605966734268d+10,
     &     1.10664711746714d+12,    -6.86502002945672d+12,
     &     6.72296642018898d+12,     5.40874613909347d+13,
     &    -1.52528855074940d+06,     1.84698063726755d+07,
     &    -6.56230075562060d+07,     8.63615989776936d+07,
     &    -4.47196626436713d+07,     3.65324542597928d+05,
     &     1.26242770223057d+07,     1.10771213080232d+08,
     &    -1.54023282764504d+08,    -1.88598260453452d+09,
     &     4.56068220874631d+09,    -5.29528863966492d+07,
     &    -2.56437172335206d+09,     1.73994175446398d+09,
     &    -4.18435340872270d+09,     4.92095274060275d+09,
     &    -4.03287285303362d+09,    -5.93038612160455d+09,
     &     2.89916071416690d+10,    -4.24495306205053d+10,
     &     1.76977389513044d+10,     2.60669356907361d+10,
     &    -7.28470128354646d+10,     1.21770226639512d+11,
     &    -1.33653525530421d+11,     7.65014523742856d+10,
     &    -7.63549066679049d+10,     6.47365381861787d+11,
     &    -8.86856416737981d+12,     1.32048154085927d+14,
     &    -1.22222865732119d+15,     4.72494669093624d+15, ! id_smallr = 6, smallr = 2.25 bohr
     &     9.51487251006047d+03,    -1.21054197466637d+04,
     &     7.85115727021079d+03,    -2.26240777138109d+03,
     &    -5.29670701290820d+03,     6.81671877395384d+03,
     &     7.21582737299122d+03,    -2.47751252210871d+03,
     &    -1.34603356839647d+04,    -8.28524606120770d+03,
     &     2.44745970614578d+04,     6.04431286559765d+04,
     &     1.61494107610818d+04,     1.51902371396955d+04,
     &     1.45523600755106d+04,    -3.79009278470092d+02,
     &     2.39695604878117d+03,     6.50986349282995d+01,
     &    -2.49784720471953d+04,     4.77600022415576d+03,
     &    -3.73960780038557d+04,    -4.76131401525951d+04,
     &    -5.08907152249500d+04,    -3.79081621864361d+04,
     &    -3.62101073833803d+04,    -1.92052947869679d+04,
     &     3.47127273570535d+04,     6.14724484221496d+04,
     &     1.48475023691002d+04,     1.13123728100032d+04,
     &    -1.76323195950401d+04,     4.00952536676205d+03,
     &    -1.18897051316821d+04,    -2.90910084294972d+03,
     &    -1.15135402283886d+04,    -1.12105790601698d+04,
     &    -5.71475014491420d+04,     7.55831135524532d+04,
     &     1.24700496039481d+05,     1.56361045452097d+05,
     &     1.65261157817526d+05,     2.52843364916456d+05,
     &    -1.11833975882921d+05,    -5.71526856113854d+05,
     &    -8.81914734194236d+04,    -1.42343598426198d+05,
     &    -1.61528633269771d+05,    -6.84373210778757d+04,
     &    -1.51100681271519d+05,    -1.39161447457047d+05,
     &     3.63198220651127d+04,    -1.92818208311984d+05,
     &     1.41989437000984d+04,     4.69189852395014d+04,
     &     2.86558345326123d+05,     2.44890704692692d+05,
     &     3.52280245191350d+05,     5.40936734293978d+05,
     &     3.25750892417143d+05,    -4.85462540868917d+05,
     &    -1.95650661660700d+05,    -3.55257153771430d+05,
     &     4.08414439005059d+05,    -2.53052860951134d+05,
     &     2.45495605611316d+04,    -1.52904130543980d+04,
     &    -6.06030904054558d+02,    -7.10104873369793d+03,
     &    -3.13061470997067d+04,     1.94009008363259d+04,
     &     1.67692306515882d+04,    -1.33557164500349d+04,
     &    -8.01631756813042d+04,    -6.00583053267660d+04,
     &     5.29337648407077d+04,     1.99666134268757d+05,
     &    -2.11668718754166d+03,     3.01558845594419d+04,
     &     1.17716520814574d+04,     8.21305498083987d+03,
     &     1.16139635694869d+04,    -7.20200310657857d+03,
     &    -3.28713634873401d+03,    -1.27902035763317d+04,
     &    -1.72885943131404d+04,    -5.44798547240884d+04,
     &    -6.43244487064844d+04,    -5.11386198049583d+04,
     &    -5.64391137735160d+04,    -4.41391648486715d+04,
     &     2.43940233932790d+04,     8.15906985318307d+04,
     &     2.64541810659243d+04,     2.19807495330759d+02,
     &    -5.88858002588283d+03,    -1.82904649062682d+03,
     &    -2.14496338035181d+04,     3.40798405101387d+04,
     &    -2.91104867766529d+04,    -8.66063910258311d+03,
     &    -1.29546733456538d+05,     4.37948505221398d+04,
     &     1.00944362118525d+05,     1.21492936490633d+05,
     &     1.86483666286669d+05,     3.37295957832236d+05,
     &    -1.57144812998483d+05,    -8.72196326845406d+05,
     &     1.69002605674250d+05,    -7.14007426330951d+04,
     &     2.92502225955160d+04,    -5.96653616883607d+04,
     &    -5.29203154537493d+04,     7.71950490919820d+04,
     &    -2.25869161636177d+05,     1.92250259882482d+05,
     &    -1.72776550710166d+05,     5.71002599493463d+04,
     &     1.32558412918697d+05,     1.18635776889110d+05,
     &     2.29330750494560d+05,     3.11072997146092d+05,
     &     1.63984004656816d+05,    -2.70907464384555d+05,
     &    -2.23744374623526d+05,     1.36614761380336d+05,
     &    -7.37621584383071d+04,    -1.85493168332088d+05,
     &     9.18340582702371d+04,    -1.11509852637364d+05,
     &     1.40502994465503d+05,     4.82630292640022d+04,
     &    -2.88340642614726d+05,     1.56567283230591d+05,
     &     2.51063735559088d+05,     3.33254377277028d+05,
     &    -1.06987548203993d+06,    -9.04540955631576d+05,
     &     6.55764259457518d+05,     3.68510708275424d+06,
     &    -1.23707958632585d+06,    -5.84285216309139d+02,
     &    -3.14837741223399d+05,     1.55314598378097d+05,
     &    -8.04730031824415d+04,    -2.91903766150806d+05,
     &     1.00908572569264d+06,    -8.15621269062220d+05,
     &     5.16535136197686d+05,    -5.49285701019539d+04,
     &    -5.12199659357073d+05,    -1.03775521629451d+05,
     &    -1.18742084049471d+06,    -1.39215509304871d+06,
     &    -1.23786599295423d+06,     1.22437742519703d+06,
     &     1.45346039129212d+06,    -6.10113659479752d+05,
     &    -1.58030052391526d+06,     5.46397652959530d+06,
     &    -6.76743562030013d+04,    -9.53339181361076d+04,
     &     1.82781114618066d+05,     3.32445627351649d+05,
     &    -6.40962269157973d+05,     5.79049296746513d+04,
     &     1.01430401679697d+05,     3.00772057555439d+05,
     &     1.05929881404490d+06,     2.98438399356197d+06,
     &    -4.18633307027281d+06,    -1.46839092992710d+07,
     &     7.19021299589279d+06,     1.28704232136556d+06,
     &     2.48756422162428d+06,     2.82106924498215d+05,
     &     3.14667526286123d+06,    -2.62907987516256d+06,
     &     3.01003982033120d+06,    -4.41661065940144d+06,
     &     3.37427210718046d+06,    -2.19639638354712d+06,
     &     8.18658179501499d+05,    -3.28296316232488d+06,
     &     5.31377606258001d+06,     1.98564941657239d+06,
     &     6.67123930788037d+06,    -5.45803001766504d+06,
     &     1.02229395419903d+07,    -4.44731313546840d+07,
     &     1.16347706674995d+08,    -1.80677813350724d+08,
     &     5.54599592028887d+04,    -5.08596766822129d+04,
     &    -3.53331328870719d+04,     1.08660337482903d+05,
     &    -1.07163025761410d+05,     3.57598658943134d+04,
     &     7.50374959809774d+03,     3.11404933818762d+05,
     &    -8.56719848577551d+05,    -2.91984744328014d+05,
     &     1.58448522982200d+06,     2.86504294258863d+06,
     &    -2.24426805377779d+06,    -4.60382261334266d+05,
     &    -6.70264601477765d+05,    -1.57614733547481d+05,
     &    -7.34062070458352d+05,     1.28693045702786d+06,
     &    -2.39435789286659d+06,     3.11267859150471d+06,
     &    -2.07773763241244d+06,     6.64189638782395d+05,
     &    -3.65852023832843d+03,     3.17801643790498d+04,
     &    -1.46495677415953d+05,    -2.74029395458975d+04,
     &     9.38040610889768d+03,    -1.50556238520670d+05,
     &     8.18452448145059d+05,    -3.42642760904873d+05,
     &    -3.68695164374142d+06,     9.24733594147077d+06,
     &    -8.43317214713195d+04,     2.56605669723925d+05,
     &    -4.61514472615555d+05,     1.59757159300529d+05,
     &    -9.25910121757708d+04,     4.55096190910857d+05,
     &    -1.78209721412275d+05,     1.03795929064392d+06,
     &     6.44134184843691d+05,     2.72855898296173d+06,
     &    -8.73968779765680d+06,    -1.03511587931790d+07,
     &     9.13784660459701d+06,     1.81985062849053d+06,
     &     3.11058160563426d+06,     5.95801906260502d+05,
     &     2.33140285840864d+06,    -6.28547071753558d+06,
     &     1.69915131304438d+07,    -2.43817175876122d+07,
     &     1.61180390446982d+07,    -3.88837195047505d+06,
     &    -2.19605104428893d+06,     4.48171362860180d+06,
     &    -4.68537074743909d+06,     2.12830590064047d+06,
     &    -5.73533100058749d+05,     2.34874239031607d+06,
     &    -1.04932770544528d+07,     2.11263175116850d+07,
     &     4.28218603150160d+07,    -1.98529351102832d+08,
     &     1.05839470731122d+05,     3.43452064966491d+04,
     &    -6.26431896268560d+04,    -3.37259470729340d+05,
     &    -4.28774458583546d+05,     1.88698109132963d+06,
     &    -9.83202836688817d+05,     1.52072570928178d+06,
     &    -9.47812719774204d+06,    -1.41145187802000d+07,
     &     3.97390312908862d+07,     3.61873915798496d+07,
     &    -3.45724416783427d+07,    -7.50680770339860d+06,
     &    -1.08253334955479d+07,    -4.18204276702211d+06,
     &     4.31459105400095d+06,     2.50019113548640d+06,
     &    -5.34775882036780d+07,     9.99887378530038d+07,
     &    -7.18750635375827d+07,     2.08996463546660d+07,
     &     8.40074763916543d+06,    -3.60739631543694d+07,
     &     4.74212580720975d+07,    -2.39760782815832d+07,
     &     1.62163658489909d+07,    -7.24043826659967d+07,
     &     1.94406354104845d+08,    -2.08423218035645d+08,
     &    -5.64713691078299d+08,     2.70735468283329d+09,
     &    -5.77883647109710d+04,    -7.92608405227715d+05,
     &     2.55116431479687d+06,    -1.36897523529383d+06,
     &    -3.67829068807010d+06,     4.91628257286565d+06,
     &    -6.32893072989524d+06,     1.21704728616145d+07,
     &     2.24294316804479d+06,     8.80909024668154d+07,
     &    -1.67798279317535d+08,    -1.16344041122537d+08,
     &     1.24564403289980d+08,     2.71012681487803d+07,
     &     3.25393683011655d+07,     3.06438234102462d+07,
     &    -8.84260784246285d+07,     1.75532917746713d+08,
     &    -1.40773769909268d+08,     1.91969791395780d+07,
     &     4.44778446548285d+07,    -3.62786892000888d+07,
     &    -4.52342823341440d+07,     2.26795783548777d+08,
     &    -3.44201717722045d+08,     1.81449595303069d+08,
     &     4.69337736639142d+07,    -6.11498780516947d+08,
     &     4.12927298017709d+09,    -2.02833142185066d+10,
     &     6.71934538758247d+10,    -1.47680003638032d+11,
     &     5.57440319768090d+05,    -2.38954887556639d+06,
     &     4.78355908250700d+06,     7.53920594663087d+05,
     &    -1.21276043725475d+07,     1.69330673208798d+07,
     &    -1.56610812177990d+07,     3.88400828400850d+07,
     &    -4.27793163007588d+07,    -3.48662030212352d+08,
     &     6.03293557624223d+08,     3.53468400977802d+08,
     &    -4.06555909705237d+08,    -1.20376845513766d+08,
     &    -1.94822210660193d+07,    -2.54668170209042d+08,
     &     5.01777380136094d+08,    -1.28673418904697d+09,
     &     2.80221394441643d+09,    -3.77705572245257d+09,
     &     2.54063151335745d+09,    -1.16359922666807d+09,
     &     1.21956560885115d+09,    -1.50300789073942d+09,
     &     1.63760191479698d+09,    -1.16938303087404d+09,
     &     2.05371930690318d+09,    -8.55457159786324d+09,
     &     2.41549301054165d+10,     1.74727620851866d+10,
     &    -6.06357645039771d+11,     2.34132810727480d+12,
     &    -6.71739259016827d+05,     3.65562666301699d+06,
     &    -1.17953667496232d+07,     1.57271065914434d+07,
     &    -1.72743636381687d+07,     3.19952853149412d+07,
     &    -4.05468437539810d+07,     3.57967947962238d+07,
     &    -1.60130346292695d+08,     1.27740453005925d+09,
     &    -1.91203274134508d+09,    -9.18661599153756d+08,
     &     1.07631308813097d+09,     5.78330648423612d+08,
     &    -3.72379933953311d+08,     1.30216588768490d+09,
     &    -1.05882333644062d+09,     4.34568666689789d+09,
     &    -1.78339211954918d+10,     3.09722009633639d+10,
     &    -2.40604915867143d+10,     1.39452275696760d+10,
     &    -1.33253636428437d+10,     9.14989118391608d+09,
     &    -5.37518909667066d+09,     4.51446579194586d+09,
     &    -1.50129065125501d+10,     1.69924239710544d+11,
     &    -1.62570853433930d+12,     7.74969138385510d+12,
     &    -6.26708251270892d+12,    -2.74305420767650d+13,
     &    -1.37539558419520d+06,     1.59898155802093d+07,
     &    -4.80186171959022d+07,     3.73423962830080d+07,
     &     8.03473879625904d+06,     2.07021469437228d+07,
     &    -1.34118266247588d+08,     2.63529508841205d+08,
     &     1.85159958570545d+08,    -2.85763441846553d+09,
     &     4.31545394647982d+09,     1.80206201969196d+09,
     &    -1.90993318814252d+09,    -1.98009742129264d+09,
     &     2.00938449583032d+09,    -3.79973887489561d+09,
     &    -1.70894817521891d+09,    -5.72813176723302d+09,
     &     6.18157277069104d+10,    -1.22772293466389d+11,
     &     9.87595131493905d+10,    -6.19142006274281d+10,
     &     6.54538431473159d+10,    -3.60133171094614d+10,
     &    -6.11279988445065d+09,     4.36275019167166d+09,
     &     1.98874751069926d+11,    -3.03932292155295d+12,
     &     2.90115486491873d+13,    -1.51902239254655d+14,
     &     1.92605617641509d+14,     2.43182755549621d+14, ! id_smallr = 7, smallr = 2.45 bohr
     &     1.04206274929206d+04,    -9.14604860618583d+03,
     &     1.35136873218920d+03,    -4.21510207881840d+03,
     &     8.97982663164064d+03,    -1.74099431209322d+04,
     &     3.35606183731877d+03,     2.37092460444146d+04,
     &    -1.11928963851341d+04,    -2.09284922789622d+04,
     &    -1.90815888713805d+04,     3.25636681460977d+04,
     &     1.10273843886676d+05,     2.45422525100582d+03,
     &     2.97417555615293d+04,    -6.48899581040134d+03,
     &     1.06001941524374d+04,    -1.57866701514656d+04,
     &     1.26589979019278d+04,    -6.57135227351392d+04,
     &     1.26297106554975d+04,    -6.82621243977560d+04,
     &    -4.94231755442178d+04,    -4.90749797530958d+04,
     &    -3.15022800851246d+04,    -2.96554463454368d+04,
     &     3.77442065505890d+04,     6.35087952558070d+04,
     &     2.50813475195732d+04,    -3.11388599361001d+03,
     &    -4.78096507483204d+03,    -1.48816794450232d+03,
     &    -1.36014807293312d+04,     2.54228815981720d+03,
     &    -1.87608810227896d+04,    -2.34422083811707d+04,
     &    -6.57131966856908d+03,    -1.02632575072567d+05,
     &     8.91348852829950d+04,     3.22190918439009d+05,
     &     2.41446557487497d+05,     2.37392380353072d+05,
     &     3.95189570991293d+05,    -1.62974968733953d+05,
     &    -1.22522013084712d+06,     4.71003606761481d+04,
     &    -3.68484598196853d+05,     4.87475760771003d+03,
     &    -2.18695954921343d+05,    -8.31470934176952d+04,
     &    -1.23893152916574d+05,     7.99915722902941d+04,
     &    -1.52064159444851d+05,     1.27965052511074d+05,
     &     2.54572922132308d+05,     3.89434514727248d+05,
     &     2.40879011462456d+05,     7.06826331207238d+05,
     &     3.17938389191571d+05,    -4.73022974695579d+05,
     &    -3.83108338233472d+05,    -4.56663051634451d+03,
     &    -1.24587737321389d+05,     1.89035538591101d+05,
     &     3.11876994506780d+04,    -2.77382571598745d+04,
     &     2.18787898621403d+04,    -9.95860952376803d+03,
     &    -2.28674606225445d+04,    -4.51675431293246d+04,
     &     1.31750192801458d+03,     7.48172547292701d+04,
     &    -3.52340142334869d+04,    -1.27652565169930d+05,
     &    -1.03054117248488d+05,     6.08103445196746d+04,
     &     3.98205208243683d+05,    -5.86238468411525d+04,
     &     5.86472803781500d+04,    -5.36319547475072d+02,
     &    -7.64908938854372d+03,     3.97536961689276d+03,
     &    -1.04705840003223d+04,    -2.29465899708710d+04,
     &    -2.25658697883925d+04,    -7.16002487977251d+04,
     &    -6.40896017103730d+04,    -7.84464937506044d+04,
     &    -4.68367784302112d+04,    -6.39208842912892d+04,
     &     2.98366753646184d+04,     9.83224269660056d+04,
     &     1.78657398273263d+04,     1.64777346140917d+04,
     &    -1.37572879532478d+04,    -8.18917209079322d+03,
     &    -2.52819006630598d+04,     1.75151388566060d+04,
     &     3.26094895370269d+03,     1.23205602262278d+04,
     &    -7.82341333188100d+04,    -1.62636220882398d+05,
     &    -9.07498769167725d+04,     3.89221499874686d+05,
     &     1.63055575094275d+05,     3.06168882036799d+05,
     &     5.94489447845408d+05,    -1.20917292715143d+05,
     &    -1.83578213422087d+06,     4.33548739739667d+05,
     &    -1.33629435324150d+05,    -2.38444174214107d+04,
     &     1.07008782855707d+05,    -1.59642606129068d+05,
     &     3.48040726597758d+04,    -4.16368744061774d+04,
     &    -3.40890452453646d+04,     9.36490759531512d+04,
     &     7.89973025874749d+04,     2.81201675003341d+05,
     &     1.50048321355298d+05,     4.63302033331658d+05,
     &     1.90544430935805d+05,    -4.00121176761626d+05,
     &     3.36278812800706d+03,    -3.39203472175922d+05,
     &     2.63271395699286d+05,    -1.56989178179520d+05,
     &     1.27747448248940d+05,    -1.10791274966364d+05,
     &     3.90786690136834d+04,     1.03079659043475d+05,
     &     2.34142022450198d+05,    -5.64376231314625d+05,
     &    -7.08266556939549d+04,     8.47699496892945d+05,
     &     7.36732566782390d+05,    -2.24196869377475d+06,
     &    -1.69849201963723d+06,     6.28718274502509d+04,
     &     7.96423780969935d+06,    -2.87342050644026d+06,
     &     1.50327262137042d+05,    -1.49777684540366d+05,
     &    -7.04052965886418d+05,     5.22984159633035d+05,
     &     1.11123735953358d+05,    -1.25335345323453d+05,
     &     3.73774271536474d+05,    -2.41221497589608d+05,
     &    -2.82645076306279d+05,    -9.50924991403924d+05,
     &    -5.02287656542990d+05,    -2.41901993438922d+06,
     &    -1.41852710788032d+06,     1.73033877846420d+06,
     &     9.22282047924930d+05,     1.92445104271128d+06,
     &    -4.58031986925861d+06,     5.22633628609894d+06,
     &    -1.28727503570306d+05,     9.30333798479056d+04,
     &    -2.33623740611166d+05,    -6.15594426753064d+04,
     &     1.55614844330271d+06,    -1.62917840546722d+06,
     &    -7.02932070510655d+05,     6.71429191629092d+05,
     &     1.18548182400966d+06,     4.81970785203182d+05,
     &     6.56435658035616d+06,    -2.74833832587603d+06,
     &    -3.20285871052147d+07,     1.67073098530463d+07,
     &     7.51915627564442d+05,     4.33355000618486d+06,
     &     1.60817474829482d+06,     3.24892142093893d+06,
     &    -4.80246338041153d+06,     4.03405066299024d+06,
     &    -2.27361551793651d+06,    -1.96423467872370d+06,
     &     6.14443214304059d+05,     9.69636121236143d+05,
     &    -2.18340435470410d+06,     9.01181702244351d+06,
     &     5.89031009111672d+06,     7.12428166451992d+06,
     &    -3.75291931707471d+07,     5.46875543322870d+07,
     &    -3.42890929057352d+07,    -3.33497726215510d+07,
     &     9.09022858863712d+04,    -9.21932481104890d+04,
     &     6.04690331035280d+04,    -1.61604987831537d+05,
     &     2.00307355732445d+05,     7.79513021500732d+04,
     &    -3.81988092796217d+05,     1.27795907180386d+05,
     &     1.03712800951201d+06,    -1.77564127264808d+06,
     &    -5.65613111500384d+05,     1.96416285698863d+06,
     &     5.70549916968373d+06,    -4.51281990831408d+06,
     &    -4.10971644492889d+05,    -1.44332377593092d+06,
     &     1.06998727795801d+05,    -1.00656242756513d+06,
     &     1.09970673458319d+06,    -5.55230923282463d+05,
     &    -2.51140738809382d+04,     4.06332086586938d+05,
     &    -1.33272762195784d+04,    -3.36165087569762d+05,
     &     6.91064740578351d+05,    -5.94975191095429d+05,
     &    -1.68447851456138d+04,    -2.76618365251326d+05,
     &     2.20539679098989d+06,    -5.51025597803615d+06,
     &     1.07518442374760d+07,    -9.38543111747333d+06,
     &    -1.08793938334112d+05,     1.00829860907089d+05,
     &     1.36839098862668d+05,    -2.80065932257275d+05,
     &    -7.96620829986311d+05,     9.66627872235460d+05,
     &     7.15149421808552d+05,    -1.88212429849762d+06,
     &     5.37077398712691d+06,    -3.24661765340911d+06,
     &     7.67781125723021d+06,    -1.65003312380575d+07,
     &    -1.77666714149022d+07,     1.59191102953097d+07,
     &     7.23716168035092d+06,     5.99918410722731d+05,
     &     2.92979024785226d+06,     4.12777394844683d+06,
     &    -4.43923564055327d+06,     3.94824957984581d+05,
     &     2.13506890243232d+06,    -1.67732391463964d+06,
     &    -1.71182343829721d+06,     2.48818974970667d+06,
     &    -3.33823915501946d+06,     1.53389990688622d+06,
     &     1.22526099337523d+06,    -6.13813382359387d+06,
     &     1.49046015767736d+07,    -1.24056573993165d+07,
     &    -3.82523123736042d+07,     1.27675002164713d+08,
     &     2.52918417760833d+05,    -2.85937671923729d+05,
     &     3.29231053322453d+05,     3.60841317017205d+05,
     &    -1.57524961132381d+06,    -1.52283728631584d+06,
     &     8.23993372267139d+06,    -1.24590843412122d+07,
     &     1.53364361932408d+07,    -2.11174835804851d+07,
     &    -3.75281895246736d+07,     9.66544560417397d+07,
     &     4.48369680652666d+07,    -4.50755867680748d+07,
     &    -4.33852887732326d+07,     6.12297987721588d+06,
     &    -2.41856815006134d+07,     9.31945536114990d+06,
     &    -5.68698370777718d+06,     1.80218966731781d+06,
     &    -5.66351420619537d+05,     3.94711873044068d+06,
     &    -3.90137390327838d+06,     1.87845287801250d+07,
     &    -2.32608809789956d+07,     1.18754936201997d+07,
     &     3.23941847319891d+06,     5.47748138976546d+07,
     &    -3.92114825740338d+08,     1.18501591553033d+09,
     &    -2.47045892682039d+09,     2.49668219595969d+09,
     &    -3.28640704990859d+05,     3.82008459290885d+05,
     &    -1.52830074838271d+06,     4.08409342448587d+06,
     &    -4.72041842251281d+05,    -1.19691330688845d+07,
     &     2.54519262825174d+07,    -5.74334907451336d+07,
     &     9.47293413859559d+07,    -6.27457985435333d+07,
     &     2.63793414520360d+08,    -4.70690309410391d+08,
     &    -9.09269857225495d+07,     1.09111358870075d+08,
     &     2.02830294014182d+08,    -1.62601414483764d+07,
     &     8.10422414972262d+07,    -1.03848508150190d+08,
     &     1.02003589966032d+08,     3.43909278977982d+06,
     &    -3.65965113756034d+07,    -8.02859916266549d+07,
     &     2.29774477650109d+08,    -3.99005472128282d+08,
     &     4.30280224115674d+08,    -2.47827638385258d+08,
     &     1.08812805421686d+08,    -2.87814493677327d+08,
     &     1.81131729963264d+08,     2.01260250220090d+08,
     &     2.67006814538925d+09,    -6.23458118768780d+09,
     &     6.51858427712505d+05,    -3.87949548046539d+05,
     &    -2.28913066168680d+06,     2.44131534836829d+06,
     &     2.21562109594827d+07,    -7.02273122854323d+07,
     &     1.37815003395484d+08,    -2.31672143700459d+08,
     &     3.30637633035655d+08,    -2.06752270835237d+08,
     &    -1.02213651405648d+09,     1.90795878729729d+09,
     &     9.78747204632041d+07,    -3.04571746786195d+08,
     &    -6.16886870713159d+08,    -1.44816328139942d+08,
     &    -1.48248260186486d+08,     3.42239707043846d+08,
     &    -1.28549795201666d+08,    -4.00281728460463d+08,
     &     3.29939318804336d+08,     4.12558831775080d+08,
     &    -1.07826895180952d+09,     1.57737991449077d+09,
     &    -1.80524255626770d+09,     1.18355557274678d+09,
     &     8.33694870014507d+08,    -1.23720029783983d+10,
     &     5.01375275731128d+10,    -7.97459317947880d+10,
     &     8.96399140037605d+10,     5.59448775260646d+10,
     &    -7.39298369985899d+05,     2.17104341706756d+05,
     &     6.29121777063382d+06,    -3.14332998042718d+07,
     &     7.40968421864410d+07,    -1.76564594790955d+08,
     &     4.24130092952244d+08,    -7.00662661070652d+08,
     &     6.71983449869477d+08,    -9.45123253870642d+08,
     &     4.42289713080810d+09,    -6.76125641697078d+09,
     &     1.13848239968532d+08,     1.04294031433164d+09,
     &     1.44592050217707d+09,     1.25926540430530d+09,
     &    -8.17774421962984d+08,     1.99959073719920d+09,
     &    -4.26372696468756d+09,     3.49213660052169d+09,
     &    -8.30347898516746d+08,     3.78740287325265d+08,
     &    -3.32431051145108d+09,     5.46730602299758d+09,
     &    -1.97142442307396d+09,    -1.92500172029854d+09,
     &     5.59427384675990d+09,    -2.79397788096921d+09,
     &    -3.23121597925172d+11,     3.05170391313235d+12,
     &    -1.64610526311110d+13,     3.66492691819482d+13,
     &     9.29086431768325d+05,    -2.89368634573997d+06,
     &     2.80983031250213d+07,    -6.50531226107304d+07,
     &    -6.14268557655715d+07,     1.90123358904071d+08,
     &     2.78115778210022d+08,    -1.39729953716074d+09,
     &     2.42996585026476d+09,    -4.71057570884230d+08,
     &    -9.44662765583472d+09,     1.58111646263446d+10,
     &    -1.92777832888883d+09,     1.18520308731812d+09,
     &    -7.28123033446976d+09,    -1.76429057602654d+08,
     &     2.05781816467235d+09,    -1.43983009613247d+10,
     &     2.46244265110307d+10,    -6.05004124935805d+09,
     &    -1.01367506695496d+10,    -7.48180682386149d+09,
     &     4.90752142918965d+10,    -8.39024444499231d+10,
     &     9.04928440730062d+10,    -7.37958381134225d+10,
     &     1.29692901098575d+11,     1.91838800644543d+11,
     &    -4.98915226682674d+12,     3.03256510104100d+13,
     &    -1.30031442676772d+14,     1.75937412979208d+14, ! id_smallr = 8, smallr = 2.65 bohr
     &     9.72044555847233d+03,    -5.14717592816509d+03,
     &     3.84828523373061d+03,    -1.92957438650008d+04,
     &     2.19924798145502d+04,    -1.41712510947253d+04,
     &    -1.94919844318782d+04,    -9.65175780502850d+03,
     &     5.00838383137181d+04,     3.44782820137938d+04,
     &    -1.44136065287130d+05,     5.20886085881124d+04,
     &     2.96225428982140d+04,     1.62012667726953d+05,
     &    -6.65842789774140d+03,     2.71224450157110d+04,
     &    -6.75823635143159d+03,    -8.61300984333321d+03,
     &    -3.48240179271780d+03,    -3.78732942231466d+04,
     &    -1.58669789407097d+04,    -6.67879055386572d+04,
     &    -5.54406318563386d+04,    -5.47435335471415d+04,
     &    -3.69169197593310d+04,    -3.47272001663219d+04,
     &     3.67295029537415d+04,     7.10890776959356d+04,
     &     2.21075889256179d+04,     8.20253444396214d+02,
     &    -1.39739585400981d+04,     1.65736271482561d+04,
     &    -1.40865911450494d+03,    -3.12795713336249d+04,
     &     1.86159464403057d+04,    -5.98701492166753d+04,
     &     5.56291572815183d+03,     7.84601892794301d+03,
     &    -2.74080343446754d+05,     3.55857146794916d+04,
     &     1.36430471581436d+06,    -1.07197585721340d+06,
     &     1.82404594838767d+06,    -2.89711751347965d+05,
     &    -2.01296800061630d+05,    -1.91537943935335d+06,
     &     5.72465800524033d+04,    -4.54220495686052d+05,
     &     1.03200163964617d+05,    -3.67692092330793d+05,
     &    -6.04355803369646d+04,     5.73194976061470d+05,
     &    -7.11792366651613d+05,     4.13087314497367d+05,
     &     1.83049129531711d+05,     4.65630923051115d+05,
     &     3.02767571259987d+05,     7.70487849904550d+05,
     &     3.67097158126992d+05,    -5.53669154277266d+05,
     &    -3.41659519450463d+05,    -1.20822275816262d+05,
     &     2.30831590471711d+04,     6.57174236958615d+04,
     &     2.80716699937914d+04,    -8.88074280270409d+03,
     &    -7.05444653968559d+03,     4.42635380748863d+04,
     &    -6.96780935141685d+04,     1.28918944476331d+04,
     &    -6.92441758690226d+04,    -1.60553596708696d+05,
     &     4.31112639282789d+05,    -3.02781422662874d+05,
     &    -1.17431780880510d+05,    -1.60335101998272d+05,
     &     1.67940545578425d+05,     5.31324641044030d+05,
     &    -5.42574136653109d+04,     4.39832143127967d+04,
     &    -1.94820400935638d+04,    -1.27726838906326d+04,
     &     2.01122885631894d+04,    -1.27629614683881d+05,
     &     4.92763641671511d+04,    -1.12412538715973d+05,
     &    -7.46944927178940d+04,    -8.71875610515178d+04/
      DATA coe_RKHS(3001:3744)/
     &    -6.28639444001137d+04,    -7.29565324712267d+04,
     &     3.33473386426839d+04,     1.08963848899261d+05,
     &     2.97352739499416d+04,     6.25471555857741d+03,
     &    -8.58304288485999d+03,    -9.36792286355168d+03,
     &    -3.30673193814560d+03,    -4.76443519423121d+04,
     &     2.33850463562844d+04,     8.11166562913259d+04,
     &    -1.15719898112373d+05,    -1.30014069907810d+05,
     &     2.69488928277226d+05,    -1.60326713579214d+06,
     &     2.49028380043364d+06,    -5.46719936963441d+05,
     &     2.44890077345800d+05,     1.41671697027183d+06,
     &    -7.98568647736552d+05,    -2.42105216902088d+06,
     &     3.17263964179425d+05,     5.50096807296087d+04,
     &    -5.24221538599623d+04,    -1.10781090723999d+05,
     &     9.68763257323777d+04,     2.64591630899188d+05,
     &    -3.47254962791394d+05,     1.87073585048355d+05,
     &     1.60665632392160d+05,     2.86424858595202d+05,
     &     2.15688634754549d+05,     5.55622561709798d+05,
     &     2.39438590128742d+05,    -4.41519365717452d+05,
     &    -1.85672912640459d+05,    -1.07317564301798d+05,
     &     2.16391543214752d+04,     6.26505493170820d+04,
     &     1.03348164817748d+05,     2.94989438942131d+04,
     &    -1.65060327334811d+04,    -2.04371230115202d+05,
     &     7.03226825473680d+05,    -3.47501498976131d+05,
     &     1.05316471278232d+06,    -5.87890454591728d+06,
     &     1.17230437851273d+07,    -1.19716166851314d+07,
     &     8.05124212629464d+06,    -1.05952936755068d+07,
     &     4.73108400511129d+06,     9.57316889656268d+06,
     &    -2.31121467949232d+06,    -1.52564145981735d+06,
     &     8.63671701676347d+05,     1.34149884504881d+04,
     &    -3.41989746514587d+06,     6.13914821547494d+06,
     &    -4.28112236334999d+06,     1.46301651563576d+06,
     &    -1.13305209407952d+06,    -1.20493908902318d+06,
     &    -4.81115755455239d+05,    -3.11603071835104d+06,
     &    -2.04901700832914d+06,     2.41008075144399d+06,
     &     1.18055646438779d+06,     2.83637473146920d+05,
     &     2.99716900641919d+06,    -5.22447658579099d+06,
     &    -2.42151517977593d+03,    -3.50276188912380d+05,
     &     4.71981059822440d+05,    -1.55083084923519d+06,
     &     2.45424751912994d+06,    -3.63525194051092d+05,
     &     2.68436554290110d+06,    -1.54948591955060d+07,
     &     1.96915968680164d+07,    -2.31058978670686d+06,
     &    -1.30793929386296d+07,     2.95730737741518d+07,
     &    -1.97473232319363d+07,    -3.88855711661338d+07,
     &     2.08071912423789d+07,     2.38540658359714d+06,
     &     7.02547320978897d+06,     4.25701862592395d+06,
     &     2.85757187560080d+06,    -2.22775752745514d+07,
     &     2.24864793697914d+07,    -1.12635749330278d+07,
     &     3.67188646958119d+05,     6.15241666273814d+06,
     &    -8.92311340963224d+06,     1.49796737354982d+07,
     &     9.46705084724984d+06,    -3.98086437806923d+06,
     &    -8.44618300988101d+06,    -2.61966474999127d+06,
     &     1.42684645119886d+07,    -2.56461910450117d+07,
     &     5.36025388588645d+04,     7.37589399048933d+04,
     &    -1.06041359918517d+05,     9.77700674194953d+04,
     &    -4.07445544067904d+05,     3.60784108183137d+05,
     &     1.11254094285205d+06,    -3.71844308505102d+06,
     &     5.09250614696182d+06,    -2.07985633029343d+06,
     &    -1.03832964716101d+06,    -2.99642596550986d+06,
     &     6.04636098013808d+06,     5.60180359320381d+06,
     &    -5.59940437896457d+06,    -1.08534082635391d+06,
     &    -1.23061514501673d+06,    -2.33416804099457d+06,
     &     3.31246383901181d+06,    -2.10001840721375d+06,
     &     4.81080250420312d+05,     3.69471410759486d+04,
     &     8.53083032408936d+05,    -1.37087897295435d+06,
     &     1.78329427158353d+06,    -1.35768107891566d+06,
     &     4.05769113861655d+04,     1.64296656918585d+05,
     &     5.29225538135136d+05,    -1.47661123117267d+06,
     &     5.73121179861272d+06,    -9.28427123036714d+06,
     &     1.01455594704989d+04,    -2.70594459727412d+05,
     &     2.18075284944116d+04,     1.63701785255697d+06,
     &    -4.04787142150042d+06,     3.17877817087885d+06,
     &    -1.94598779763366d+06,     1.50239063005612d+05,
     &     1.33782464809337d+07,    -2.80621048876872d+07,
     &     3.29619311097265d+07,    -4.95690242172476d+06,
     &    -3.27759300253362d+07,    -1.43497414677803d+07,
     &     2.29505011687946d+07,     7.09441792080521d+06,
     &    -1.15037173003594d+06,     2.97070854151261d+07,
     &    -3.75918807531545d+07,     2.45752858412573d+07,
     &    -7.94515047759746d+06,     1.15112356303461d+06,
     &    -5.08224033731485d+06,     7.30559496649521d+06,
     &    -1.09554379974544d+07,     6.63840016477591d+06,
     &    -5.35986308282253d+05,    -4.83265579205514d+05,
     &     4.67753735290590d+05,    -6.29075591522995d+06,
     &     3.19142367511452d+07,    -5.78085486122486d+07,
     &     6.16581450605128d+04,     5.45459067875531d+05,
     &    -9.28960545935468d+05,     2.09776403188996d+06,
     &    -1.84246694323327d+06,    -3.37663735941939d+05,
     &    -1.05040369264260d+07,     2.64208061630357d+07,
     &    -4.45849967883408d+07,     1.14721158607337d+08,
     &    -1.76485551587048d+08,     1.54680628239958d+07,
     &     1.92148787403208d+08,     1.86281259148366d+07,
     &    -7.79083240831099d+07,    -4.40382875454670d+07,
     &    -2.68857648205580d+06,    -6.57515786886586d+07,
     &     7.06271600258354d+07,    -5.62702808290356d+07,
     &     3.28376937060457d+07,    -1.04019371815402d+07,
     &     1.60069107783896d+07,    -1.74527129869349d+07,
     &     3.30550206249897d+07,    -1.97467002972384d+07,
     &     3.82244094179127d+06,     2.15496491465505d+06,
     &    -2.37395320437583d+07,     1.79855689352381d+08,
     &    -6.09614314468593d+08,     3.07207367178081d+08,
     &     1.63285098684335d+05,    -1.55473464786531d+06,
     &     3.29661257826990d+06,    -1.38600333101208d+07,
     &     3.85515462457193d+07,    -4.94913236276571d+07,
     &     4.09814603014869d+07,    -6.82424005079673d+07,
     &     1.27952820495318d+08,    -1.77493879052877d+08,
     &     2.22481330190533d+08,     3.71445010990640d+08,
     &    -1.09355042182586d+09,     1.30780109910042d+08,
     &     1.80233298415705d+08,     4.20232081067297d+08,
     &    -4.32285386276749d+08,     4.30809842340473d+08,
     &     5.15698803920425d+08,    -1.25512929330870d+09,
     &     8.55066594775735d+08,    -3.32792198682382d+08,
     &     2.72752034601732d+08,    -3.58365689332269d+08,
     &     4.13179377306290d+08,    -2.84733963540427d+08,
     &     1.07887246813201d+08,    -2.19356904859986d+08,
     &     6.94968170111614d+08,    -1.76198730686650d+09,
     &     4.84505094708177d+09,    -1.05006319050409d+10,
     &    -4.88053644465749d+03,     1.75532460105083d+06,
     &     2.56387419775887d+06,    -2.96433765707854d+07,
     &     6.90145890219916d+07,    -3.56528451136348d+07,
     &     1.06044146972151d+07,    -2.09673979884959d+08,
     &     6.71277337389768d+08,    -1.17415023257602d+09,
     &     1.45272391150105d+09,    -3.44589410158659d+09,
     &     5.28449647961938d+09,    -1.11051384138976d+09,
     &    -4.16846265547333d+08,    -1.49230592499976d+09,
     &     1.03802406967480d+09,    -6.98727151259930d+07,
     &    -4.91727847999000d+09,     9.15199085140067d+09,
     &    -6.59494410727685d+09,     2.61862753440113d+09,
     &    -1.99547029203908d+09,     3.45654032968210d+09,
     &    -5.76322466262388d+09,     4.33184729803710d+09,
     &    -1.28373445721452d+09,     3.81564700387578d+08,
     &    -5.11106253633565d+09,     8.52424522859801d+10,
     &    -7.80781732132389d+11,     3.00045193283276d+12,
     &     4.15001168952721d+05,    -3.12316585651893d+06,
     &    -4.17514074395616d+06,     5.90200933424806d+07,
     &    -2.18051374229687d+08,     3.75910710354967d+08,
     &    -3.08961354649420d+08,    -6.34446014575161d+07,
     &    -1.24506973559258d+09,     6.94760725987761d+09,
     &    -1.22637178365284d+10,     1.76787029229210d+10,
     &    -2.02383878250508d+10,     4.41744661134254d+09,
     &     2.30806317590382d+09,     7.43242893368271d+08,
     &     6.63704323125095d+09,    -5.79141541548887d+09,
     &    -3.62664009765472d+09,     1.03926106869789d+10,
     &    -6.01600257235217d+09,     2.13968985438809d+09,
     &    -3.78906260505081d+08,    -1.58177487251655d+10,
     &     4.11254591957347d+10,    -3.43240656336536d+10,
     &     1.30059059591586d+10,    -6.52122082609106d+09,
     &    -8.91056108369834d+10,     1.05289428434914d+12,
     &    -6.69008355816733d+12,     1.55728978902083d+13,
     &    -7.55717381468374d+05,     8.89734080035271d+06,
     &    -5.39451820983505d+07,     3.48320105600135d+08,
     &    -9.38839856315153d+08,     5.62302456479449d+08,
     &     1.91233039108858d+09,    -5.90116696775888d+09,
     &     1.32844742930494d+10,    -2.81863148150755d+10,
     &     4.09106260471570d+10,    -5.13875620577740d+10,
     &     5.19510078180948d+10,    -1.03251512010222d+10,
     &     1.50785737722978d+09,    -3.50939785512815d+10,
     &     4.85759372639594d+10,    -6.96634057702831d+10,
     &     1.04020589860937d+11,    -1.14251340081436d+11,
     &     6.89528613725225d+10,    -4.28070961353122d+10,
     &     4.04584383999941d+10,     4.34566415008956d+10,
     &    -1.92372399353795d+11,     1.79323297913549d+11,
     &    -8.29521345922714d+10,     1.00619743688137d+11,
     &    -1.44550002396410d+11,    -6.60761850553953d+12,
     &     1.20951488757394d+14,    -6.43581334905807d+14, ! id_smallr = 9, smallr = 2.85 bohr
     &     1.07045308506576d+04,    -4.92634593488883d+03,
     &     1.92128873075423d+03,     3.22174925500901d+03,
     &    -3.68268912065506d+04,     3.94701485269034d+04,
     &    -2.64870155030266d+02,    -1.13669202764497d+05,
     &     2.33996962468574d+04,     2.97012059973824d+05,
     &    -4.63256156883592d+05,     3.34762451411022d+05,
     &    -2.65686462713573d+05,     2.22493267547699d+05,
     &     1.28771371312493d+05,     1.61944143563463d+04,
     &     2.33970494685265d+04,    -4.50271705584136d+04,
     &    -3.77766642701829d+03,     2.63073778086173d+04,
     &    -9.40799094029331d+04,    -4.57946799962688d+04,
     &    -7.12530161794732d+04,    -6.50839062964239d+04,
     &    -2.99517592915735d+04,    -5.08855317816419d+04,
     &     3.59827760573479d+04,     7.50355596361055d+04,
     &     2.09287472768485d+04,     3.13067185738868d+03,
     &    -6.12109098205913d+04,     1.16845193512827d+05,
     &    -4.89382419807663d+02,    -3.39397040366603d+04,
     &     4.43348736483969d+03,    -5.14238454705958d+03,
     &    -9.34186567571975d+04,     4.44002210789572d+04,
     &    -6.63678060318215d+04,    -1.32511224083054d+05,
     &    -6.44312109682316d+05,     2.44043506213062d+06,
     &    -8.40551071196746d+03,    -1.95881441574714d+05,
     &     1.70909245093881d+06,    -1.92583427906429d+06,
     &    -1.84721612261678d+06,    -4.12486657433019d+05,
     &    -7.05428158842033d+04,    -2.68229723182455d+05,
     &    -1.50429882392690d+05,     4.88452498628350d+05,
     &    -4.52090332788360d+05,     3.34131380175525d+05,
     &     2.54024707302666d+05,     6.51950574237615d+05,
     &     8.69808032610355d+04,     1.05193574078310d+06,
     &     3.85656612604191d+05,    -5.88807026140518d+05,
     &    -4.09706017938220d+05,    -5.60179046859834d+04,
     &    -1.38451315697011d+05,     2.63835392803484d+05,
     &     3.46976978682864d+04,    -1.45243001143688d+04,
     &     6.84649123656527d+03,    -1.24621178246974d+04,
     &     8.12714653810400d+04,    -1.32633594173763d+05,
     &    -1.32575658999885d+04,     7.38832824954160d+04,
     &    -7.47907024338985d+05,     1.45506236074758d+06,
     &    -1.22437437342897d+06,     3.03612232548181d+05,
     &    -4.16013379682511d+05,     4.91288386664078d+05,
     &     5.10803723003891d+05,     2.78593858897374d+04,
     &    -5.40199537928349d+04,     8.27724754279000d+04,
     &    -1.92883624649857d+05,     1.19444541941198d+05,
     &    -1.40103340763636d+05,    -8.10858600904296d+04,
     &    -1.07400633442335d+05,    -1.14577693682396d+05,
     &    -4.55160370955968d+04,    -1.09091529030660d+05,
     &     3.94493728711078d+04,     1.22689339434522d+05,
     &     3.67727461514396d+04,    -1.24903619783857d+03,
     &     1.00199622791842d+04,    -3.23118083418747d+04,
     &    -6.24416362625767d+03,    -5.62361541203612d+04,
     &     4.32976583123600d+04,    -6.27507560857438d+04,
     &     2.22003732172824d+05,    -3.27624581812998d+05,
     &     4.09445873730697d+04,     3.95160537799278d+04,
     &    -2.58791664087516d+06,     5.75629127889164d+06,
     &    -4.70623979238943d+06,     5.07239369345306d+06,
     &    -1.27420408234469d+06,    -6.83125081160719d+05,
     &    -2.89307555904950d+06,     1.22103945880948d+05,
     &     8.62062914325855d+04,     1.25349054883717d+05,
     &    -2.03703388921034d+05,     2.24430077148546d+05,
     &    -2.59377227763097d+03,     1.14076826533794d+05,
     &     2.06002138088564d+05,     4.18160958293901d+05,
     &     1.68128887193059d+05,     7.33998810996148d+05,
     &     2.93271455352190d+05,    -5.24247881709811d+05,
     &    -2.18638608378067d+05,    -7.43143102893482d+04,
     &    -1.03955824167562d+05,     1.65028863576589d+05,
     &     1.56529571410920d+05,     1.35141834845134d+04,
     &    -5.86533521931747d+04,     1.75016964279687d+05,
     &    -6.48423619442930d+05,     1.77944264908454d+06,
     &    -1.45539947865480d+06,     3.89772552460514d+06,
     &    -1.13835565896113d+07,     9.13050000758744d+06,
     &     9.36552981161084d+06,    -2.37299926655735d+07,
     &     8.31680953625438d+06,     1.61773485866866d+06,
     &     1.04069535989376d+07,    -8.98063756204915d+05,
     &    -1.23368307272406d+06,    -6.00055601519123d+05,
     &    -1.05834859306726d+06,     1.67548385346722d+06,
     &    -1.03222165818152d+06,     3.48907853718638d+04,
     &    -7.83721650312069d+05,    -1.03478268416861d+06,
     &    -2.00369806138929d+06,    -3.01887773194401d+06,
     &    -2.92462977787600d+06,     3.10910382916958d+06,
     &     1.44569982094501d+06,     1.56464150941275d+05,
     &     2.49429717332579d+06,    -2.29257988348994d+06,
     &    -6.03767225215683d+04,    -3.35757683971530d+05,
     &     2.34799953024539d+05,     2.16616698905661d+05,
     &    -3.37608764491626d+06,     6.98114949225136d+06,
     &    -6.32545905485612d+06,     1.85438227539297d+07,
     &    -5.34427551753023d+07,     7.25758561649849d+07,
     &    -6.20193845658430d+07,     5.21046205551427d+07,
     &    -1.08005623384379d+06,    -3.70064682049799d+07,
     &    -2.29118846491395d+07,     1.08045698181938d+07,
     &     1.04395687594526d+07,     4.03127396513923d+06,
     &    -3.61028710813029d+06,     2.16567784024467d+07,
     &    -2.22976562786127d+07,     5.42681814995032d+06,
     &    -2.50134960726325d+06,    -7.59578538575634d+06,
     &     1.89533362292560d+07,    -1.17673512276274d+06,
     &     1.92995037868584d+07,    -1.04065066730651d+07,
     &    -5.83481160676614d+06,    -3.09151808149221d+06,
     &     2.28900072092214d+07,    -7.13392945803949d+07,
     &     1.17772015995389d+05,     2.94688035766887d+04,
     &    -1.07683062623717d+05,     6.60819530908928d+04,
     &     2.35163450845098d+05,    -1.41674897514877d+06,
     &     2.47293854626612d+06,    -3.60118094920557d+06,
     &     7.06690509914111d+06,    -1.01441195822880d+07,
     &     1.04929605170533d+07,    -5.02481506758169d+06,
     &    -9.63054541220774d+06,     1.99310536700429d+07,
     &    -5.13898747518143d+06,    -2.03795529877412d+05,
     &    -5.95233765259180d+06,     5.04947156761179d+05,
     &     2.37384854564350d+06,    -5.60278359977212d+06,
     &     3.55018018687538d+06,    -3.87372142950960d+05,
     &    -4.33366005779738d+04,     1.66292313974171d+06,
     &    -2.85296227928806d+06,     1.47956508018997d+06,
     &    -8.56063008173793d+05,     6.91681060561849d+05,
     &     4.02500171403132d+05,    -3.86111515085124d+05,
     &     3.25030993183530d+05,    -3.12019617841146d+05,
     &    -6.43600250291079d+04,    -2.99229328388303d+05,
     &     4.47334253405828d+05,    -9.60480391441334d+05,
     &     4.54825991116006d+06,    -1.11458775863257d+07,
     &     1.24804075796718d+07,    -1.66457559973603d+07,
     &     4.10808912496945d+07,    -8.42714505311698d+07,
     &     1.36931412237084d+08,    -1.52937806458652d+08,
     &     1.52858408915664d+08,    -1.66442648238468d+08,
     &     6.68246318376980d+07,    -7.01830293894865d+06,
     &     2.23456312624516d+07,     9.21413083581797d+06,
     &     1.95566032652335d+06,    -1.27489196028768d+07,
     &     1.24552952774888d+07,    -7.43627715728980d+06,
     &     4.28777899812465d+06,    -1.20801133945160d+07,
     &     1.51619490966254d+07,    -9.40008758634630d+06,
     &     3.85533663037773d+06,    -2.03970198731615d+06,
     &     2.39442623453784d+06,    -3.10666306214062d+07,
     &     1.04092840808910d+08,    -1.06448675260858d+08,
     &     3.42487661258370d+05,     3.14622980376477d+05,
     &    -4.41724723711726d+05,    -3.28782853649626d+05,
     &     2.89738364406817d+06,     1.57767177266688d+06,
     &    -1.72242099081939d+07,     1.93318800354856d+07,
     &    -6.14207071023351d+07,     2.12297501944039d+08,
     &    -3.57238849561818d+08,     4.15345748827689d+08,
     &    -6.08798315606751d+08,     8.32467296548388d+08,
     &    -3.84859594841622d+08,     6.00208117481835d+07,
     &    -1.40306173222445d+08,     6.57450933072673d+07,
     &    -1.34440046827143d+08,     8.59466877695391d+07,
     &    -1.33124840974162d+07,     1.17406447676602d+06,
     &     7.43048363835494d+06,     3.52063095799700d+06,
     &     1.05838121643449d+07,    -5.85654281672492d+06,
     &     3.98279633412202d+06,    -1.83249821452490d+07,
     &     4.50609439332947d+07,     5.22218540110837d+07,
     &    -3.88789342149617d+08,     4.52440068145183d+08,
     &    -2.35076863518645d+05,    -9.79436158295618d+05,
     &     4.90627217012977d+04,     5.53114852495680d+06,
     &    -3.73127638345691d+07,     1.10465750246926d+08,
     &    -1.27090510152293d+08,    -5.14307698320921d+07,
     &     3.58268106789120d+08,    -6.38643614417781d+08,
     &     9.27764662516252d+08,    -1.02507341391930d+09,
     &     1.94399872468592d+09,    -3.41215277863277d+09,
     &     1.79314280580773d+09,    -3.92982432093006d+08,
     &     7.78635298742982d+08,    -3.97988353369556d+08,
     &    -5.08530857483242d+07,     1.28045687715914d+09,
     &    -1.32874013907248d+09,     5.36351326915826d+08,
     &    -4.15171888002915d+08,     6.38379908302905d+08,
     &    -1.04657009294764d+09,     7.02233364666986d+08,
     &    -2.42634016373929d+08,     2.14270327528426d+08,
     &     4.87776057177489d+08,    -5.35553625968246d+09,
     &     1.29602502482984d+10,    -1.36488362852680d+10,
     &     8.39258747260060d+05,     1.52018271794717d+06,
     &    -1.27684819119470d+06,     1.19526359452162d+07,
     &    -8.07658887834965d+07,     1.65031083378984d+08,
     &    -6.84305456044053d+07,     2.36618561820244d+08,
     &    -1.20858305481646d+09,     1.41399204782961d+09,
     &     1.24888477223974d+09,    -4.57558692900679d+09,
     &    -2.48062253885696d+08,     9.89558493009753d+09,
     &    -6.60841951633598d+09,     2.59140435848979d+09,
     &    -3.17624683597151d+09,    -7.25575888111861d+08,
     &    -4.37117020739088d+08,     4.94209316973091d+09,
     &    -5.33414593917831d+09,     2.05220767513896d+09,
     &     9.85329962782616d+08,    -6.35783792156605d+09,
     &     1.24263476457013d+10,    -9.01422187503887d+09,
     &     3.07986599470434d+09,    -4.36898354267139d+09,
     &     2.10494594939064d+10,    -1.28029646247232d+11,
     &     4.90046883999979d+11,    -4.97859928694928d+11,
     &    -5.43943889845147d+05,    -3.14852588126601d+06,
     &     1.76999234102287d+06,    -2.26802555546021d+07,
     &     1.94506064281596d+08,    -6.56921491862979d+08,
     &     5.33900431237600d+08,     3.69990963449523d+09,
     &    -1.79105906971585d+10,     3.94270531207947d+10,
     &    -5.36177661315521d+10,     5.50443900227181d+10,
     &    -2.78307330842792d+10,    -1.63810645861148d+10,
     &     1.24047201950339d+10,     2.39423191922335d+09,
     &    -5.88922479012160d+09,     2.44349610201482d+10,
     &    -1.87156888239395d+10,    -9.11236400753899d+09,
     &     2.07769893775161d+10,    -9.48280052649906d+09,
     &    -8.48589992963073d+09,     4.86146035283123d+10,
     &    -9.81420453960359d+10,     7.44371240137385d+10,
     &    -2.30781478633190d+10,     2.14549434514658d+10,
     &    -2.23134143974438d+11,     1.79597582028546d+12,
     &    -6.38571512228888d+12,     5.51895965485135d+12,
     &     1.41468469586395d+06,     4.21767341207874d+05,
     &     1.96979206236864d+07,    -1.69506874411167d+08,
     &     1.09981426500287d+09,    -2.87237717329583d+09,
     &     2.35902305308496d+09,    -5.09939153709889d+08,
     &     8.40260941079563d+09,    -3.78918093051561d+10,
     &     9.20257216673170d+10,    -1.42925407961829d+11,
     &     9.49182293847050d+10,     2.62201763519032d+10,
     &    -3.23455414042411d+10,     2.70271437278597d+10,
     &    -1.04899384575720d+11,     1.54540018671944d+11,
     &    -3.52485992206419d+10,    -2.67958669856593d+11,
     &     3.04760172529868d+11,    -1.32939481060086d+11,
     &     1.04900112795205d+11,    -2.39990738389295d+11,
     &     4.56849545915060d+11,    -3.49217132942667d+11,
     &     9.52748491630969d+10,    -1.08896588644084d+11,
     &     2.34732673320838d+12,    -2.06147379952838d+13,
     &     1.13374034455934d+14,    -3.93167031841571d+14/
      END
! *********************************************************************
!                     Reproduce 2D Vint Ending
! *********************************************************************