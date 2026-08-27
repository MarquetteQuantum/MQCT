      subroutine jacobi_poly ( n, alpha, beta, x, cx )
! Taken from Numerical Recipes

c*********************************************************************72
c
cc JACOBI_POLY evaluates the Jacobi polynomials at X.
c
c  Differential equation:
c
c    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
c
c  Recursion:
c
c    P(0,ALPHA,BETA,X) = 1,
c
c    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
c
c    P(N,ALPHA,BETA,X)  = 
c      ( 
c        (2*N+ALPHA+BETA-1) 
c        * ((ALPHA^2-BETA^2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X) 
c        * P(N-1,ALPHA,BETA,X)
c        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
c      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
c
c  Restrictions:
c
c    -1 < ALPHA
c    -1 < BETA
c
c  Norm:
c
c    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA 
c      * P(N,ALPHA,BETA,X)^2 dX 
c    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
c      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
c
c  Special values:
c
c    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.  
c    Note that polynomials 0 through N will be computed.
c
c    Input, double precision ALPHA, one of the parameters defining the Jacobi
c    polynomials, ALPHA must be greater than -1.
c
c    Input, double precision BETA, the second parameter defining the Jacobi
c    polynomials, BETA must be greater than -1.
c
c    Input, double precision X, the point at which the polynomials are 
c    to be evaluated.
c
c    Output, double precision CX(0:N), the values of the first N+1 Jacobi
c    polynomials at the point X.
c
      implicit none

      integer n

      double precision alpha
      double precision beta
      double precision cx(0:n)
      double precision c1
      double precision c2
      double precision c3
      double precision c4
      integer i
      double precision r_i
      double precision x

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
        write ( *, '(a,g14.6)' ) 
     &    '  Illegal input value of ALPHA = ', alpha
        write ( *, '(a)' ) '  But ALPHA must be greater than -1.'
        stop
      end if
 
      if ( beta .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
        write ( *, '(a,g14.6)' ) 
     &    '  Illegal input value of BETA = ', beta
        write ( *, '(a)' ) '  But BETA must be greater than -1.'
        stop
      end if
  
      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      cx(1) = ( 1.0d00 + 0.5d00 * ( alpha + beta ) ) * x 
     &  + 0.5d00 * ( alpha - beta )
	 
       do i = 2, n

        r_i = dble ( i ) 

        c1 = 2.0D00 * r_i * ( r_i + alpha + beta ) 
     &    * ( 2.0D00 * r_i - 2.0D00 + alpha + beta )

        c2 = ( 2.0D00 * r_i - 1.0D00 + alpha + beta ) 
     &    * ( 2.0D00 * r_i  + alpha + beta ) 
     &    * ( 2.0D00 * r_i - 2.0D00 + alpha + beta )

        c3 = ( 2.0D00 * r_i - 1.0D00 + alpha + beta ) 
     &    * ( alpha + beta ) * ( alpha - beta )

        c4 = - 2.0D00 * ( r_i - 1.0D00 + alpha ) 
     &    * ( r_i - 1.0D00 + beta ) 
     &    * ( 2.0D00 * r_i + alpha + beta )

        cx(i) = ( ( c3 + c2 * x ) * cx(i-1) + c4 * cx(i-2) ) / c1

      end do

      return
      end      	  	  
  
      FUNCTION  sign_v(a)
      IMPLICIT NONE	  
      REAL*8 a,sign_v
      IF(a.ge.0d0) sign_v = 1d0
      IF(a.lt.0d0) sign_v = -1d0	  
      END FUNCTION sign_v	  
     
      SUBROUTINE  gauleg(ngp, xabsc, weig)
! Taken from Numerical Recipes
      IMPLICIT NONE	  
      INTEGER  i, j, m
      REAL*8  p1, p2, p3, pp, z, z1
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL*8, INTENT(OUT) :: xabsc(ngp), weig(ngp)
      REAL*8  :: EPS, M_PI
      PARAMETER (EPS=3.0d-15)       	!EPS is the relative precision
      PARAMETER (M_PI=dacos(-1d0))      ! Pi value


	   m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

	     do i = 1, m				! Loop over the desired roots */

     		z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON'S method   */
100     	p1 = 1.0d0
        	p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

        	do j = 1, ngp
           	p3 = p2
           	p2 = p1
           	p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
        	enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
        	pp = ngp*(z*p1-p2)/(z*z-1.0d0)
        	z1 = z
        	z = z1 - p1/pp             ! Newton's Method  */

        	if (dabs(z-z1) .gt. EPS) GOTO  100

      	xabsc(i) =  - z                    	! Roots will be bewteen -1.0 & 1.0 */
      	xabsc(ngp+1-i) =  + z                	! and symmetric about the origin  */
      	weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
      	weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

      end do     ! i loop

      END SUBROUTINE gauleg
	  
 
      FUNCTION plgndr(l,m,x)
! Taken from Numerical Recipes
      INTEGER l,m
      REAL*8 plgndr,x
      INTEGER i,ll
      REAL*8 fact,pll,pmm,pmmp1,somx2
c	  write(*,*) x
      if(m.lt.0. or. m.gt.l .or. abs(x).gt.1.) then
      write(*,*) m,l,abs(x)
      stop 'bad arguments in plgndr'
      endif
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END
	  
	  
      REAL*8 FUNCTION CG(J1,J2,J3,M1,M2,M3)
      INTEGER J1,J2,J3,M1,M2,M3
      REAL*8 W3JS
      EXTERNAL  W3JS
      CG = W3JS(J1*2,J2*2,J3*2,M1*2,M2*2,-M3*2)
     & *sqrt(2d0*J3+1d0)*(-1)**(-j1+j2-m3)
      END FUNCTION CG
	  
      REAL*8 FUNCTION CG_HF(J1,J2,J3,M1,M2,M3)
      REAL*8 J1,J3,M1,M3
	  INTEGER J2,M2,round
      INTEGER J1d,J2d,J3d,M1d,M2d,M3d	  
      REAL*8 W3JS
      EXTERNAL  W3JS,round
	  J1d = round(j1*2)
	  J2d = J2*2 !!!! EXPANSION
	  J3d = round(j3*2)
	  M1d = round(m1*2)
	  M2d = M2*2 !!! EXPANSION
	  M3d = round(m3*2)		  
      CG_HF = W3JS(J1d,J2d,J3d,M1d,M2d,-M3d)
     & *sqrt(2d0*J3+1d0)*(-1)**(-round(j1+j2-m3))!!!!! work on it
      END FUNCTION CG_HF	  
	   
       FUNCTION W3JS(J1,J2,J3,M1,M2,M3)
       USE FACTORIAL
       IMPLICIT NONE
       INTEGER Z,ZMIN,ZMAX
       REAL*8 W3JS,CC,CC1,CC2,DENOM
       INTEGER M1,M2,M3,J1,J2,J3
       INTEGER IA,IB,IC,JSUM,ID,IE,IF
       INTEGER IH,IG
       W3JS=0.0
       IF(M1+M2+M3.NE.0) GOTO 1000
       IA=J1+J2
       IF(J3.GT.IA) GOTO 1000
       IB=J1-J2
       IF(J3.LT.IABS(IB)) GOTO 1000
       JSUM=J3+IA
       IC=J1-M1
       ID=J2-M2
       IF(MOD(JSUM,2).NE.0) GOTO 1000
       IF(MOD(IC,2).NE.0) GOTO 1000
       IF(MOD(ID,2).NE.0) GOTO 1000
       IF(IABS(M1).GT.J1) GOTO 1000
       IF(IABS(M2).GT.J2) GOTO 1000
       IF(IABS(M3).GT.J3) GOTO 1000
       IE=J3-J2+M1
       IF=J3-J1-M2
       ZMIN=MAX0(0,-IE,-IF)
       IG=IA-J3
       IH=J2+M2
       ZMAX=MIN0(IG,IH,IC)
       CC=0.0
       DO 200 Z=ZMIN,ZMAX,2
       DENOM=EXP(LOGFACT(Z/2)+LOGFACT((IG-Z)/2)+LOGFACT((IC-Z)/2)
     1 +LOGFACT((IH-Z)/2)+LOGFACT((IE+Z)/2)+LOGFACT((IF+Z)/2))
      IF(MOD(Z,4).NE.0) DENOM=-DENOM
      CC=CC+1.0/DENOM
200   CONTINUE
      CC1=EXP(LOGFACT(IG/2)+LOGFACT((J3+IB)/2)+LOGFACT((J3-IB)/2)
     1 -LOGFACT((JSUM+2)/2))
      CC2=EXP(LOGFACT((J1+M1)/2)+LOGFACT(IC/2)+LOGFACT(IH/2)
     1 +LOGFACT(ID/2)+LOGFACT((J3-M3)/2)+LOGFACT((J3+M3)/2))
      CC=CC*SQRT(CC1*CC2)
      IF(MOD(IB-M3,4).NE.0) CC=-CC
      W3JS=CC
1000  RETURN
      END

! Bikram Start: Feb 2022
! This part is to use recursive method of computing Wigner 3j symbol

      REAL*8 FUNCTION CG_bikram(J1,J2,J3,M1,M2,M3,check_fact)
      INTEGER J1,J2,J3,M1,M2,M3
	  integer i, counter, jmin_recur, jmax_recur
      real*8 bk3j, w3j(J2+J3+1)
	  logical check_fact
      REAL*8 W3JS
      EXTERNAL  W3JS
	  
	  if(.not.check_fact) then
      call Wigner3j(w3j,jmin_recur,jmax_recur, J2,J3,M1,M2,-M3)

      counter = 0
      do i = jmin_recur, jmax_recur
      counter = counter + 1
      if(i.eq.J1) then
	  bk3j = w3j(counter)
      exit
	  end if
      end do
	  
      CG_bikram = bk3j*sqrt(2d0*J3+1d0)*(-1)**(-j1+j2-m3)
	  
	  else if(check_fact) then
      CG_bikram = W3JS(J1*2,J2*2,J3*2,M1*2,M2*2,-M3*2)
     & *sqrt(2d0*J3+1d0)*(-1)**(-j1+j2-m3)
	  
	  else
	  print*, "Something is wrong in computing Clebsh-Gordan COeffs"
	  stop
	  
	  end if
      END FUNCTION CG_bikram

      subroutine Wigner3j(w3j, jmin, jmax, j2, j3, m1, m2, m3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !	This subroutine will calculate the Wigner 3j symbols
    !
    !		j  j2 j3
    !		m1 m2 m3
    !
    !	for all allowable values of j. The returned values in the array j are 
    !	calculated only for the limits
    !
    !		jmin = max(|j2-j3|, |m1|)
    !		jmax = j2 + j3
    !
    !	To be non-zero, m1 + m2 + m3 = 0. In addition, it is assumed that all j and m are 
    !	integers. Returned values have a relative error less than ~1.d-8 when j2 and j3 
    !	are less than 103 (see below). In practice, this routine is probably usable up to 165.
    !
    !	This routine is based upon the stable non-linear recurence relations of Luscombe and 
    !	Luban (1998) for the "non classical" regions near jmin and jmax. For the classical 
    !	region, the standard three term recursion relationship is used (Schulten and Gordon 1975). 
    !	Note that this three term recursion can be unstable and can also lead to overflows. Thus 
    !	the values are rescaled by a factor "scalef" whenever the absolute value of the 3j coefficient 
    !	becomes greater than unity. Also, the direction of the iteration starts from low values of j
    !	to high values, but when abs(w3j(j+2)/w3j(j)) is less than one, the iteration will restart 
    !	from high to low values. More efficient algorithms might be found for specific cases 
    !	(for instance, when all m's are zero).
    !
    !	Verification: 
    !
    !	The results have been verified against this routine run in quadruple precision.
    !	For 1.e7 acceptable random values of j2, j3, m2, and m3 between -200 and 200, the relative error
    !	was calculated only for those 3j coefficients that had an absolute value greater than 
    !	1.d-17 (values smaller than this are for all practical purposed zero, and can be heavily 
    !	affected by machine roundoff errors or underflow). 853 combinations of parameters were found
    !	to have relative errors greater than 1.d-8. Here I list the minimum value of max(j2,j3) for
    !	different ranges of error, as well as the number of times this occured
    !	
    !	1.d-7 < error  <=1.d-8 = 103	# = 483
    !	1.d-6 < error <= 1.d-7 =  116	# = 240
    !	1.d-5 < error <= 1.d-6 =  165	# = 93
    !	1.d-4 < error <= 1.d-5 = 167	# = 36
    !
    !	Many times (maybe always), the large relative errors occur when the 3j coefficient 
    !	changes sign and is close to zero. (I.e., adjacent values are about 10.e7 times greater 
    !	in magnitude.) Thus, if one does not need to know highly accurate values of the 3j coefficients
    !	when they are almost zero (i.e., ~1.d-10) then this routine is probably usable up to about 160.
    !
    !	These results have also been verified for parameter values less than 100 using a code
    !	based on the algorith of de Blanc (1987), which was originally coded by Olav van Genabeek, 
    !	and modified by M. Fang (note that this code was run in quadruple precision, and
    !	only calculates one coefficient for each call. I also have no idea if this code
    !	was verified.) Maximum relative errors in this case were less than 1.d-8 for a large number
    !	of values (again, only 3j coefficients greater than 1.d-17 were considered here).
    !	
    !	The biggest improvement that could be made in this routine is to determine when one should
    !	stop iterating in the forward direction, and start iterating from high to low values. 
    !
    !	Calling parameters
    !		IN	
    !			j2, j3, m1, m2, m3 	Integer values.
    !		OUT	
    !			w3j			Array of length jmax - jmin + 1.
    !			jmin, jmax		Minimum and maximum values
    !						out output array.
    !	Dependencies: None
    !	
    !	Written by Mark Wieczorek August (2004)
    !
    !	August 2009: Based on the suggestions of Roelof Rietbroek, the calculation of RS has been slightly
    !	modified so that division by zero will not cause a run time crash (this behavior depends on how the 
    !	compiler treats IEEE floating point exceptions). These values were never used in the original code 
    !	when this did occur.
    !
    !	Copyright (c) 2005-2009, Mark A. Wieczorek
    !	All rights reserved.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		implicit none
		integer, intent(in) ::	j2, j3, m1, m2, m3
		integer, intent(out) ::	jmin, jmax
		real*8, intent(out) ::	w3j(j2+j3+1)
		real*8 :: wnmid, wpmid, scalef, denom, rs(j2+j3+1),
     & wl(j2+j3+1), wu(j2+j3+1), xjmin, yjmin, yjmax, zjmax, xj, zj
		integer :: j, jnum, jp, jn, k, flag1, flag2, jmid


		if (size(w3j) < j2+j3+1) then
		  print*, "Error --- Wigner3j"
		  print*, "W3J must be dimensioned (J2+J3+1) where J2 and
     & J3 are ", j2, j3
		  print*, "Input array is dimensioned ", size(w3j)
		  stop
		endif

		w3j = 0.0d0

		flag1 = 0
		flag2 = 0

		scalef = 1.0d3

		jmin = max(abs(j2-j3), abs(m1))
		jmax = j2 + j3
		jnum = jmax - jmin + 1

		if (abs(m2) > j2 .or. abs(m3) > j3) then
		  return
		elseif (m1 + m2 + m3 /= 0) then
		  return
		elseif (jmax < jmin) then
		  return
		endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !	Only one term is present
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		if (jnum == 1) then
		  w3j(1) = 1.0d0 / sqrt(2.0d0*jmin+1.0d0)
		  if ( (w3j(1) < 0.0d0 .and. (-1)**(j2-j3+m2+m3) > 0) .or.
     & (w3j(1) > 0.0d0 .and. (-1)**(j2-j3+m2+m3) < 0) ) 
     & w3j(1) = -w3j(1)
		  return	
		endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! 	Calculate lower non-classical values for [jmin, jn]. If the second term
    !	can not be calculated because the recursion relationsips give rise to a
    !	1/0, then set flag1 to 1.  If all m's are zero, then this is not a problem 
    !	as all odd terms must be zero.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		rs = 0.0d0
		wl = 0.0d0

		xjmin = x(jmin)
		yjmin = y(jmin)

		if (m1 == 0 .and. m2 == 0 .and. m3 == 0) then		! All m's are zero

		  wl(jindex(jmin)) = 1.0d0
		  wl(jindex(jmin+1)) = 0.0d0
		  jn = jmin+1

		elseif (yjmin == 0.0d0) then				! The second terms is either zero

		  if (xjmin == 0.0d0) then			! or undefined
			flag1 = 1
			jn = jmin
		  else
			wl(jindex(jmin)) = 1.0d0
			wl(jindex(jmin+1)) = 0.0d0
			jn = jmin+1
		  endif

		elseif ( xjmin * yjmin >= 0.0d0) then			! The second term is outside of the 
		  ! non-classical region 
		  wl(jindex(jmin)) = 1.0d0
		  wl(jindex(jmin+1)) = -yjmin / xjmin
		  jn = jmin+1

		else							! Calculate terms in the non-classical region

		  rs(jindex(jmin)) = -xjmin / yjmin

		  jn = jmax
		  do j=jmin + 1, jmax-1, 1
			denom =  y(j) + z(j)*rs(jindex(j-1))
			xj = x(j)
			if (abs(xj) > abs(denom) .or. xj * denom >= 0.0d0 .or. 
     & denom == 0.0d0) then
			  jn = j-1
			  exit
			else
			  rs(jindex(j)) = -xj / denom
			endif

		  enddo

		  wl(jindex(jn)) = 1.0d0

		  do k=1, jn - jmin, 1
			wl(jindex(jn-k)) = wl(jindex(jn-k+1)) * rs(jindex(jn-k))
		  enddo

		  if (jn == jmin) then					! Calculate at least two terms so that
			wl(jindex(jmin+1)) = -yjmin / xjmin		! these can be used in three term
			jn = jmin+1					! recursion

		  endif

		endif

		if (jn == jmax) then					! All terms are calculated

		  w3j(1:jnum) = wl(1:jnum)
		  call normw3j
		  call fixsign

		  return

		endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! 	Calculate upper non-classical values for [jp, jmax].
    !	If the second last term can not be calculated because the
    !	recursion relations give a 1/0, then set flag2 to 1. 
    !	(Note, I don't think that this ever happens).
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		wu = 0.0d0

		yjmax = y(jmax)
		zjmax = z(jmax)

		if (m1 == 0 .and. m2 == 0 .and. m3 == 0) then

		  wu(jindex(jmax)) = 1.0d0
		  wu(jindex(jmax-1)) = 0.0d0
		  jp = jmax-1

		elseif (yjmax == 0.0d0) then

		  if (zjmax == 0.0d0) then
			flag2 = 1
			jp = jmax
		  else
			wu(jindex(jmax)) = 1.0d0
			wu(jindex(jmax-1)) = - yjmax / zjmax
			jp = jmax-1
		  endif

		elseif (yjmax * zjmax >= 0.0d0) then

		  wu(jindex(jmax)) = 1.0d0
		  wu(jindex(jmax-1)) = - yjmax / zjmax
		  jp = jmax-1

		else
		  rs(jindex(jmax)) = -zjmax / yjmax

		  jp = jmin
		  do j=jmax-1, jn, -1
			denom = y(j) + x(j)*rs(jindex(j+1))
			zj = z(j)
			if (abs(zj) > abs(denom) .or. zj * denom >= 0.0d0 .or. 
     & denom == 0.0d0) then
			  jp = j+1
			  exit
			else
			  rs(jindex(j)) = -zj / denom
			endif
		  enddo	

		  wu(jindex(jp)) = 1.0d0

		  do k=1, jmax - jp, 1
			wu(jindex(jp+k)) = wu(jindex(jp+k-1))*rs(jindex(jp+k))
		  enddo

		  if (jp == jmax) then
			wu(jindex(jmax-1)) = - yjmax / zjmax
			jp = jmax-1
		  endif

		endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! 	Calculate classical terms for [jn+1, jp-1] using standard three
    ! 	term rercusion relationship. Start from both jn and jp and stop at the
    ! 	midpoint. If flag1 is set, then perform the recursion solely from high to
    ! 	low values. If flag2 is set, then perform the recursion solely from low to high.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		if (flag1 == 0) then

		  jmid = (jn + jp)/2

		  do j=jn, jmid - 1, 1			
			wl(jindex(j+1)) = - (z(j)*wl(jindex(j-1)) +y(j)
     & *wl(jindex(j))) / x(j)

			if (abs(wl(jindex(j+1))) > 1.0d0) then				! watch out for overflows.
			  wl(jindex(jmin):jindex(j+1)) = 
     & wl(jindex(jmin):jindex(j+1)) / scalef
			endif

			if (abs(wl(jindex(j+1)) / wl(jindex(j-1))) < 1.0d0 			! if values are decreasing
     & .and. wl(jindex(j+1)) /= 0.0d0) then				! then stop upward iteration
			  jmid = j+1						! and start with the downward
		  exit							! iteration.
			endif
		  enddo

        if (jmid <= jmin) jmid = jmin + 1
		  wnmid = wl(jindex(jmid))
		  if (jmid > jmin .and. abs(wnmid/wl(jindex(jmid-1))) < 1.d-6
     &   .and. wl(jindex(jmid-1)) /= 0.0d0) then
		  	wnmid = wl(jindex(jmid-1))
		  	jmid = jmid - 1
		  endif


		  do j=jp, jmid+1, -1
			wu(jindex(j-1)) = - (x(j)*wu(jindex(j+1)) + y(j)
     & *wu(jindex(j)) ) / z(j)
			if (abs(wu(jindex(j-1))) > 1.0d0) then
			  wu(jindex(j-1):jindex(jmax)) = 
     & wu(jindex(j-1):jindex(jmax)) / scalef
			endif

		  enddo

		  wpmid = wu(jindex(jmid))

      ! rescale two sequences to common midpoint

		  if (jmid == jmax) then
			w3j(1:jnum) = wl(1:jnum)
		  elseif (jmid == jmin) then
			w3j(1:jnum) = wu(1:jnum)
		  else
			w3j(1:jindex(jmid)) = wl(1:jindex(jmid)) * wpmid / wnmid 
			w3j(jindex(jmid+1):jindex(jmax)) = 
     & wu(jindex(jmid+1):jindex(jmax))
		  endif

		elseif (flag1 == 1 .and. flag2 == 0) then	! iterature in downward direction only

		  do j=jp, jmin+1, -1
			wu(jindex(j-1)) = - (x(j)*wu(jindex(j+1)) +
     & y(j)*wu(jindex(j)) ) / z(j)
			if (abs(wu(jindex(j-1))) > 1) then
			  wu(jindex(j-1):jindex(jmax)) = 
     & wu(jindex(j-1):jindex(jmax)) / scalef
			endif
		  enddo

		  w3j(1:jnum) = wu(1:jnum)

		elseif (flag2 == 1 .and. flag1 == 0) then	! iterature in upward direction only

		  do j=jn, jp-1, 1
			wl(jindex(j+1)) = - (z(j)*wl(jindex(j-1)) +
     & y(j)*wl(jindex(j))) / x(j)
			if (abs(wl(jindex(j+1))) > 1) then
			  wl(jindex(jmin):jindex(j+1)) = 
     & wl(jindex(jmin):jindex(j+1))/ scalef
			endif
		  enddo

		  w3j(1:jnum) = wl(1:jnum)

		elseif (flag1 == 1 .and. flag2 == 1) then

		  print*, "Fatal Error --- Wigner3j"
		  print*, "Can not calculate function for input values,
     & both flag1 and flag 2 are set."
		  stop
		endif


		call normw3j
		call fixsign


		  contains

			integer function jindex(j)
				integer :: j
				jindex = j-jmin+1
			end function jindex

			real*8 function a(j)
			  integer :: j
			  a = (dble(j)**2 - dble(j2-j3)**2) * (dble(j2+j3+1)**2 
     & - dble(j)**2) * (dble(j)**2-dble(m1)**2)
			  a = sqrt(a)
		  end function a

		  real*8 function y(j)
			integer :: j
			y = -dble(2*j+1) * 
     & ( dble(m1) * (dble(j2)*dble(j2+1) - 
     & dble(j3)*dble(j3+1) ) - dble(m3-m2)*dble(j)*dble(j+1) )
		  end function y

		real*8 function x(j)	
		  integer :: j
		  x = dble(j) * a(j+1)
			end function x

			real*8 function z(j)
			  integer :: j
			  z = dble(j+1)*a(j)
		  end function z

		  subroutine normw3j
			  real*8:: norm
			  integer j

			  norm = 0.0d0
			  do j = jmin, jmax
				norm = norm + dble(2*j+1) * w3j(jindex(j))**2
			  enddo

			  w3j(1:jnum) = w3j(1:jnum) / sqrt(norm)

		  end subroutine normw3j

		  subroutine fixsign

			  if ( (w3j(jindex(jmax)) < 0.0d0 .and. 
     & (-1)**(j2-j3+m2+m3) > 0) .or. 
     & (w3j(jindex(jmax)) > 0.0d0 .and. 
     & (-1)**(j2-j3+m2+m3) < 0) ) then
				w3j(1:jnum) = -w3j(1:jnum)
			  endif

		  end subroutine fixsign


		end subroutine Wigner3j

! Bikram End.

      FUNCTION gammln(xx)
      REAL*8 gammln,xx
!      Returns the value ln[Γ(xx)] for xx > 0.
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
!Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
!accuracy is good enough.
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     * 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     * -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do  j=1,6
      y=y+1.d0
      ser=ser+cof(j)/y
      enddo 
      gammln=tmp+log(stp*ser/x)
      return
      END	FUNCTION gammln  
	  
	  
      subroutine lf_function ( m, n, alpha, x, cx )

!*****************************************************************************80
!
!! LF_FUNCTION evaluates the Laguerre function Lf(n,alpha,x).
!
!  Recursion:
!
!    Lf(0,ALPHA,X) = 1
!    Lf(1,ALPHA,X) = 1+ALPHA-X
!
!    Lf(N,ALPHA,X) = (2*N-1+ALPHA-X)/N * Lf(N-1,ALPHA,X) 
!                      - (N-1+ALPHA)/N * Lf(N-2,ALPHA,X)
!
!  Restrictions:
!
!    -1 < ALPHA
!
!  Special values:
!
!    Lf(N,0,X) = L(N,X).
!    Lf(N,ALPHA,X) = LM(N,ALPHA,X) for ALPHA integral.
!
!  Norm:
!
!    Integral ( 0 <= X < +oo ) exp ( - X ) * Lf(N,ALPHA,X)^2 dX
!    = Gamma ( N + ALPHA + 1 ) / N!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order function to compute.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.  -1 < ALPHA is required.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(1:M,0:N), the functions of 
!    degrees 0 through N evaluated at the points X.
!
      implicit none

      integer ( kind = 4 ) m
      integer ( kind = 4 ) n

      real ( kind = 8 ) alpha
      real ( kind = 8 ) cx(1:m,0:n)
      integer ( kind = 4 ) i
      real ( kind = 8 ) x(1:m)

      if ( alpha <= -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LF_FUNCTION - Fatal error!'
        write ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
        write ( *, '(a)' ) '  but ALPHA must be greater than -1.'
        stop
      end if
 
      if ( n < 0 ) then
        return
      end if

      cx(1:m,0) = 1.0D+00

      if ( n == 0 ) then
        return
      end if

      cx(1:m,1) = 1.0D+00 + alpha - x(1:m)

      do i = 2, n
        cx(1:m,i) = (
     & ( real ( 2 * i - 1, kind = 8 ) + alpha - x(1:m) ) * cx(1:m,i-1)  
     &  + ( real ( - i + 1, kind = 8 ) - alpha)
     &	* cx(1:m,i-2) )
     & / real (  i,  kind = 8 )
      end do

      return
      end

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
! Taken from Numerical Recipes
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
! x1 < x2 < ... < xN, and given values yp1 and ypn for the first derivative of the interpolating
! function at points 1 and n, respectively, this routine returns an array y2(1:n) of
! length n which contains the second derivatives of the interpolating function at the tabulated
! points xi. If yp1 and/or ypn are equal to 1x10^30 or larger, the routine is signaled to set
! the corresponding boundary condition for a natural spline, with zero second derivative on
! that boundary. Parameter: NMAX is the largest anticipated value of n.
!
! Please notice that the program spline is called only once to process an entire tabulated
! function in arrays xi and yi. Once this has been done, the values of the interpolated function
! for any value of x are obtained by calls (as many as desired) to a separate routine splint.
!
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=600000)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
c      REAL*8, ALLOCATABLE :: u(:)
c      ALLOCATE(u(n))	  
      if (yp1.gt..99e30) then
      y2(1)=0.
      u(1)=0.
      else
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
  
      u(i)=(6.*((y(i+1)-y(i))/(x(i+
     * 1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     * u(i-1))/p
c      PRINT*, i, u(i),p,x(i+1),x(i)
11       continue
      if(ypn.gt..99e30) then
      qn=0.
      un=0.
      else
      qn=0.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
12       continue
c      DEALLOCATE(u)
      return
      END


      SUBROUTINE splint(xa,ya,y2a,n,x,y,y_prime)
! Taken from Numerical Recipes
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
! xai 's in order), and given the array y2a(1:n), which is the output from the "spline" subroutine,
! and given a value of x, this routine returns a cubic-spline interpolated value y at x.
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n),y_prime
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1       if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
      khi=k
      else
      klo=k
      endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
      WRITE(*,*) xa(n),x,n
      WRITE(*,*) "DISTANCE",xa	  
      STOP 'bad xa input in splint'
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     * 2)/6.
      y_prime=-ya(klo)/h+ya(khi)/h+(-(3*a**2-1d0)*y2a(klo)/h
     & +(3*b**2-1d0)*y2a(khi)/h) 
     * *(h**2)/6.0
      return
      END
	  

      SUBROUTINE splint_point(xa,ya,y2a,n,x,y,y_prime,n_p)
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
! xai 's in order), and given the array y2a(1:n), which is the output from the "spline" subroutine,
! and given a value of x, this routine returns a cubic-spline interpolated value y at x.
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n),y_prime
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=n_p
      khi=n_p+1
1       if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
      khi=k
      else
      klo=k
      endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
      WRITE(*,*) xa(n),x,n
      STOP 'bad xa input in splint'
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     * 2)/6.
      y_prime=-ya(klo)/h+ya(khi)/h+(-(3*a**2-1d0)*y2a(klo)/h
     & +(3*b**2-1d0)*y2a(khi)/h) 
     * *(h**2)/6.0
      return
      END	  

      FUNCTION rand0(idum)  ! ran2 Numerical Recipes
      implicit none
      INTEGER*4 idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      ReaL(8) rand0,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      rand0=min(AM*iy,RNMX)
   21 format(f17.12,100i11)
      return
      END function rand0
c*! SUUPOR ROUTIES TAKEN FROM NUMERICAL REICPIES       
	   
       SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
! Taken from Numerical Recipes
       implicit none
       INTEGER n
       REAL*8 h,x,dydx(n),y(n),yout(n)  
       EXTERNAL derivs  
       INTEGER i  
       REAL*8 h6,hh,xh
       REAL*8, ALLOCATABLE :: dym(:), dyt(:), yt(:)
 
       IF (n .LT. 1) THEN
         WRITE(*,*) 'rk4: invalid system size n=', n
         STOP 'rk4 invalid n'
       ENDIF
       ALLOCATE(dym(n), dyt(n), yt(n))

       hh=h*0.5d0  
       h6=h/6.d0  
       xh=x+hh  
 
       do 11 i=1,n  
         yt(i)=y(i)+hh*dydx(i)  
 11    continue  
 
       call derivs(xh,yt,dyt)  
 
       do 12 i=1,n  
         yt(i)=y(i)+hh*dyt(i)  
 12    continue  
 
       call derivs(xh,yt,dym)  
 
       do 13 i=1,n  
         yt(i)=y(i)+h*dym(i)  
         dym(i)=dyt(i)+dym(i)  
 13    continue  
 
       call derivs(x+h,yt,dyt)  
 
       do 14 i=1,n  
         yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))  
 14    continue  
       DEALLOCATE(dym, dyt, yt)
       return  
       END  


c-------------------------------------------------------------------
      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL*8 
     * eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),SAFE1,SAFE2,
     * REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=50000
     * ,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25d0,SAFE2=.7d0,
     * REDMAX=1.d-5,REDMIN=.7d0,TINY=1.d-30,SCALMX=.1d0)
CU    USES derivs,mmid,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL*8 eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,xnew,
     * a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),ysav(NMAX),
     * yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      EXTERNAL derivs
      DATA first/.true./,epsold/-1./
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x)pause 'step size underflow in bsstep'
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1./(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1./err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.e35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END

       SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs,
     & hmin)
! Taken from Numerical Recipes
       USE VARIABLES
       implicit none
       INTEGER n
       REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n),hmin  
       EXTERNAL derivs  
CU     USES derivs,rkck  
       INTEGER i  
       REAL*8 errmax,h,htemp,xnew,SAFETY,PGROW,  
     & PSHRNK,ERRCON  
       REAL*8, ALLOCATABLE :: yerr(:), ytemp(:)
       PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,
     & ERRCON=1.89d-4)  

       IF (n .LT. 1) THEN
         WRITE(*,*) 'rkqs: invalid system size n=', n
         STOP 'rkqs invalid n'
       ENDIF
       ALLOCATE(yerr(n), ytemp(n))

       h=htry  
 
 1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)  
	   if(h.le.hmin) then
		 hdid=h  
         x=x+h  
         do 43 i=1,n  
           y(i)=ytemp(i)  
 43      continue  
 
         DEALLOCATE(yerr, ytemp)
         return  
	   endif
       errmax=0.d0  
 
       IF (INCL_AMPTDS) THEN
         do 998 i=1,n  
           errmax=max(errmax,abs(yerr(i)/yscal(i)))  
 998     continue
       ELSE
         do 999 i=n-7,n-2  
           errmax=max(errmax,abs(yerr(i)/yscal(i)))  
 999     continue  
       ENDIF

       errmax=errmax/eps  

       if(errmax.gt.1.d0)then  
         htemp=SAFETY*h*(errmax**PSHRNK)  
         h=sign(max(abs(htemp),0.1d0*abs(h)),h)  
         xnew=x+h  
*         if(xnew.eq.x) PRINT*, 'stepsize underflow in rkqs'  
         goto 1  
       else  
         if(errmax.gt.ERRCON)then  
           hnext=SAFETY*h*(errmax**PGROW)  
         else  
           hnext=5.d0*h  
         endif  
 	 	 hdid=h  
         x=x+h  
         do 12 i=1,n  
           y(i)=ytemp(i)  
 12      continue  
 
         DEALLOCATE(yerr, ytemp)
         return  
       endif  
       END  


c------------------------------------------------------------------


       SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
! Taken from Numerical Recipes
       implicit none
       INTEGER n
       REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)  
       EXTERNAL derivs  
CU     USES derivs  
       INTEGER i  
       REAL*8 A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     & B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6  
       REAL*8, ALLOCATABLE :: ak2(:), ak3(:), ak4(:), ak5(:), ak6(:),
     & ytemp(:)
       PARAMETER (A2=.2d0,A3=.3d0,A4=.6,A5=1.d0,
     & A6=.875d0,B21=.2d0,B31=3.d0/40.d0,  
     & B32=9.d0/40.d0,B41=.3d0,
     & B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52=2.5d0,  
     & B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,
     & B62=175.d0/512.d0,  
     & B63=575.d0/13824.d0,B64=44275.d0/110592.d0,
     & B65=253.d0/4096.d0,C1=37.d0/378.d0,  
     & C3=250.d0/621.d0,C4=125.d0/594.d0,
     & C6=512.d0/1771.d0,DC1=C1-2825.d0/27648.d0,  
     & DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,
     & DC5=-277.d0/14336.d0,  
     & DC6=C6-.25d0)  

       IF (n .LT. 1) THEN
         WRITE(*,*) 'rkck: invalid system size n=', n
         STOP 'rkck invalid n'
       ENDIF
       ALLOCATE(ak2(n), ak3(n), ak4(n), ak5(n), ak6(n), ytemp(n))
 
       do 11 i=1,n  
         ytemp(i)=y(i)+B21*h*dydx(i)  
 11    continue  
 
       call derivs(x+A2*h,ytemp,ak2)  
 
       do 12 i=1,n  
         ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))  
 12    continue  
 
       call derivs(x+A3*h,ytemp,ak3)  
 
       do 13 i=1,n  
         ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))  
 13    continue  
 
       call derivs(x+A4*h,ytemp,ak4)  
 
       do 14 i=1,n  
         ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))  
 14    continue  
 
       call derivs(x+A5*h,ytemp,ak5)  
 
       do 15 i=1,n  
         ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+  
     & B65*ak5(i))  
 15    continue  
 
       call derivs(x+A6*h,ytemp,ak6)  
 
       do 16 i=1,n  
         yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))  
 16    continue  
 
       do 17 i=1,n  
         yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*  
     & ak6(i))  
 
 17    continue  
       DEALLOCATE(ak2, ak3, ak4, ak5, ak6, ytemp)
       return  
       END  


c-------------------------------------------------------------------


      subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *	rkqs,hnext,fail_odeint)
! Taken from Numerical Recipes and modified by Bikramaditya Mandal
      implicit none
      logical fail_odeint,last_step,end_prop
      integer nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      double precision eps,h1,hmin,x1,x2,ystart(nvar),TINY
      external derivs,rkqs,bsstep
      parameter(NMAX=100000,KMAXX=200,TINY=1.d-5)
      integer i,kmax,kount,nstp
      double precision dxsav,h,hdid,hnext,x,xsav
      double precision, allocatable :: dydx(:), y(:), yscal(:)
      double precision xp(KMAXX),yp(NMAX,KMAXX)
      double precision hk4	  
      common /path/ kmax,kount,dxsav,xp,yp

      IF (nvar .LT. 1) THEN
        WRITE(*,*) 'odeint: invalid nvar=', nvar
        STOP 'odeint invalid nvar'
      ENDIF
      ALLOCATE(dydx(nvar), y(nvar), yscal(nvar))

      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      kmax = 0	! MY OWN SETUP
!      hmin = max(hmin,abs(x2-x1)/(MAXSTP-1d0)) 
      MAXSTP  = int(abs(x2-x1)/hmin)  +1   	  
	  MAXSTP=10000				!Bikram
c	  do
      do 11 i=1,nvar
      y(i)=ystart(i)
11      continue
c		enddo
      if (kmax.gt.0) xsav=x-2.d0*dxsav
      fail_odeint = .FALSE.
      last_step = .FALSE.	  
c	  do
      do 16 nstp=1,MAXSTP
      CALL TIME_CHECKER_ODEINT(end_prop)
      IF(end_prop) THEN
      fail_odeint = .TRUE.	  
	  RETURN
      ENDIF	  
      call derivs(x,y,dydx)
      if(nstp.eq.MAXSTP) last_step = .TRUE.
      if(last_step) goto 345	  
c	  do
      do 12 i=1,nvar
      yscal(i)=dabs(y(i))+dabs(h*dydx(i))+TINY
12      continue
c		enddo
      if(kmax.gt.0)then
      if(abs(x-xsav).gt.dabs(dxsav))then
      if(kount.lt.kmax-1)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 13 i=1,nvar
      yp(i,kount)=y(i)
13      continue
c		enddo
      xsav=x
      endif
      endif
      endif
      if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x
	  
      call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,hmin)

      if(hdid.eq.h)then
      nok=nok+1
      else
      nbad=nbad+1
      endif
345   if(last_step) then
      hk4  = sign(abs(x2-x),x2-x1)
      call rk4(y,dydx,nvar,x,hk4,y,derivs)
      x = x2+sign(TINY,x2-x1)	  
      endif	  
      if((x-x2)*(x2-x1).ge.0.d0)then
c	  do
      do 14 i=1,nvar
      ystart(i)=y(i)
14      continue
c		enddo
      if(kmax.ne.0)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 15 i=1,nvar
      yp(i,kount)=y(i)
15      continue
c		enddo
      endif

      return
      endif
      if(abs(hnext).lt.hmin) then
      hnext = sign(hmin,x2-x1)	  
      endif	  
      h=hnext
16      continue
c		enddo
      fail_odeint = .TRUE.
      return
      end

*-----------------------------------------------------------------------

! Bikram Start Oct 2019:
      subroutine bk_int(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *	rkqs,hnext,fail_odeint,f_time,cut_r,nmbr_r,
     &  vibration_cnt,numb_oscl_prds,
     & nmbr_phi,period_cnt,numb_orb_prds)
! This subroutine is written by Bikramaditya Mandal
      implicit none
      logical fail_odeint,last_step,end_prop
      integer nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      double precision eps,h1,hmin,x1,x2,ystart(nvar),TINY
      external derivs,rkqs,bsstep
      parameter(NMAX=100000,KMAXX=200,TINY=1.d-5)
      integer i,kmax,kount,nstp	
	  integer nmbr_r,vibration_cnt,numb_oscl_prds							!Bikram
	  integer nmbr_phi,period_cnt,numb_orb_prds								!Bikram
      double precision dxsav,h,hdid,hnext,x,xsav
      double precision, allocatable :: dydx(:), y(:), yscal(:)
      double precision f_time,cut_r,tmp_R2,tmp_R1							!Bikram
      double precision xp(KMAXX),yp(NMAX,KMAXX)
      double precision hk4
	  real*8, parameter :: pi=4.0d0*datan(1.0d0)
      common /path/ kmax,kount,dxsav,xp,yp

      IF (nvar .LT. 1) THEN
        WRITE(*,*) 'bk_int: invalid nvar=', nvar
        STOP 'bk_int invalid nvar'
      ENDIF
      ALLOCATE(dydx(nvar), y(nvar), yscal(nvar))
	  
! Bikram Start:	  
	  tmp_R2=0.d0
	  tmp_R1=0.d0
	  vibration_cnt=0
	  period_cnt=0
! Bikram End.
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      kmax = 0	! MY OWN SETUP
!      MAXSTP  = int(abs(x2-x1)/hmin)  +1   	  
	  MAXSTP=10000
c	  do
      do 11 i=1,nvar
      y(i)=ystart(i)
11      continue
c		enddo
      if (kmax.gt.0) xsav=x-2.d0*dxsav
      fail_odeint = .FALSE.
      last_step = .FALSE.	  
c	  do
      do 16 nstp=1,MAXSTP
      CALL TIME_CHECKER_ODEINT(end_prop)
      IF(end_prop) THEN
      fail_odeint = .TRUE.	  
	  RETURN
      ENDIF	  
      call derivs(x,y,dydx)
      if(nstp.eq.MAXSTP) last_step = .TRUE.
      if(last_step) goto 345	  
c	  do
      do 12 i=1,nvar
      yscal(i)=dabs(y(i))+dabs(h*dydx(i))+TINY
12      continue
c		enddo
      if(kmax.gt.0)then
      if(abs(x-xsav).gt.dabs(dxsav))then
      if(kount.lt.kmax-1)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 13 i=1,nvar
      yp(i,kount)=y(i)
13      continue
c		enddo
      xsav=x
      endif
      endif
      endif
      if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x
      call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,hmin)
	  
! Bikram start Oct 2019:
	  if(y(nmbr_r).gt.cut_r) then 
	  f_time=x
c	  do
      do 114 i=1,nvar
      ystart(i)=y(i)
114      continue
c		enddo
      if(kmax.ne.0)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 115 i=1,nvar
      yp(i,kount)=y(i)
115      continue
c		enddo
      endif
	  return
	  endif
 
! Check if the molecule is trapped and vibrating:
      if(nstp.gt.2) then
	  if((tmp_R2-tmp_R1).lt.0.00d0 .and. (tmp_R1-y(nmbr_r)).gt.0.00d0) 
     & vibration_cnt=vibration_cnt + 1
	  endif
	  tmp_R2=tmp_R1
      tmp_R1=y(nmbr_r)
	  IF(vibration_cnt.ge.numb_oscl_prds) THEN
	  f_time=x
c	  do
      do 116 i=1,nvar
      ystart(i)=y(i)
116      continue
c		enddo
      if(kmax.ne.0)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 117 i=1,nvar
      yp(i,kount)=y(i)
117      continue
c		enddo
      endif
	  return
      ENDIF
	  
! Check if molecule is orbiting:
      period_cnt = int(abs(y(nmbr_phi))/2d0/pi)
      IF(numb_orb_prds.le.period_cnt) THEN
	  f_time=x
c	  do
      do 118 i=1,nvar
      ystart(i)=y(i)
118      continue
c		enddo
      if(kmax.ne.0)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 119 i=1,nvar
      yp(i,kount)=y(i)
119      continue
c		enddo
      endif
	  return
      ENDIF 
! Bikram End.
	  
! Bikram End

      if(hdid.eq.h)then
      nok=nok+1
      else
      nbad=nbad+1
      endif
345   if(last_step) then
      hk4  = sign(abs(x2-x),x2-x1)
      call rk4(y,dydx,nvar,x,hk4,y,derivs)
      x = x2+sign(TINY,x2-x1)	  
      endif	  
      if((x-x2)*(x2-x1).ge.0.d0)then
c	  do
      do 14 i=1,nvar
      ystart(i)=y(i)
14      continue
c		enddo
      if(kmax.ne.0)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 15 i=1,nvar
      yp(i,kount)=y(i)
15      continue
c		enddo
      endif

      return
      endif
      if(abs(hnext).lt.hmin) then
      hnext = sign(hmin,x2-x1)	  
      endif	  
      h=hnext
16      continue
c		enddo
      fail_odeint = .TRUE.
      return
      end

*-----------------------------------------------------------------------
! Bikram End.

! Bikram Start May 2020:
      subroutine bk_int_adia(ystart, nvar, x1, x2, eps, h1, hmin,
     &  nok, nbad, derivs, rkqs, hnext, fail_odeint, f_time, cut_r,
     & vibration_cnt, numb_oscl_prds,
     & period_cnt, numb_orb_prds,
     & bkt, bkn, bkr, bkdr, bkp, bkdp)
! This subroutine is written by Bikramaditya Mandal
      implicit none
      logical fail_odeint,last_step,end_prop
      integer nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      double precision eps,h1,hmin,x1,x2,ystart(nvar),TINY
      external derivs,rkqs,bsstep
      parameter(NMAX=100000,KMAXX=200,TINY=1.d-5)
      integer i,kmax,kount,nstp	
	  integer vibration_cnt,numb_oscl_prds							!Bikram
	  integer period_cnt,numb_orb_prds								!Bikram
      double precision dxsav,h,hdid,hnext,x,xsav
      double precision, allocatable :: dydx(:), y(:), yscal(:)
      double precision f_time,cut_r,tmp_R2,tmp_R1							!Bikram
      double precision xp(KMAXX),yp(NMAX,KMAXX)
      double precision hk4
	  real*8, parameter :: pi=4.0d0*datan(1.0d0)
	  integer bkn, fu
	  real*8 bkt(bkn), bkr(bkn), bkdr(bkn), bkp(bkn), bkdp(bkn)
	  real*8 rqr,rqp,yprr,yprp, ratmin, ratmax, ratm
	  CHARACTER(LEN=100) bkpr
      common /path/ kmax,kount,dxsav,xp,yp

      IF (nvar .LT. 1) THEN
        WRITE(*,*) 'bk_int_adia: invalid nvar=', nvar
        STOP 'bk_int_adia invalid nvar'
      ENDIF
      ALLOCATE(dydx(nvar), y(nvar), yscal(nvar))
	  
! Bikram Start:	
	  ratmin = h1
	  ratmax = h1
	  tmp_R2=0.d0
	  tmp_R1=0.d0
	  vibration_cnt=0
	  period_cnt=0
! Bikram End.
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      kmax = 0	! MY OWN SETUP
!      MAXSTP  = int(abs(x2-x1)/hmin)  +1   	  
	  MAXSTP=10000  
c	  do
      do 11 i=1,nvar
      y(i)=ystart(i)
11      continue
c		enddo
      if (kmax.gt.0) xsav=x-2.d0*dxsav
      fail_odeint = .FALSE.
      last_step = .FALSE.	  
c	  do
      do 16 nstp=1,MAXSTP
      CALL TIME_CHECKER_ODEINT(end_prop)
      IF(end_prop) THEN
      fail_odeint = .TRUE.	  
	  RETURN
      ENDIF	  
      call derivs(x,y,dydx)
      if(nstp.eq.MAXSTP) last_step = .TRUE.
      if(last_step) goto 345	  
c	  do
      do 12 i=1,nvar
      yscal(i)=dabs(y(i))+dabs(h*dydx(i))+TINY
12      continue
c		enddo
      if(kmax.gt.0)then
      if(abs(x-xsav).gt.dabs(dxsav))then
      if(kount.lt.kmax-1)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 13 i=1,nvar
      yp(i,kount)=y(i)
13      continue
c		enddo
      xsav=x
      endif
      endif
      endif
      if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x
	  if(hnext.lt.ratmin) ratmin = hnext
	  if(hnext.gt.ratmax) ratmax = hnext
	  
      call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,hmin)
	  
! Bikram start Oct 2019:
	  call splint(bkt, bkr, bkdr, bkn, x, rqr, yprr)
	  call splint(bkt, bkp, bkdp, bkn, x, rqp, yprp)
	  
	  if(rqr.gt.cut_r) then 
	  f_time=x
c	  do
      do 114 i=1,nvar
      ystart(i)=y(i)
114      continue
c		enddo
      if(kmax.ne.0)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 115 i=1,nvar
      yp(i,kount)=y(i)
115      continue
c		enddo
      endif
	  ratm = abs(ratmax)/abs(ratmin)
	  return
	  endif
 
! Check if the molecule is trapped and vibrating:
      if(nstp.gt.2) then
	  if((tmp_R2-tmp_R1).lt.0.00d0 .and. (tmp_R1-rqr).gt.0.00d0) 
     & vibration_cnt=vibration_cnt + 1
	  endif
	  tmp_R2=tmp_R1
      tmp_R1=rqr
	  IF(vibration_cnt.ge.numb_oscl_prds) THEN
	  f_time=x
c	  do
      do 116 i=1,nvar
      ystart(i)=y(i)
116      continue
c		enddo
      if(kmax.ne.0)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 117 i=1,nvar
      yp(i,kount)=y(i)
117      continue
c		enddo
      endif
	  ratm = abs(ratmax)/abs(ratmin)
	  return
      ENDIF
	  
! Check if molecule is orbiting:
      period_cnt = int(abs(rqp)/2d0/pi)
      IF(numb_orb_prds.le.period_cnt) THEN
	  f_time=x
c	  do
      do 118 i=1,nvar
      ystart(i)=y(i)
118      continue
c		enddo
      if(kmax.ne.0)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 119 i=1,nvar
      yp(i,kount)=y(i)
119      continue
c		enddo
      endif
	  ratm = abs(ratmax)/abs(ratmin)
	  return
      ENDIF 
! Bikram End.
	  
! Bikram End

      if(hdid.eq.h)then
      nok=nok+1
      else
      nbad=nbad+1
      endif
345   if(last_step) then
      hk4  = sign(abs(x2-x),x2-x1)
      call rk4(y,dydx,nvar,x,hk4,y,derivs)
      x = x2+sign(TINY,x2-x1)	  
      endif	  
      if((x-x2)*(x2-x1).ge.0.d0)then
c	  do
      do 14 i=1,nvar
      ystart(i)=y(i)
14      continue
c		enddo
      if(kmax.ne.0)then
      kount=kount+1
      xp(kount)=x
c	  do
      do 15 i=1,nvar
      yp(i,kount)=y(i)
15      continue
c		enddo
      endif

      return
      endif
      if(abs(hnext).lt.hmin) then
      hnext = sign(hmin,x2-x1)	  
      endif	  
      h=hnext
16      continue
c		enddo
      fail_odeint = .TRUE.
      return
      end

*-----------------------------------------------------------------------
! Bikram End.

      subroutine rk2(y0,n,x,h,y2,y3,derivs)
! Taken from Numerical Recipes
      implicit none
      integer n
      real*8 y0(n),y2(n),y3(n),x,h
      external derivs
      call derivs(x,y0,y3)
      y2=y0+(h/2)*y3
      call derivs(x+h/2,y2,y3)
      y0=y0+y3*h
      end
*-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      function scalpr(x,y)
      implicit none
      real*8 x(3),y(3),scalpr
      scalpr=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
      end function scalpr
      subroutine h_polynomial_value ( m, n, x, p )

!*****************************************************************************80
!
!! H_POLYNOMIAL_VALUE evaluates H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!  Differential equation:
!
!    Y'' - 2 X Y' + 2 N Y = 0
!
!  First terms:
!
!      1
!      2 X
!      4 X^2     -  2
!      8 X^3     - 12 X
!     16 X^4     - 48 X^2     + 12
!     32 X^5    - 160 X^3    + 120 X
!     64 X^6    - 480 X^4    + 720 X^2    - 120
!    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
!    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
!    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
!   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Norm:
!
!    Integral ( -oo < X < oo ) exp ( - X^2 ) * H(N,X)^2 dX
!    = sqrt ( PI ) * 2^N * N!
!
!    H(N,X) = (-1)^N * exp ( X^2 ) * dn/dXn ( exp(-X^2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Larry Andrews,
!    Special Functions of Mathematics for Engineers,
!    Second Edition, 
!    Oxford University Press, 1998.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) P(M,0:N), the values of the first N+1 Hermite
!    polynomials at the point X.
!
      implicit none

      integer ( kind = 4 ) m
      integer ( kind = 4 ) n

      integer ( kind = 4 ) j
      real ( kind = 8 ) p(m,0:n)
      real ( kind = 8 ) x(m)

      if ( n < 0 ) then
      return
      end if

      p(1:m,0) = 1.0D+00

      if ( n == 0 ) then
      return
      end if

      p(1:m,1) = 2.0D+00 * x(1:m)
 
      do j = 2, n
      p(1:m,j) = 2.0D+00 * x(1:m) * p(1:m,j-1) 
     & - 2.0D+00 * real ( j - 1, kind = 8 ) * p(1:m,j-2)
      end do
 
      return
      end

      subroutine r8_fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, s )

!*****************************************************************************80
!
!! R8_FEHL takes one Fehlberg fourth-fifth order step (double precision).
!
!  Discussion:
!
!    This routine integrates a system of NEQN first order ordinary differential
!    equations of the form
!      dY(i)/dT = F(T,Y(1:NEQN))
!    where the initial values Y and the initial derivatives
!    YP are specified at the starting point T.
!
!    The routine advances the solution over the fixed step H and returns
!    the fifth order (sixth order accurate locally) solution
!    approximation at T+H in array S.
!
!    The formulas have been grouped to control loss of significance.
!    The routine should be called with an H not smaller than 13 units of
!    roundoff in T so that the various independent arguments can be
!    distinguished.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Erwin Fehlberg,
!    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
!    NASA Technical Report R-315, 1969.
!
!    Lawrence Shampine, Herman Watts, S Davenport,
!    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
!    SIAM Review,
!    Volume 18, pages 376-411, 1976.
!
!  Parameters:
!
!    Input, external F, a user-supplied subroutine to evaluate the
!    derivatives Y'(T), of the form:
!
!      subroutine f ( t, y, yp )
!      real ( kind = 8 ) t
!      real ( kind = 8 ) y(*)
!      real ( kind = 8 ) yp(*)
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations to be integrated.
!
!    Input, real ( kind = 8 ) Y(NEQN), the current value of the
!    dependent variable.
!
!    Input, real ( kind = 8 ) T, the current value of the independent
!    variable.
!
!    Input, real ( kind = 8 ) H, the step size to take.
!
!    Input, real ( kind = 8 ) YP(NEQN), the current value of the
!    derivative of the dependent variable.
!
!    Output, real ( kind = 8 ) F1(NEQN), F2(NEQN), F3(NEQN), F4(NEQN),
!    F5(NEQN), derivative values needed for the computation.
!
!    Output, real ( kind = 8 ) S(NEQN), the estimate of the solution at T+H.
!
      implicit none

      integer ( kind = 4 ) neqn

      real ( kind = 8 ) ch
      external f
      real ( kind = 8 ) f1(neqn)
      real ( kind = 8 ) f2(neqn)
      real ( kind = 8 ) f3(neqn)
      real ( kind = 8 ) f4(neqn)
      real ( kind = 8 ) f5(neqn)
      real ( kind = 8 ) h
      real ( kind = 8 ) s(neqn)
      real ( kind = 8 ) t
      real ( kind = 8 ) y(neqn)
      real ( kind = 8 ) yp(neqn)

      ch = h / 4.0D+00

      f5(1:neqn) = y(1:neqn) + ch * yp(1:neqn)

      call f ( t + ch, f5, f1 )

      ch = 3.0D+00 * h / 32.0D+00

      f5(1:neqn) = y(1:neqn)+ch*(yp(1:neqn)+3.0D+00*f1(1:neqn))

      call f(t+3.0D+00*h/8.0D+00,f5,f2)

      ch = h/2197.0D+00

      f5(1:neqn) = y(1:neqn) + ch* 
     & ( 1932.0D+00*yp(1:neqn) 
     & + ( 7296.0D+00*f2(1:neqn)-7200.0D+00*f1(1:neqn))
     & )

      call f(t + 12.0D+00*h/13.0D+00,f5,f3)

      ch = h/4104.0D+00

      f5(1:neqn) = y(1:neqn) + ch* 
     & ( 
     & (8341.0D+00*yp(1:neqn) - 845.0D+00*f3(1:neqn)) 
     & + (29440.0D+00*f2(1:neqn) - 32832.0D+00*f1(1:neqn) ) 
     & )

      call f(t+h,f5,f4 )

      ch = h/20520.0D+00

      f1(1:neqn) = y(1:neqn) + ch*
     & ( 
     &   ( -6080.0D+00 * yp(1:neqn) 
     &   + ( 9295.0D+00 * f3(1:neqn) - 5643.0D+00 * f4(1:neqn) ) 
     &   ) 
     & + ( 41040.0D+00 * f1(1:neqn) - 28352.0D+00 * f2(1:neqn) ) 
     & )

      call f( t + h/2.0D+00, f1, f5 )
!
!  Ready to compute the approximate solution at T+H.
!
      ch = h/7618050.0D+00

      s(1:neqn) = y(1:neqn) + ch*
     & ( 
     & ( 902880.0D+00 * yp(1:neqn) 
     &  + ( 3855735.0D+00 * f3(1:neqn) - 1371249.0D+00 * f4(1:neqn) ) ) 
     & + ( 3953664.0D+00 * f2(1:neqn) + 277020.0D+00 * f5(1:neqn) ) 
     & )

      return
      end
	  
      subroutine vect_product(a,b,c)
      REAL*8 a(3), b(3), c(3)
      a(1) = b(2)*c(3) - b(3)*c(2)
      a(2) = b(3)*c(1) - b(1)*c(3)
      a(3) = b(1)*c(2) - b(2)*c(1) 	  
      end subroutine vect_product	  
      
      SUBROUTINE ERROR_SIGNALING(error_send_recv,part_of_input)
! This subroutine is updated by Bikramaditya Mandal
      USE ERRORS
      IMPLICIT NONE
      INTEGER error_send_recv,part_of_input
      SELECT CASE(part_of_input)
      CASE(1)
!!!!  BASIS	  
      SELECT CASE(error_send_recv)
      CASE(2)	  
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "LEVELS_FILE="	  
      CASE(3)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "NMB_OF_CHNLS="
      CASE(8)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "BE="
      CASE(45)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "IDENT="
      CASE(52)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "CS="
      CASE(53)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "PRINT_STATES="
      CASE(54)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "EXCLUDE_STATES="	  
!      CASE(77)
!      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
!     & " NN_APPROX="	
      END SELECT
 !!! END BASIS
 !!! SYSTEM STARTS
      CASE(2)
      SELECT CASE(error_send_recv)
      CASE(5)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "PROPAGATOR="
      CASE(14)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "DFLC_FNCT="	  
      CASE(20)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "SHOW_ORBIT="		  
      CASE(25)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "START_FROM_CHECK_FILE="
      CASE(27)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "MONTE_CARLO="
      CASE(28)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "CALC_ELAST="
! Bikram Start Nov 2019:
	  CASE(37)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "BIKRAM_INT="
! Bikram End.
      END SELECT
!!! END SYSTEM
!!! POTNETIAL
      CASE(3)
      SELECT CASE(error_send_recv)
      CASE(3)	  
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "R_UNITS="
      CASE(11)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "MTRX_READ="
      CASE(12)	 
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "PROG_RUN="
      CASE(13)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "ATOMIC_DIST_COORD="
      CASE(14)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "MAKE_GRID_FILE="
      CASE(15)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "CALC_EXPANSION="	  
      CASE(22)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "R_GRID_FROM_FILE="	  
      CASE(27)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "PRINT_MIJ="	  
      CASE(28)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "TERMS_ONFLY="	  
      CASE(29)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "UNFORMAT="
      CASE(30)
      PRINT*, "CHECK INPUT VARIABLE FOR THE KEY WORD",
     & "PRINT_ELAST="	 
      END SELECT	  
      END SELECT	  
    	  
      STOP	  
      END SUBROUTINE ERROR_SIGNALING

      SUBROUTINE TIME_CHECKER_ODEINT(RETURN_VAL)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE CHECK_FILE
      USE MPI_TASK_TRAJECT	  
      USE INTEGRATOR	  
      USE MPI
      USE MPI_DATA
      IMPLICIT NONE
      LOGICAL	RETURN_VAL,belongs
      INTEGER i,origin	
      RETURN_VAL = .FALSE.	 	  
      IF(mpi_task_defined) THEN
      IF(int(myid/mpi_task_per_traject)*mpi_task_per_traject.eq.myid)
     & THEN
      RETURN_VAL = .FALSE.	  
      IF(write_check_file_defined) THEN
      TIME_DOTRAJECT = MPI_Wtime()
      IF((TIME_DOTRAJECT-TIME_BEGIN).gt.60d0*TIME_MIN_CHECK) THEN
      RETURN_VAL = .TRUE.
      INTERUPT_PROPAGATION = RETURN_VAL	  
      ENDIF      
      ENDIF

      ENDIF
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	   
      origin = 0		  
      CALL MPI_Bcast(INTERUPT_PROPAGATION, 1,
     & MPI_LOGICAL, origin, comms(i),ierr_mpi)
      CALL MPI_Bcast(RETURN_VAL, 1,
     & MPI_LOGICAL, origin, comms(i),ierr_mpi)	 
      ENDIF
      ENDDO
      ELSE
      RETURN_VAL = .FALSE.	  
      IF(write_check_file_defined) THEN
      TIME_DOTRAJECT = MPI_Wtime()
      IF((TIME_DOTRAJECT-TIME_BEGIN).gt.60d0*TIME_MIN_CHECK) THEN
      RETURN_VAL = .TRUE.
      INTERUPT_PROPAGATION = RETURN_VAL	  
      
      ENDIF      
      ENDIF
      ENDIF
	  
	  
      END SUBROUTINE TIME_CHECKER_ODEINT	  
	  
	  
      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      INTEGER iest,nv,IMAX,NMAX
      REAL*8 xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50000)
      INTEGER j,k1
      REAL*8 delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END	  
	  
      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      INTEGER nstep,nvar,NMAX
      REAL*8 htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=50000)
      INTEGER i,n
      REAL*8 h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END

      !This subroutine checks whether a variable overflows a 32-bit integer. 
	  !I added it for troubleshooting.
	  SUBROUTINE MQCT_ASSIGN_MPICNT(mpicnt, true_count, tag)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: mpicnt					! 32 Bit Output value
      INTEGER*8, INTENT(IN) :: true_count				! 64 Bit Output value
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tag		! Identifier for logging
      INTEGER*8 :: max_mpi_cnt							
      max_mpi_cnt = INT(HUGE(mpicnt), 8)				! The upper bound 2^31
      IF (true_count .LT. 0_8) THEN						! Check negative count
        IF (PRESENT(tag)) THEN
          WRITE(*,'(A,A)') 'MPICNT negative ', TRIM(tag)
        ELSE
          WRITE(*,'(A)') 'MQCT_ASSIGN_MPICNT: negative count'
        ENDIF
        STOP 'MQCT_ASSIGN_MPICNT'
      ENDIF
      IF (true_count .GT. max_mpi_cnt) THEN				! Check Integer Overflow 
        IF (PRESENT(tag)) THEN
          WRITE(*,'(A,A,1X,I0)') 'MPICNT too large ',
     &     TRIM(tag), true_count
        ELSE
          WRITE(*,'(A,I0)') 'MPICNT too large ',
     &     true_count
        ENDIF
        STOP 'MQCT_ASSIGN_MPICNT'
      ENDIF
      mpicnt = INT(true_count)							! Convert to int-32
      END SUBROUTINE MQCT_ASSIGN_MPICNT

      !This subroutine broadcasts int-32 variables
	  SUBROUTINE MQCT_BCAST_I8(scalar, root, comm, ierr)
      USE MPI
      IMPLICIT NONE
      INTEGER*8, INTENT(INOUT) :: scalar
      INTEGER, INTENT(IN) :: root, comm
      INTEGER, INTENT(OUT) :: ierr
      INTEGER :: parts(2)
      INTEGER*8 :: mask32
      mask32 = 4294967295_8								! (2^32-1) Fortran counts from 1 to 32
      parts(1) = INT(IAND(scalar, mask32))				! Holds the first 32-bit of DATA
      parts(2) = INT(ISHFT(scalar, -32))				! Holds the second 32-bit of DATA
      CALL MPI_BCAST(parts, 2, MPI_INTEGER, root, comm, ierr) !Broadcasting
      scalar = IAND(INT(parts(1), 8), mask32)			! Restoring initial value
     &     + ISHFT(IAND(INT(parts(2), 8), mask32), 32)
      END SUBROUTINE MQCT_BCAST_I8

      !This subroutine finds the largest int-32 value
	  SUBROUTINE MQCT_MPI_MAX_CNT8(max_cnt)
      IMPLICIT NONE
      INTEGER :: mpicnt
      INTEGER*8, INTENT(OUT) :: max_cnt
      max_cnt = INT(HUGE(mpicnt), 8)
      END SUBROUTINE MQCT_MPI_MAX_CNT8

      !This subroutine is needed to broadcast large REAL*8 arrays (root to all ranks)
	  SUBROUTINE MQCT_MPI_BCAST_R8(buf, n_total, root, comm, ierr, tag)
      USE MPI
      IMPLICIT NONE
      
	  REAL*8 buf(*)										! Array that contains REAL*8 values
	  INTEGER*8, INTENT(IN) :: n_total					! Total number elements in array
      INTEGER, INTENT(IN) :: root, comm					! Broadcasting root rank and MPI communicator
      INTEGER, INTENT(OUT) :: ierr						! MPI error
      
	  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tag
      
	  INTEGER*8 offset, n_left, n_chunk, max_cnt        ! Loop counters and chunk sizing
      INTEGER mpicnt, pos                               ! Chunk size in 32-bit and array offset index
      
	  IF (n_total .LE. 0_8) RETURN						! If the total number is zero return
      
	  CALL MQCT_MPI_MAX_CNT8(max_cnt)
      offset = 0_8
      
	  DO WHILE (offset .LT. n_total)					! Cycle through an array
        n_left = n_total - offset						
        n_chunk = n_left
        IF (n_chunk .GT. max_cnt) n_chunk = max_cnt		! Limit the chunk size to the maximum allowed for int-32
        pos = INT(offset + 1_8)							! Index where to start broadcast in buf
        CALL MQCT_ASSIGN_MPICNT(mpicnt, n_chunk, tag)
        CALL MPI_BCAST(buf(pos), mpicnt, MPI_REAL8, root, comm, ierr)
        offset = offset + n_chunk						! Jump to the next chunk, skipping the offset in between
      ENDDO
      
	  END SUBROUTINE MQCT_MPI_BCAST_R8
	  
	  !This subroutine is needed to SEND large REAL*8 arrays (one sender to one receiver)
	  !Similar to MQCT_MPI_BCAST_R8
      SUBROUTINE MQCT_MPI_SEND_R8(buf, n_total, dest, stag, comm,
     & ierr, tag)
      USE MPI
      IMPLICIT NONE
      REAL*8 buf(*)
      INTEGER*8, INTENT(IN) :: n_total
      INTEGER, INTENT(IN) :: dest, stag, comm
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tag
      INTEGER*8 offset, n_left, n_chunk, max_cnt
      INTEGER mpicnt, pos
      IF (n_total .LE. 0_8) RETURN
      CALL MQCT_MPI_MAX_CNT8(max_cnt)
      offset = 0_8
      DO WHILE (offset .LT. n_total)
        n_left = n_total - offset
        n_chunk = n_left
        IF (n_chunk .GT. max_cnt) n_chunk = max_cnt
        pos = INT(offset + 1_8)
        CALL MQCT_ASSIGN_MPICNT(mpicnt, n_chunk, tag)
        CALL MPI_SEND(buf(pos), mpicnt, MPI_REAL8, dest, stag,
     &   comm, ierr)
        offset = offset + n_chunk
      ENDDO
      END SUBROUTINE MQCT_MPI_SEND_R8

      !This subroutine is needed to RECEIVE large REAL*8 arrays (Receive-partner)
	  !Similar to MQCT_MPI_BCAST_R8
	  SUBROUTINE MQCT_MPI_RECV_R8(buf, n_total, src, stag, comm,
     & status, ierr, tag)
      USE MPI
      IMPLICIT NONE
      REAL*8 buf(*)
      INTEGER*8, INTENT(IN) :: n_total
      INTEGER, INTENT(IN) :: src, stag, comm
      INTEGER status(*)
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tag
      INTEGER*8 offset, n_left, n_chunk, max_cnt
      INTEGER mpicnt, pos
      IF (n_total .LE. 0_8) RETURN
      CALL MQCT_MPI_MAX_CNT8(max_cnt)
      offset = 0_8
      DO WHILE (offset .LT. n_total)
        n_left = n_total - offset
        n_chunk = n_left
        IF (n_chunk .GT. max_cnt) n_chunk = max_cnt
        pos = INT(offset + 1_8)
        CALL MQCT_ASSIGN_MPICNT(mpicnt, n_chunk, tag)
        CALL MPI_RECV(buf(pos), mpicnt, MPI_REAL8, src, stag,
     &   comm, status, ierr)
        offset = offset + n_chunk
      ENDDO
      END SUBROUTINE MQCT_MPI_RECV_R8

      !This subroutine is needed to broadcast large integer arrays (root to all ranks)
	  !Similar to MQCT_MPI_BCAST_R8
	  SUBROUTINE MQCT_MPI_BCAST_I4(buf, n_total, root, comm, ierr, tag)
      USE MPI
      IMPLICIT NONE
      INTEGER buf(*)
      INTEGER*8, INTENT(IN) :: n_total
      INTEGER, INTENT(IN) :: root, comm
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tag
      INTEGER*8 offset, n_left, n_chunk, max_cnt
      INTEGER mpicnt, pos
      IF (n_total .LE. 0_8) RETURN
      CALL MQCT_MPI_MAX_CNT8(max_cnt)
      offset = 0_8
      DO WHILE (offset .LT. n_total)
        n_left = n_total - offset
        n_chunk = n_left
        IF (n_chunk .GT. max_cnt) n_chunk = max_cnt
        pos = INT(offset + 1_8)
        CALL MQCT_ASSIGN_MPICNT(mpicnt, n_chunk, tag)
        CALL MPI_BCAST(buf(pos), mpicnt, MPI_INTEGER, root, comm,
     &   ierr)
        offset = offset + n_chunk
      ENDDO
      END SUBROUTINE MQCT_MPI_BCAST_I4

      !This subroutine is needed to broadcast large INTEGER*8 arrays (root to all ranks)
	  SUBROUTINE MQCT_MPI_BCAST_I8(buf, n_total, root, comm, ierr, tag)
      USE MPI
      IMPLICIT NONE
      INTEGER*8 buf(*)
      INTEGER*8, INTENT(IN) :: n_total
      INTEGER, INTENT(IN) :: root, comm
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tag
      INTEGER*8 offset, n_left, n_chunk, max_cnt
      INTEGER mpicnt, pos
      IF (n_total .LE. 0_8) RETURN
      CALL MQCT_MPI_MAX_CNT8(max_cnt)
      offset = 0_8
      DO WHILE (offset .LT. n_total)
        n_left = n_total - offset
        n_chunk = n_left
        IF (n_chunk .GT. max_cnt) n_chunk = max_cnt
        pos = INT(offset + 1_8)
        CALL MQCT_ASSIGN_MPICNT(mpicnt, n_chunk, tag)
        CALL MPI_BCAST(buf(pos), mpicnt, MPI_INTEGER8, root, comm,
     &   ierr)
        offset = offset + n_chunk
      ENDDO
      END SUBROUTINE MQCT_MPI_BCAST_I8

      ! Gathers one INTEGER*8 value from each rank to recvbuf on root
	  SUBROUTINE MQCT_MPI_GATHER_I8(sendbuf, recvbuf, root, comm, ierr)
      USE MPI
      IMPLICIT NONE
      INTEGER*8, INTENT(IN) :: sendbuf
      INTEGER*8 recvbuf(*)
      INTEGER, INTENT(IN) :: root, comm
      INTEGER, INTENT(OUT) :: ierr
      CALL MPI_GATHER(sendbuf, 1, MPI_INTEGER8, recvbuf, 1,
     & MPI_INTEGER8, root, comm, ierr)
      END SUBROUTINE MQCT_MPI_GATHER_I8

      !This subroutine is needed to broadcast large logical arrays (root to all ranks)
	  !Similar to MQCT_MPI_BCAST_R8
	  SUBROUTINE MQCT_MPI_BCAST_LOG(buf, n_total, root, comm, ierr, tag)
      USE MPI
      IMPLICIT NONE
      LOGICAL buf(*)
      INTEGER*8, INTENT(IN) :: n_total
      INTEGER, INTENT(IN) :: root, comm
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tag
      INTEGER*8 offset, n_left, n_chunk, max_cnt
      INTEGER mpicnt, pos
      IF (n_total .LE. 0_8) RETURN
      CALL MQCT_MPI_MAX_CNT8(max_cnt)
      offset = 0_8
      DO WHILE (offset .LT. n_total)
        n_left = n_total - offset
        n_chunk = n_left
        IF (n_chunk .GT. max_cnt) n_chunk = max_cnt
        pos = INT(offset + 1_8)
        CALL MQCT_ASSIGN_MPICNT(mpicnt, n_chunk, tag)
        CALL MPI_BCAST(buf(pos), mpicnt, MPI_LOGICAL, root, comm,
     &   ierr)
        offset = offset + n_chunk
      ENDDO
      END SUBROUTINE MQCT_MPI_BCAST_LOG

!Divides range into balanced indices across MPI ranks  
      SUBROUTINE MQCT_MPI_CHUNK(rank_id, nproc, n, lo, hi,
     & if_n_lt_np_empty)
      IMPLICIT NONE
      
	  ! Input parameters: rank ID, number of processors, and total size
	  INTEGER, INTENT(IN) :: rank_id, nproc
      INTEGER*8, INTENT(IN) :: n
      
	  ! Output parameters: start (lo) and end (hi) indices
	  INTEGER*8, INTENT(OUT) :: lo, hi
      
	  ! Flag to control when total size is less than the number of processors
	  LOGICAL, INTENT(IN) :: if_n_lt_np_empty
      
	  ! Local variables for chunk size, remainder, and int64 conversions
	  INTEGER*8 chunk, rem, rid
      INTEGER np8
      
	  ! Convert inputs to 64-bit integers to prevent overflow
	  np8 = INT(nproc, 8)
      rid = INT(rank_id, 8)
      
	  ! Return an empty range if number of processors or total size is invalid
	  IF (nproc .LT. 1 .OR. n .LT. 1_8) THEN
        lo = 1_8
        hi = 0_8
        RETURN
      ENDIF
      
	  ! Return an empty range size < ranks
	  IF (if_n_lt_np_empty .AND. n .LT. np8) THEN
        lo = 1_8
        hi = 0_8
        RETURN
      ENDIF
      
	  ! Calculate chunk size per rank and remainder
	  chunk = n / np8
      rem = n - chunk*np8
      
	  ! Indices computation assuming even distribution
	  lo = rid*chunk + 1_8
      hi = (rid + 1_8)*chunk
	  
	  
      IF (rid .LT. rem) THEN
        lo = lo + rid					! Shift start ind by rank id
        hi = hi + rid + 1_8				! Shift end ind by rank id plus one
      ELSE
        lo = lo + rem					! Shift start ind by remainder
        hi = hi + rem					! Shift end ind by remainder
      ENDIF
      
	  ! Check the upper bound does not exceed total items n
	  IF (hi .GT. n) hi = n
	  
	  ! Normalize empty ranges ( to avoid negative element count further in the loop)
	  ! Normalizion high = low - 1 makes sure that loop will be zero and code will not crush
	  ! To understand study an example hi=0 and lo=10 for unnormalized and normalizeed
      IF (lo .GT. hi) hi = lo - 1_8
      END SUBROUTINE MQCT_MPI_CHUNK

      ! For given global index k, which two states (st1, st2) does it refer to
	  SUBROUTINE GLOBAL_K_TO_PAIR(k, st1, st2)
      
	  USE VARIABLES
      IMPLICIT NONE
      
	  ! Input parameters: global index
	  INTEGER*8, INTENT(IN) :: k
      
	  ! Output parameters: state pair indices
	  INTEGER, INTENT(OUT) :: st1, st2
      
	  ! Local Offset and Local Column Count 
	  INTEGER*8 loc, nloc						 				
      
	  ! Check if the index matrix is distributed across processes
	  IF (ind_mat_local_defined) THEN							
        
		! Convert global index k to local index loc
		loc = k - ind_mat_k0 + 1_8								
        nloc = SIZE(ind_mat, 2, 8)
        
		! Verify local index falls within valid range
		IF (loc .LT. 1_8 .OR. loc .GT. nloc) THEN				
          WRITE(*,'(A,I0,A,I0,A,I0)') 'GLOBAL_K_TO_PAIR: k=', k,
     &     ' outside local [', ind_mat_k0, ', ...]'
          STOP 'GLOBAL_K_TO_PAIR'
        ENDIF
        
		! Read state pair from the local chunk
		st1 = ind_mat(1, loc)
        st2 = ind_mat(2, loc)
      
	  ELSE
        
		! Read state pair directly
		st1 = ind_mat(1, k)
        st2 = ind_mat(2, k)
      ENDIF
      
	  END SUBROUTINE GLOBAL_K_TO_PAIR

      ! If we are given a state pair (st1, st2) then the this subroutine returns the global coupling index k
	  ! We use it at propagation step 
	  SUBROUTINE PAIR_TO_GLOBAL_K(st1, st2, k)
      USE VARIABLES
      IMPLICIT NONE
      
	  ! Input parameters: state pair indices
	  INTEGER, INTENT(IN) :: st1, st2
      
	  ! Output parameters: global index
	  INTEGER*8, INTENT(OUT) :: k
      
	  ! Loop counter and allowance flag
	  INTEGER st2w
      LOGICAL allowed
	  
	  ! Check for cumulative array allocation
      IF (.NOT. allocated(cum_global_k)) THEN						
        WRITE(*,'(A)') 'PAIR_TO_GLOBAL_K: cum_global_k not allocated'
        STOP 'PAIR_TO_GLOBAL_K'
      ENDIF
	  
	  ! We return the k value before st1
      k = cum_global_k(st1)										 
      
	  ! Iterate through pairs (st1, st2) up to st1
	  DO st2w = 1, st1
        CALL PAIR_IS_ALLOWED(st1, st2w, allowed)				! Subroutine which checks if the pair is allowed or not
        IF (.NOT. allowed) CYCLE								! If not allowed we skip
        k = k + 1_8												! k = k + 1
        IF (st2w .EQ. st2) RETURN								! Return the last value of k
      ENDDO
      
	  ! Write an error and abort
	  WRITE(*,'(A,2I8)') 'PAIR_TO_GLOBAL_K: pair not found ', st1, st2
      STOP 'PAIR_TO_GLOBAL_K'
      END SUBROUTINE PAIR_TO_GLOBAL_K

      ! There are situations when code crashed so I added it to abort
	  ! I will fix these errors step-by-step
	  SUBROUTINE GUARD_DISTRIB_IND_MAT(feature)
      USE VARIABLES
      USE MPI_DATA
      USE MPI
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: feature
      IF (.NOT. dist_ind_mat_write) RETURN
      IF (MYID.eq.0) THEN
        WRITE(*,'(A,A,A)') 'ERROR: ', TRIM(feature),
     &   ' not supported with distributed ind_mat (dist_ind_mat_write)'
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)
      CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr_mpi)
      END SUBROUTINE GUARD_DISTRIB_IND_MAT

      ! Check if a transition allowed 
	  SUBROUTINE PAIR_IS_ALLOWED(st1, st2, allowed)
      USE VARIABLES
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: st1, st2
      LOGICAL, INTENT(OUT) :: allowed
      allowed = .FALSE.
      IF (st1 .LT. 1 .OR. st2 .LT. 1) RETURN					! st1 < 1 or st2 < 1 ERROR
      IF (identical_particles_defined) THEN						! If identical
        IF (parity_state(st1) .NE. parity_state(st2)) RETURN	! If they are not equal then we retrun FALSE
      ENDIF
      IF (m12(st1) .EQ. m12(st2)) allowed = .TRUE.
      END SUBROUTINE PAIR_IS_ALLOWED

      ! Total size calculation in parallel
	  SUBROUTINE INIT_IND_MAT_LAYOUT(nstates)
      USE VARIABLES
      USE MPI_DATA
      USE MPI
      IMPLICIT NONE
      
	  ! Input parameters: Number of states
	  INTEGER, INTENT(IN) :: nstates						
      
	  ! Local loop counters, bounds, and process variables
	  INTEGER st1_lo, st1_hi, nloc, p, st1, st2, nrow
      LOGICAL allowed
      INTEGER*8 st1_lo8, st1_hi8
      
	  ! Seg(:) stores chunk of the total size
	  INTEGER, ALLOCATABLE :: seg(:)
	  
      ! Recvcounts states the size of seg(:) on each processor for MPI_ALLGATHERV; displs where to allocate seg(:) chunks for MPI_ALLGATHERV 
	  INTEGER, ALLOCATABLE :: recvcounts(:), displs(:)		
      
	  ! Reallocate row count and cumulative
	  IF (allocated(nrow_for_st1)) DEALLOCATE(nrow_for_st1)
      IF (allocated(cum_global_k)) DEALLOCATE(cum_global_k)
      ALLOCATE(nrow_for_st1(nstates))
      nrow_for_st1 = 0										! Fill nrow_for_st1 with 0
      
	  ! Partition total states (st1) across MPI ranks
	  CALL MQCT_MPI_CHUNK(MYID, nproc, INT(nstates, 8), st1_lo8,
     & st1_hi8, .FALSE.)
      st1_lo = INT(st1_lo8)
      st1_hi = INT(st1_hi8)
      nloc = st1_hi - st1_lo + 1							! Local chunk size 
      
	  ! Allocate local buffer and count allowed pair
	  IF (nloc .LT. 1) THEN									! If the chunk size is less then 0 -> error
        ALLOCATE(seg(1))
        seg(1) = 0
        nloc = 1
      ELSE
        ALLOCATE(seg(nloc))									! Local seg(:) allocated 
        DO st1 = st1_lo, st1_hi
          nrow = 0
          DO st2 = 1, st1
            CALL PAIR_IS_ALLOWED(st1, st2, allowed)
            IF (allowed) nrow = nrow + 1
          ENDDO
          seg(st1-st1_lo+1) = nrow
        ENDDO
      ENDIF
      
	  ! Gather local row counts from all ranks into global nrow_for_st1 array
	  IF (nproc .GT. 1) THEN
        ALLOCATE(recvcounts(nproc), displs(nproc))
        
		! Recompute chunk bounds and receive counts for each rank
		DO p = 0, nproc - 1									! It probably could have been more efficient but I just recalculate to draw boundaries for MPI_ALLGATHERV
          CALL MQCT_MPI_CHUNK(p, nproc, INT(nstates, 8), st1_lo8,
     &     st1_hi8, .FALSE.)
          st1_lo = INT(st1_lo8)
          st1_hi = INT(st1_hi8)
          recvcounts(p+1) = st1_hi - st1_lo + 1
          IF (recvcounts(p+1) .LT. 1) recvcounts(p+1) = 1
        ENDDO
        
		! Buffer displacements for MPI_ALLGATHERV
		displs(1) = 0
        DO p = 2, nproc
          displs(p) = displs(p-1) + recvcounts(p-1)
        ENDDO
        
		! Gather row counts into nrow_for_st1 across all ranks
		CALL MPI_ALLGATHERV(seg, nloc, MPI_INTEGER, nrow_for_st1,
     &   recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr_mpi)
        DEALLOCATE(recvcounts, displs)
      ELSE
        IF (nloc .GE. nstates) THEN							! Additional safeguard
          nrow_for_st1 = seg(1:nstates)
        ENDIF
      ENDIF
      DEALLOCATE(seg)										! All the seg(:) are deallocated
      ALLOCATE(cum_global_k(nstates+1))						! cum_global_k is cumulative of values inside nrow_for_st1
      cum_global_k(1) = 0_8									
      DO st1 = 1, nstates									! Loop to find cumulative
        cum_global_k(st1+1) = cum_global_k(st1)
     &   + INT(nrow_for_st1(st1), 8)
      ENDDO
      total_size = cum_global_k(nstates+1)
      END SUBROUTINE INIT_IND_MAT_LAYOUT

      !This subroutine makes that each rank only iterates the st1 values needed for its (k_st;k_fn) slice when building local ind_mat
	  SUBROUTINE FILL_LOCAL_IND_MAT(k_st, k_fn, nfilled)
      USE VARIABLES
      IMPLICIT NONE
      INTEGER*8, INTENT(IN) :: k_st, k_fn				! Starting and ending index
      INTEGER*8, INTENT(OUT) :: nfilled					! Number of elements filled
      INTEGER*8 nloc, k_run, loc, k_q					
      INTEGER st1, st2, st1_lo, st1_hi, lo, hi, mid, nst
      LOGICAL allowed, found
      nloc = k_fn - k_st + 1_8							! Range length for allocation
      IF (nloc .LE. 0_8) THEN							! Error if nloc < 0 
        nfilled = 0
        IF (allocated(ind_mat)) DEALLOCATE(ind_mat)
        ind_mat_local_defined = .FALSE.
        RETURN
      ENDIF
      IF (allocated(ind_mat)) DEALLOCATE(ind_mat)
      ALLOCATE(ind_mat(2, nloc))						! Allocation for the local index matrix for current range size
      ind_mat = 0
      nst = SIZE(nrow_for_st1)
      
	  ! Binary search cum_global_k to locate starting state st1_lo containing k_st
	  k_q = k_st
      lo = 1
      hi = nst
      found = .FALSE.
      DO WHILE (lo .LE. hi)
        mid = (lo + hi)/2
        IF (cum_global_k(mid+1) .LT. k_q) THEN
          lo = mid + 1
        ELSE IF (cum_global_k(mid) .GE. k_q) THEN
          hi = mid - 1
        ELSE
          st1_lo = mid
          found = .TRUE.
          EXIT
        ENDIF
      ENDDO
      IF (.NOT. found) st1_lo = lo
      
	  ! Binary search cum_global_k to locate ending state st1_hi containing k_fn
	  k_q = k_fn
      lo = 1
      hi = nst
      found = .FALSE.
      DO WHILE (lo .LE. hi)
        mid = (lo + hi)/2
        IF (cum_global_k(mid+1) .LT. k_q) THEN
          lo = mid + 1
        ELSE IF (cum_global_k(mid) .GE. k_q) THEN
          hi = mid - 1
        ELSE
          st1_hi = mid
          found = .TRUE.
          EXIT
        ENDIF
      ENDDO
      IF (.NOT. found) st1_hi = lo
      k_run = 0_8
      
	  ! Fill ind_mat for states in the global range
	  DO st1 = 1, SIZE(nrow_for_st1)					! Iteration over state rows containing target range
        IF (st1 .LT. st1_lo) THEN						! Fastforward past prior states 
          k_run = cum_global_k(st1+1)
          CYCLE
        ENDIF
        IF (st1 .GT. st1_hi) EXIT						! Additional safeguard
        k_run = cum_global_k(st1)						! Initial k for the following chunk 
        DO st2 = 1, st1
          CALL PAIR_IS_ALLOWED(st1, st2, allowed)		! Check selection rules
          IF (.NOT. allowed) CYCLE
          k_run = k_run + 1_8
          IF (k_run .LT. k_st .OR. k_run .GT. k_fn) CYCLE 	! Skip elements outside requested range
          loc = k_run - k_st + 1_8						! Allocation
          ind_mat(1, loc) = st1
          ind_mat(2, loc) = st2
        ENDDO
      ENDDO
      nfilled = nloc
      ind_mat_k0 = k_st
      ind_mat_local_defined = .TRUE.					! We state that ind_mat array was calculated in parallel
      END SUBROUTINE FILL_LOCAL_IND_MAT
