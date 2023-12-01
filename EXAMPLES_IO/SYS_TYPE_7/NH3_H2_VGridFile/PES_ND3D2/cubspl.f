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
 
 
