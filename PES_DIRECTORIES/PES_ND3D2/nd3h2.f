c  Calculates ND3-H2 potential by shifting the NH3 origin in the NH3-H2 potential
c
      subroutine nd3h2pot(r1, thx1, phx, thpx, phpx, y)
      implicit none
      real*8 mN, mH, mD, zeta, zeta1, delta1, rNH
      real*8 x, z, delta, rho, pi
      real*8 r, thx, phx, thpx, phpx, y
      real*8 r1, thx1
      parameter (mN = 14.0030740052d0, mH = 1.0078250321d0)
      parameter (mD = 2.014101778d0)
      pi = dacos(-1d0)
c     rNH = 1.9098839d0
c     rho = 112.15*pi/180d0
c zie berekende NH3-H2 potentiaal in MNRAS 399, 425 (2009)
c and Mol. Phys. 102, 2297 (2004)
      rNH = 1.9512d0        !
      rho = 111.49*pi/180d0
      zeta = mN/(3*mH+mN)
c distance between N atom and center of mass of NH3
      delta = rNH*dcos(rho)*(1d0-zeta)
      zeta1 = mN/(3*mD+mN)
c distance between N atom and center of mass of ND3
      delta1 = rNH*dcos(rho)*(1d0-zeta1)  ! more negative than delta
      x = r1*dsin(thx1*pi/180d0)
      z = r1*dcos(thx1*pi/180d0) + delta1 - delta
      r = dsqrt(x*x + z*z)
      thx = dacos(z/r)*180d0/pi
      call nh3h2pot(r, thx, phx, thpx, phpx, y)
      return
      end
c     This program is a slightly modified version of fit5d.f
c     from Valiron et al for H2O-H2
c     It allows to obtain an angular description of the NH3-H2 potential
c     at a distance supplied by the user.
c     The current default is phi=phi'=0
c     The output is in fort.8
c
c     A. Faure Jan 2009
c Input: R in bohr, angles in degrees
c
      subroutine nh3h2pot(r, thx, phx, thpx, phpx, y)
      implicit none
      integer ndim, i, j
      parameter(ndim=190)
      integer icntrl, nlam, lam(4,ndim), p1, q1, p2, p
      character dir*80
      real*8 R, COEF(ndim), rdummy,
     $     thx, phx, thpx, phpx, tnormed, y, pi, factor
      logical first
      data first/.true./
      save
      Pi=acos(-1.0d0)

c-------------------------------- INITIALIZATION
      if (first) then
        icntrl=-1
        dir='/mmfs1/home/4168mandalb/MQCT/v20_21/PES_ND3D2/FIT120/'
        rdummy=0d0
        do i=1,ndim
           coef(i)=0d0
        enddo
        call nh3h2_5d(icntrl, dir, ndim, nlam, lam, rdummy, coef)
        write(6,*) 'NLAM', nlam
        first = .false.
        return
      endif

c-------------------------------- PES EVALUATION

      icntrl=0
      call nh3h2_5d(icntrl, dir, ndim, nlam, lam, r, coef)
      y=0
      do i=1,nlam
         p1=lam(1,i)
         q1=lam(2,i)
         p2=lam(3,i)
         p =lam(4,i)

         y=y+coef(i)*tnormed(p1,q1,p2,p,thx,phx,thpx,phpx)

      enddo

      end

      subroutine NH3H2_5d(ICNTRL, DIR, NDIM, NLAM, LAM, RD, P)
c
c     -------------------------------------------------------------------
c     Radial fit of CCSD(T) PES for 5-D NH3-H2
c     Combines expansion coefficients for two PES.
c       d: aug-cc-pVDZ+bf PES, 3000 geometries
c       t: aug-cc-pVTZ+bf PES, 1000 geometries
c     Expansion coefficients are read for
c       1) d expansion
c       2) dt-d expansion, where dt is a CBS extrapolation od d and t PES
c     Original source by PV
c     Adapted for present NH3-H2 PES by PV, march & april 2006.
c     -------------------------------------------------------------------
c
c  ARGUMENT LIST
c
c  ICNTRL      (input)   -1  mandatory initialization
c                         0  to evaluate P(R)
c  DIR         (input)   slash/ terminated directory containing data   (-1)
c  NDIM        (input)   dimension of LAM(4,*) and P(*)
c  NLAM        (output)  number of terms of the expansion
c  LAM(4,NDIM) (output)  (4,NDIM) p1, q1, p2, p expansion terms        (-1)
c  RD          (input)   intermolecular distance (a.u.)                (0)
c  P(NDIM)     (output)  coefficients for angular expansion in cm-1    (0)
c              properly interpolated or extrapolated at distance R
c              assuming normalized expansion terms. Normalization
c              should follow MOLSCAT's conventions, for the
c              corresponding case.
c
c     These coefficients come in the ordering indicated in LAM.  This
c     ordering reflects the list of terms given in the input file
c     mesh.dat and can thus be easily changed.  Terms can easily be
c     reordered, removed or added in mesh.dat provided they form a
c     subset of the original expansion in the fit*.dat files for
c     expansion 1.
c
c     The routine performs numerous consistency checks.
c
c  RELATED FILES
c
c     mesh.dat, fit*.dat, in directory DIR.
c
      implicit none
      integer icntrl, ndim, nlam, lam(4,ndim), np1
      character*(*) dir
      real*8 rd, p(ndim)
      integer ndx, ntx, nttx
      parameter (ndx=200, ntx=150, nttx=150)
      integer ndist, id,
     $     Ngeom1, Nfunc1, Ngeom1x, Nfunc1x,
     $     Ngeom2, Nfunc2, Ngeom2x, Nfunc2x,
     $     i, it, nt1, nt2, k,
     $     l1x(4,ntx), l1(4,ntx), l2x(4,ntx), l2(4,ntx), nex
      integer ind, ind2
      real*8 fdist, Dist(ndx),
     $     Sinv1(ndx), RMS1(ndx), EEST1(ndx),
     $     Sinv2(ndx), RMS2(ndx), EEST2(ndx)
      real*8 coef1(ndx,ntx), coef2(ndx,ntx),
     $     cc(ndx,nttx), c1(ndx,nttx), c2(ndx,nttx), c3(ndx,nttx)
      integer ntt, itt, in, ll(4,nttx), indll1(nttx), indll2(nttx)
      logical chk(ntx), error
      integer n0000, kleft, kright
      real*8 a(nttx), alpha(nttx), dleft, dright, hinvl, hinvr
      real*8 c(nttx), b(nttx), v0000, ccinv0000, fact
      character cfit1*20, cfit2*20, comment*132, name*132
      integer le
      real*8 h, r, rinv, aniso, tres, thres
      real*8 f, u, fu, Pi, Pis2
c     funct(u)=(1-cos(u*Pi))/2 ; f(u)=funct(funct(u))
      f(u)=(1-cos((1-cos(u*Pi))*Pis2))*0.5d0
      save

      if (ICNTRL.eq.0) goto 1000

c-------------------------------------------------------------------
c Initialization (ICNTRL=-1)
c-------------------------------------------------------------------
      if (ICNTRL.ne.-1) then
         write(6,*) 'invalid ICNTRL, expected -1'
         call die
      endif
      Pi=acos(-1.0d0)
      Pis2=Pi/2

c Set directory containing the fitted and mesh data (slash terminated)

      le=len_trim(dir)

c Open mesh.dat

c      name=dir(:le)//'mesh-para.dat'
      name=dir(:le)//'mesh.dat'
      open(unit=1, file=name, status='old')

      read(1,*) comment

c Read expansion data in fit*.dat files

      read(1,*) ndist
      if (ndist.gt.ndx) then
         write(6,*) 'increase ndx'
         call die
      endif

      do id=1,ndist

         read(1,*) Dist(id), cfit1, cfit2
         if (id.gt.1) then
            if (Dist(id).le.Dist(id-1)) then
               write(6,*) 'non monotonic Dist'
               call die
            endif
         endif

c-------------------------------------------------------------------
c  Expansion 1
c-------------------------------------------------------------------

         name=dir(:le)//cfit1
         open(unit=2, file=name, status='old')

c Skip two comment lines
         read(2,*)
         read(2,*)

c Read data at dist(id)
         read(2,*) fdist, Ngeom1, Nfunc1,
     $        Sinv1(id), RMS1(id), EEST1(id)
         if (dist(id).ne.fdist) then
            write(6,*) 'Inconsistent distance in fit*.dat'
            write(6,*) fdist, dist(id)
            call die
         endif
         if (id.eq.1) then
            Ngeom1x=Ngeom1
            Nfunc1x=Nfunc1
            nt1=Nfunc1
            if (Nfunc1.gt.ntx) then
               write(6,*) 'increase ntx'
               call die
            endif
         else
            if (Ngeom1x.ne.Ngeom1) then
               write(6,*) 'inconsistent Ngeom1'
               call die
            endif
            if (Nfunc1x.ne.Nfunc1) then
               write(6,*) 'inconsistent Nfunc1'
               call die
            endif
         endif

         do it=1,nt1
            read(2,*) (l1(i,it),i=1,4), coef1(id,it)
         enddo

         if (id.eq.1) then
            do it=1,nt1
               do i=1,4
                  l1x(i,it)=l1(i,it)
               enddo
            enddo
         else
            k=0
            do it=1,nt1
               do i=1,4
                  k=k+abs(l1x(i,it)-l1(i,it))
               enddo
            enddo
            if (k.ne.0) then
               write(6,*) 'inconsistent expansion terms'
               call die
            endif
         endif

         close(unit=2)

c-------------------------------------------------------------------
c  Expansion 2
c-------------------------------------------------------------------

         name=dir(:le)//cfit2
         open(unit=2, file=name, status='old')

c Skip two comment lines
         read(2,*)
         read(2,*)

c Read data at dist(id)
         read(2,*) fdist, Ngeom2, Nfunc2,
     $        Sinv2(id), RMS2(id), EEST2(id)
         if (dist(id).ne.fdist) then
            write(6,*) 'Inconsistent distance in ref.DDDD.dat'
            write(6,*) fdist, dist(id)
            call die
         endif
         if (id.eq.1) then
            Ngeom2x=Ngeom2
            Nfunc2x=Nfunc2
            nt2=Nfunc2
            if (Nfunc2.gt.ntx) then
               write(6,*) 'increase ntx'
               call die
            endif
         else
            if (Ngeom2x.ne.Ngeom2) then
               write(6,*) 'inconsistent Ngeom2'
               call die
            endif
            if (Nfunc2x.ne.Nfunc2) then
               write(6,*) 'inconsistent Nfunc2'
               call die
            endif
         endif

         do it=1,nt2
            read(2,*) (l2(i,it),i=1,4), coef2(id,it)
         enddo

         if (id.eq.1) then
            do it=1,nt2
               do i=1,4
                  l2x(i,it)=l2(i,it)
               enddo
            enddo
         else
            k=0
            do it=1,nt2
               do i=1,4
                  k=k+abs(l2x(i,it)-l2(i,it))
               enddo
            enddo
            if (k.ne.0) then
               write(6,*) 'inconsistent expansion terms'
               call die
            endif
         endif

         close(unit=2)

         if (id.eq.1)
     $        write(6,'(a7,2(a6,a4,3a8))') 'R  ',
     $        'Geo1', 'Nf', 'Sinv', 'rms', 'eest',
     $        'Geo2', 'Nf', 'Sinv', 'rms', 'eest'

         write(6,'(f7.3,2(i6,i4,3f8.3))') fdist,
     $        Ngeom1, Nfunc1, Sinv1(id), RMS1(id), EEST1(id),
     $        Ngeom2, Nfunc2, Sinv2(id), RMS2(id), EEST2(id)

      enddo

c____________________________________________________________________

c Select final expansion terms and corresponding pointers to original data

      read(1,*) ntt
      if (ntt.gt.nttx) then
         write(6,*) 'increase nttx'
         call die
      endif
      read(1,*) ((ll(i,itt), i=1,4), itt=1,ntt)

      do itt=1,nttx
         indll1(itt)=0
         indll2(itt)=0
      enddo

c Match terms for expansion 1
      do it=1,ntx
         chk(it)=.false.
      enddo
      do itt=1,ntt
         do it=1,nt1
            if (       l1(1,it).eq.ll(1,itt)
     $           .and. l1(2,it).eq.ll(2,itt)
     $           .and. l1(3,it).eq.ll(3,itt)
     $           .and. l1(4,it).eq.ll(4,itt) ) then
               if (chk(it)) then
                  write(6,*) 'l1=', (l1(i,it),i=1,4)
                  write(6,*) 'll=', (ll(i,itt),i=1,4)
                  write(6,*) 'itt=', itt, '     it=', it
                  write(6,*) 'error -- duplicate term in expansion 1'
                  call die
               else
                  indll1(itt)=it
                  chk(it)=.true.
               endif
            endif
         enddo
      enddo

c Match terms for expansion 2
      do it=1,ntx
         chk(it)=.false.
      enddo
      do itt=1,ntt
         do it=1,nt2
            if (       l2(1,it).eq.ll(1,itt)
     $           .and. l2(2,it).eq.ll(2,itt)
     $           .and. l2(3,it).eq.ll(3,itt)
     $           .and. l2(4,it).eq.ll(4,itt) ) then
               if (chk(it)) then
                  write(6,*) 'l2=', (l2(i,it),i=1,4)
                  write(6,*) 'll=', (ll(i,itt),i=1,4)
                  write(6,*) 'itt=', itt, '     it=', it
                  write(6,*) 'error -- duplicate term in expansion 2'
                  call die
               else
                  indll2(itt)=it
                  chk(it)=.true.
               endif
            endif
         enddo
      enddo

c Test all requested expansion terms are present in the data
      error=.false.
      do itt=1,ntt
         if (indll1(itt).eq.0 .and. indll2(itt).eq.0) then
            error=.true.
            write(6,*)
     $           'missing pointer for term ll=', (ll(i,itt),i=1,4)
         endif
      enddo
      if (error) call die

c Accumulate expansions 1 and 2 values into cc
      do itt=1,ntt
         do id=1,ndist
            cc(id,itt)=0
         enddo
         it=indll1(itt)
         if (it.ne.0) then
            do id=1,ndist
               cc(id,itt)=cc(id,itt)+coef1(id,it)
            enddo
         endif
         it=indll2(itt)
         if (it.ne.0) then
            do id=1,ndist
               cc(id,itt)=cc(id,itt)+coef2(id,it)
            enddo
         endif
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
      if (n0000.eq.0) then
         write(6,*) 'n0000 term missing in expansion'
         call die
      endif

c Retrieve and set up long range terms as C/R**beta

c All terms above thres at both dist(ndist-1) and dist(ndist) are
c extrapolated by a power law obtained from these two largest R.

      do itt=1,ntt
         b(itt)=0
         c(itt)=0
      enddo
      thres=1.d-1
      nex=0
      do itt=1, ntt
         if (abs(cc(ndist-1,itt)).ge.thres
     $        .and. abs(cc(ndist,itt)).ge.thres) then
            nex=nex+1
            b(itt)=log(cc(ndist-1,itt)/cc(ndist,itt))
     $           /log(dist(ndist)/dist(ndist-1))
            c(itt)=cc(ndist,itt)*dist(ndist)**b(itt)
            if (nex.eq.1) then
               write(6,*)
               write(6,*) 'Long range extrapolation:'
               write(6,'(a4,a16,a14,a12)')
     $           'itt', 'Vxxxx  ', 'Coeff', '1/R Power'
            endif
            write(6,'(i4,4x,4i3,1pe14.3,0pf12.3)')
     $           itt, (ll(i,itt),i=1,4), c(itt), b(itt)
         endif
      enddo
      write(6,'(i8,a)') nex, ' terms extrapolated at long range'

c Initialize extrapolation at short range

      write(6,*)
      write(6,*) 'Monitor exponential short range extrapolation'
      tres=3000
      write(6,'(1x,a,1pe12.2,a,0pf8.3)') 'Select terms above TRES =',
     $     tres, '  at R=', dist(1)
      write(6,'(a4,a16,a14,a12)')
     $     'itt', 'Vxxxx  ', 'Coeff', 'Exp fact.'
      do itt=1,ntt
         if (abs(cc(2,itt)) .lt. 1d-10) then
            alpha(itt)=-1
         else if (cc(1,itt)/cc(2,itt) .lt. 1d-10) then
            alpha(itt)=-1
         else
            alpha(itt)=log(cc(1,itt)/cc(2,itt))/(dist(2)-dist(1))
         endif
         a(itt)=cc(1,itt)*exp(alpha(itt)*dist(1))
         if (abs(cc(1,itt)).gt.tres) then
            write(6,'(i4,4x,4i3,1pe14.3,0pf12.3)')
     $           itt, (ll(i,itt),i=1,4), a(itt), alpha(itt)
         endif
      enddo
      write(6,*)
      write(6,*)'    Dist      Isotropic   Sigma(Aniso)    Ratio'
      r=2
      do
         aniso=0
         do itt=1,ntt
            if (itt.ne.n0000 .and. abs(cc(1,itt)).gt.tres) then
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
     $     '==> Anisotropic terms may be large at shortest range;'
      write(6,*)
     $     '==> Blind extrapolation of anisotropic terms is DANGEROUS;'
      write(6,*)
     $     '==> Set up a more secure proportional extrapolation scheme.'

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
c r=dist(ndist)+1 is added to cc.
c
      np1=ndist+1
      if (np1.gt.ndx) then
         write(6,*) 'increase ndx'
         call die
      endif
      dist(np1)=dist(ndist)+1
      do itt=1,ntt
         cc(np1,itt)=c(itt)*(1/dist(np1))**b(itt)
      enddo

c Initialize spline coefficients for original data

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

c-------------------------------------------------------------------

c return proper data
      if (ntt.gt.ndim) then
         write(6,*) 'increase NDIM'
         call die
      endif
      nlam=ntt
      do itt=1,ntt
         do i=1,4
            lam(i,itt)=ll(i,itt)
         enddo
      enddo
      close(unit=1)
      write(6,*)
      write(6,*) 'Initialization done.'
      write(6,*)
      return


c-------------------------------------------------------------------
c Interpolation or extrapolation of the PES expansion (ICNTRL=0)
c-------------------------------------------------------------------

 1000 if (ntt.gt.ndim) then
         write(6,*) 'invalid NDIM'
         call die
      endif
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

      end

      subroutine die
      stop
      end
