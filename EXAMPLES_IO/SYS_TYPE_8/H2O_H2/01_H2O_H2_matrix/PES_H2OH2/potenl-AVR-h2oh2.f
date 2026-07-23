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

      if (ICNTRL.eq.0) goto  1000                                                                                                                               

c-------------------------------------------------------------------
c Initialization (ICNTRL=-1)
c-------------------------------------------------------------------
      if (ICNTRL.ne.-1) stop 'invalid ICNTRL, expected -1'
      dir = './H2Ov000_H2v0/'
      OPEN(UNIT=6,FILE="PESH2OH2_VALIRON_OUT.out")
      corr_9d_flag = .true.
      r12_flag = .true.
c      read(*,POTL)
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
      write(*,*) 'p(1)'
      NDIM = NTT

c MLD end change

      write(6,*)
      write(6,*) 'Initialization done.'
      write(6,*)
      return      


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


