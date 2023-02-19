c
c      subroutine poten(R,oa,ob,val,mode)
c cccccccccccc
c
c The following header relevant for the ntpot=161 version with
c derivatives. forceab(3) vector is the force on B due to A,
c torque on A due to B is torquea and torque on B due to A is torqueb.
c
      subroutine poten(R,oa,ob,val,forceab,torquea,torqueb,mode)
      implicit real*8 (a-h,o-z)
c
      parameter (maxb=500,maxp=140000)
      parameter (nsitemax=42,ntypemax=25,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=3,maxpar2=13)
      parameter (mxl=100)
c
      dimension oa(3), ob(3), values(20), valder(20)
      dimension torquea(3), torqueb(3), forceab(3)
c Only for the purpose of checking.....
      dimension valuesp(20),valuesm(20)
c
c
c The following common to pass the vector of linear parameters...
c
      common/leastsq/dat(maxp),a(maxp,maxb),a0(maxp),g(maxp),c(maxb),
     1 ata(maxb*maxb),b(maxb),works(2*maxb),iworks(maxb),np,nb
      common/sites/ sitea(3,nsitemax), siteb(3,nsitemax),nsitea,nsiteb
      common/temp_sites/ siteat(3,nsitemax), sitebt(3,nsitemax),
     1                   tranA(3,3), tranB(3,3)
      common/fit/rpt(6,maxp),dat0(maxp),sigma(maxp),iweight,ntpot,
     1 idonl, iopt
      common/type_used/ itypu(ntypemax,ntypemax)
      common/types/ itypea(nsitemax),itypeb(nsitemax)
      common/misc/ R_0,isyst,npowers,numlin
      data bohr2a /0.529177249d0/
c
c If mode = -1 -- read input
c
      if(mode.eq.-1) then
c       call data
       call data1
c
c Stop for a moment...
c
       if(idonl.eq.1) then
        write(*,*)'Before'
        read(5,*)numlin
	write(*,*)'numlin',numlin
        do i=1,numlin
         read(5,*) c(i)
        end do
       endif
       call fct(40)
       call rdmom
       if(R_0.ne.0.d0) then
c read the overall asymptotics data...
        
        call rdcoef
        call rddamp
       endif
       return
      endif
c
c Compute rotated and shifted Cartesians...
c
        call ang2cart(R,oa,ob)
c        write(6,*)' Rotation matrix for A:'
c        do irb=1,3
c         write(6,*)(tranA(irb,jrb),jrb=1,3)
c        end do
c        write(6,*)' Rotation matrix for B:'
c        do irb=1,3
c         write(6,*)(tranB(irb,jrb),jrb=1,3)
c        end do
c         write(6,*)'Rotated Cartesians:'
c         write(6,*)'Monomer A:'
c         do i=1,nsitea
c          write(6,*)(siteat(k,i),k=1,3)
c         end do
c         write(6,*)'Monomer B:'
c         do i=1,nsiteb
c          write(6,*)(sitebt(k,i),k=1,3)
c         end do
c
c Calculate just the Janzen model part and exit...
c
c        call janzen(e1pol,e1exch)
c        val = e1exch
c        return
c
c Zero out the used-type matrix
c
        do k1=1,ntypemax
        do k2=1,ntypemax
         itypu(k1,k2) = 0
        end do
        end do
c
c Loop over sites in A and B
c
      iii = 1
      val = 0.d0
c initialize torques and the force of A exerted on B
      do kk=1,3
       torquea(kk) = 0.d0
       torqueb(kk) = 0.d0
       forceab(kk) = 0.d0
      end do
c Add the FQ induction part (force and torques updated too)
c      call polder(R,siteat,sitebt,fcind,forceab,torquea,torqueb)
c      call polderit(siteat,sitebt,fcind)
c      val = val + fcind
c
      do ia=1,nsitea
      do ib=1,nsiteb
c initialize the derivative (separately for each pair)
       der = 0.d0
       valp = 0.d0
c
c Compute the distance between the two sites
c
       rij = dist(siteat(1,ia),sitebt(1,ib))
c       write(6,*)
c       write(6,*)'Site of A',ia,' Site of B',ib
c       write(6,*)'Distance (Anstr):',rij
c
c Compute the "basis functions" values pertaining to the
c required for of the potential. numt is the number of terms returned...
c
c       call potparts(ntpot,rij,ia,ib,numt,values)
       call potparts(ntpot,rij,ia,ib,numt,values,valder)
c      write(*,*)ntpot,rij,numt,values(numt)
c       write(6,*)'dist,En:',ia,ib,rij/bohr2a,values(1)*1.59360d-3
c
c If abs(ntpot)>10, the last values is the fixed part.
c
       numt0 = numt
       if(abs(ntpot).gt.10) then
        numt0 = numt-1
        valp = valp + values(numt)
        der = der + valder(numt)
       endif
c
c decide which basis function it is (basis function is determined
c by the pair of site types and the position in values vector).
c iii is the first basis function for this pair of types. It is
c updated into the next available position for the next pair of types.
c
       itu = itypu(itypea(ia),itypeb(ib))
       if(itu.eq.0) then
        itypu(itypea(ia),itypeb(ib)) = iii
        itypu(itypeb(ib),itypea(ia)) = iii
        itu = iii
        iii = iii + numt0
       endif
       do i = 1,numt0
        itu1 = itu+i-1
        valp = valp + c(itu1)*values(i)
        der = der + c(itu1)*valder(i)  
c        write(6,*)'A_ab term:',values(i)*c(itu1)
       end do
c       write(6,*)'Total site-site contribution:',valp
       val = val + valp
c       write(*,*)rij,itypea(ia),itypeb(ib),ia,ib,itu,numt0,val
       
       
       
       
c der contains the derivative of v_ab over r_ab
c Calculate the contributions to the torques.
       call caltorque(torquea,torqueb,forceab,der,ia,ib,rij,R)
      end do
      end do
      
      
c      call ind_iter(siteat,nsitea,sitebt,nsiteb,R,energy)
      
c      write(6,*)'Total energy',val,energy,val+energy
      
c      val=val+energy
      
c transform the torques to the molecular frame
      call matvec2(tranA,torquea)
      call matvec2(tranB,torqueb)
c print out the torques
c      write(6,*)'TorqueA:',torquea
c      write(6,*)'TorqueB:',torqueb
c      write(6,*)'ForceAB:',forceab
c      fab = 0.d0
c      do i=1,3
c       fab = fab + forceab(i)**2
c      end do
c      fab = dsqrt(fab)
c      write(6,*)'Forceab: length:',fab
c
c Don't forget to add asymptotics when the fitting part tested...
c RB test: add the dipole-dipole induction term....
c Note that the derivatives of this are NOT included....
c       call dipind(siteat(1,1),sitebt(1,1),oa,ob,0,0,fcind)
c       val = val + fcind
      return
      end
c
c maxp - maximum nuber of data points 
c nsitemax - maximum number of sites in one monomer
c ntypemax - maximum number of site types
c maxb - maximum nuber of "basis functions", must be >=
c        0.5*nsitemax*(nsitemax+1)*number_of_terms_in_the_pair_potential
c maxpar1 - maximum number of types of one-site-type nonl parameters
c maxpar2 - maximum number of types of two-site-type nonl parameters
c
c nopar - number of optimized one-site-type nonlinear parameters
c noparab - number of optimized two-site-type nonlintypesnear parameters
c
c nparm - number of nonzero one-site-type nonlinear parameters
c nparab - number of nonzero two-site-type nonlintypesnear parameters
c
      subroutine data
      implicit real*8 (a-h,o-z)
c
      parameter (maxb=500,maxp=140000)
      parameter (nsitemax=42,ntypemax=25,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=3,maxpar2=13)
      parameter (mxl=100)
      parameter (maxn=6*nsitemax+maxpar1*ntypemax+maxpar2*ntypemax2)
c
      character*8 label(9)
      common/npar/param(maxpar1,ntypemax),
     1 parab(maxpar2,ntypemax,ntypemax),nparm,nparab
      common/sites/ sitea(3,nsitemax), siteb(3,nsitemax),nsitea,nsiteb
      common/types/ itypea(nsitemax),itypeb(nsitemax)
      common/optdrv/ iosita(2,3*nsitemax),iositb(2,3*nsitemax),
     1 iopar(2,ntypemax),ioparab(3,ntypemax2),nosa,nosb,nopar,noparab
c
      common/fit/rpt(6,maxp),dat0(maxp),sigma(maxp),iweight,ntpot,
     1 idonl, iopt
c
      common/leastsq/dat(maxp),a(maxp,maxb),a0(maxp),g(maxp),c(maxb),
     1 ata(maxb*maxb),b(maxb),works(2*maxb),iworks(maxb),np,nb
      COMMON/EOPT/TOLF,TOLR,ANOISE,TIMLIM,LINP,LOUT,ALFST,IPR,NPAR,
     &            SAFETL
      common/inou/ icon(80)
      common/grdfit/ rlow(maxn),rhigh(maxn),npts(maxn)
      common/damp/ ldum(5,mxl),cdum(maxb),ndum
      data bohr2a /0.529177249d0/
c
c Zero out the A and B site position optimization indicators...
c
      do i=1,3*nsitemax
      do j=1,2
       iosita(j,i) = 0
       iositb(j,i) = 0
      end do
      end do 
c
c read in the sites of monomer A
c
      nosa = 0
      read(5,*) nsitea
      do i=1,nsitea
       read(5,*) itypea(i),(sitea(j,i),j=1,3),k1,k2,k3
c transform to Angstroms...
       do j=1,3
        sitea(j,i) = bohr2a*sitea(j,i)
       end do
       if(k1.ne.0) then
        nosa = nosa + 1
        iosita(1,nosa) = i
        iosita(2,nosa) = 1
       endif
       if(k2.ne.0) then
        nosa = nosa + 1
        iosita(1,nosa) = i
        iosita(2,nosa) = 2
       endif
       if(k3.ne.0) then
        nosa = nosa + 1
        iosita(1,nosa) = i
        iosita(2,nosa) = 3
       endif
      end do
c
c read in the sites of monomer B
c
      nosb = 0
      read(5,*) nsiteb
      do i=1,nsiteb
       read(5,*) itypeb(i),(siteb(j,i),j=1,3),k1,k2,k3
c transform to Angstroms...
       do j=1,3
        siteb(j,i) = bohr2a*siteb(j,i)
       end do
       if(k1.ne.0) then
        nosb = nosb + 1
        iositb(1,nosb) = i
        iositb(2,nosb) = 1
       endif
       if(k2.ne.0) then
        nosb = nosb + 1
        iositb(1,nosb) = i
        iositb(2,nosb) = 2
       endif
       if(k3.ne.0) then
        nosb = nosb + 1
        iositb(1,nosb) = i
        iositb(2,nosb) = 3
       endif
      end do
c
c Zero out the nonlinear parameters matrices
c
      do i=1,ntypemax
      do k=1,maxpar1
       param(k,i) = 0.d0
      end do
      end do
      do i=1,ntypemax
      do j=1,ntypemax
      do k=1,maxpar2
       parab(k,i,j) = 0.d0
      end do
      end do
      end do
c
c Zero out the nonl params optimization indicators...
c
      do i=1,ntypemax
      do k=1,2
       iopar(k,i) = 0
      end do
      end do
      do i=1,ntypemax2
      do k=1,3
       ioparab(k,i) = 0
      end do
      end do
c
c Read in the one-site nonlinear parameters and opt indicators..
c
      nopar = 0
      read(5,*) nparm
      do i=1,nparm
       read(5,*) ityp, inumpar, val, indopt
       param(inumpar,ityp) = val
c Convert charges to the proper units 
       if(inumpar.eq.1) param(1,ityp) = 18.22262373d0*param(1,ityp)
       if(indopt.ne.0) then
        nopar = nopar + 1
        iopar(1,nopar) = ityp
        iopar(2,nopar) = inumpar
       endif
      end do
c
c Read in the two-site nonlinear parameters and opt indicators..
c
      noparab = 0
      read(5,*) nparab
      do i=1,nparab
       read(5,*) ityp1,ityp2,inumpar, val, indopt
       parab(inumpar,ityp1,ityp2) = val  
       parab(inumpar,ityp2,ityp1) = val
       if(indopt.ne.0) then
        noparab = noparab + 1
        ioparab(1,noparab) = ityp1
        ioparab(2,noparab) = ityp2  
        ioparab(3,noparab) = inumpar
       endif
      end do
c
c Read the type of the potential form to be used, inear coeff. reading
c indicator, and the optimization parameters....
c
      read(5,*) ntpot, idonl, iopt, iweight
      read(5,*) TOLF,TOLR,ANOISE,TIMLIM,LINP,LOUT,ALFST,IPR,
     1            SAFETL,icon(2)
      if(iopt.eq.2) then
c
c read in the grid upper boundaries and mesh...
c
       read(5,*) ndum
       do i=1,ndum
        read(5,*)rhigh(i),npts(i)
       end do
      endif
c
      entry output_old 
c
c Print the read-in parameters....
c
      write(6,*)'            Input parameters'
c
      write(6,*)'Monomer A has',nsitea,' sites:'
      write(6,*)'type         x         y          z'
      do i=1,nsitea
       write(6,10) itypea(i), (sitea(j,i),j=1,3)
      end do
 10   format(2x,i4,3e20.7)
      write(6,*)nosa,' opitmized coordinates for monomer A:'
      do i=1,nosa
       write(6,*)'for site',iosita(1,i),' coordinate',iosita(2,i)
      end do
c
      write(6,*)'Monomer B has',nsiteb,' sites:'
      write(6,*)'type         x         y          z'
      do i=1,nsiteb
       write(6,10) itypeb(i), (siteb(j,i),j=1,3)
      end do
      write(6,*)nosb,' opitmized coordinates for monomer B:'
      do i=1,nosb
       write(6,*)'for site',iositb(1,i),' coordinate',iositb(2,i)
      end do
c
      write(6,*)'One-site nonlinear parameters:'
      do k=1,maxpar1
       write(6,*)'Parameter',k
       do i=1,ntypemax
c Convert charges back to atomic units.. BE CAREFUL, !!!!!
        if(k.eq.1) param(k,i) = param(k,i)/18.22262373d0
        indopt = 0
        do l=1,nopar
         if((iopar(1,l).eq.i).and.(iopar(2,l).eq.k)) indopt = 1
        end do
        if(param(k,i).ne.0.d0.or.indopt.eq.1) then
        if(indopt.eq.0) then
         write(6,*)'Type',i,param(k,i),' fixed'
        else
         write(6,*)'Type',i,param(k,i),' optimized'
        endif
        endif
       end do
      end do
c
      write(6,*)'Two-site nonlinear parameters:'
      do k=1,maxpar2
       write(6,*)'Parameter',k
       do i=1,ntypemax
       do j=1,i
        indopt = 0
        do l=1,noparab
         if((ioparab(1,l).eq.i).and.(ioparab(2,l).eq.j).
     1 and.(ioparab(3,l).eq.k)) indopt = 1
         if((ioparab(2,l).eq.i).and.(ioparab(1,l).eq.j).
     1 and.(ioparab(3,l).eq.k)) indopt = 1 
        end do
        if(parab(k,i,j).ne.0.d0.or.indopt.eq.1) then
        if(indopt.eq.0) then
         write(6,*)'Types',i,j,parab(k,i,j),' fixed'
        else
         write(6,*)'Types',i,j,parab(k,i,j),' optimized'
        endif
        endif
       end do
       end do
       end do
c
       write(6,*)'Potential type used:',ntpot
       write(6,*)'damping parameters'
       do i=1,ndum
        write(6,13)(ldum(k,i),k=1,5),cdum(i)
       end do
c
 13     format(5i4,e20.10)
       return
       end
c
c For a given set of Euler angles (in radians) and the intermonomer separation R
c compute Cartesian coordinates of all sites. Assume that the initial
c coordinates of the monomer X sites are given in the center of mass system 
c of monomer X. Output in common/temp_sites/
c
      subroutine ang2cart(R,oa,ob)
      implicit real*8 (a-h,o-z)
c
      parameter (maxb=500,maxp=140000)
      parameter (nsitemax=42,ntypemax=25,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=3,maxpar2=13)
      parameter (mxl=100)
c
      dimension trand(3,3), tranc(3,3),oa(3),ob(3)
c
      common/sites/ sitea(3,nsitemax), siteb(3,nsitemax),nsitea,nsiteb
c
      common/temp_sites/ siteat(3,nsitemax), sitebt(3,nsitemax),
     1                   trand,tranc
c
c
c Build the rotation matrix for monomer A...
c
c      ad = 0.d0
      ad = oa(1)
      bd = oa(2)
      gd = oa(3)
      sad = dsin(ad)
      sbd = dsin(bd)
      sgd = dsin(gd)
      cad = dcos(ad)
      cbd = dcos(bd)
      cgd = dcos(gd)
      trand(1,1) = cad*cbd*cgd - sad*sgd
      trand(1,2) = -cad*cbd*sgd - sad*cgd
      trand(1,3) = cad*sbd
      trand(2,1) = sad*cbd*cgd + cad*sgd
      trand(2,2) = -sad*cbd*sgd + cad*cgd
      trand(2,3) = sad*sbd
      trand(3,1) = -sbd*cgd
      trand(3,2) = sbd*sgd
      trand(3,3) = cbd
c
c Build the rotation matrix for monomer B...
c
c      ac = ob(1)-oa(1)
      ac = ob(1)
      bc = ob(2)
      gc = ob(3)
      sac = dsin(ac)
      sbc = dsin(bc)
      sgc = dsin(gc)
      cac = dcos(ac)
      cbc = dcos(bc)
      cgc = dcos(gc)
      tranc(1,1) = cac*cbc*cgc - sac*sgc
      tranc(1,2) = -cac*cbc*sgc - sac*cgc
      tranc(1,3) = cac*sbc
      tranc(2,1) = sac*cbc*cgc + cac*sgc
      tranc(2,2) = -sac*cbc*sgc + cac*cgc
      tranc(2,3) = sac*sbc
      tranc(3,1) = -sbc*cgc
      tranc(3,2) = sbc*sgc
      tranc(3,3) = cbc
c
c Transform the coordinates of monomer A sites...
c
      do i=1,nsitea
       call matvec1(trand,sitea(1,i),siteat(1,i))
      end do
c
c Transform the coordinates of monomer B sites...
c
      do i=1,nsiteb
       call matvec1(tranc,siteb(1,i),sitebt(1,i))
      end do
c
c Shift monomer B by R in the positive z direction...
c
      do i=1,nsiteb
       sitebt(3,i) = sitebt(3,i) + R
      end do
c
      return
      end
c
c ----- Multiply vetor v by the matrix a. Store in u.
c
       subroutine matvec1(a,v,u)
       implicit real*8 (a-h,o-z)
       dimension a(3,3),v(3),u(3)
c
       n = 3
       do i=1,n
        u(i) = 0.d0
        do j=1,n
         u(i) = u(i) + a(i,j)*v(j)
        end do
       end do
       return
       end
c
c ----- Multiply vetor v by the transposed matrix a. Store in v again.
c
       subroutine matvec2(a,v)
       implicit real*8 (a-h,o-z)
       dimension a(3,3),v(3),u(3)
c
       n = 3
       do i=1,n
        u(i) = 0.d0
        do j=1,n
         u(i) = u(i) + a(j,i)*v(j)
        end do
       end do
       do i=1,n
        v(i) = u(i)
       end do
       return
       end
c
      function dist(a,b)
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3)
      ttt = 0.d0
      do i=1,3
       diff = a(i)-b(i)
       ttt = ttt + diff*diff
      end do
      dist = dsqrt(ttt)
      return
      end
c
c      subroutine potparts(ntpot,rij,ia,ib,numt,values)
      subroutine potparts(ntpot,rij,ia,ib,numt,values,valder)
      implicit real*8 (a-h,o-z)
c
      parameter (maxb=500,maxp=140000)
      parameter (nsitemax=42,ntypemax=25,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=3,maxpar2=13)
      parameter (mxl=100)
c
      dimension values(1),valder(1)
      common/npar/param(maxpar1,ntypemax),
     1 parab(maxpar2,ntypemax,ntypemax),nparm,nparab
      common/types/ itypea(nsitemax),itypeb(nsitemax)
      common/inout_pot/ inoutp
      common/radder/ der
      common/misc/ R_0,isyst,npowers 
      data a0 /0.529177249d0/
      data au2cal /627.510d0/
c*****************************************************
c ntpot =11 : LJ + charge-electrostatics (like Murthy 5q potential)
c*****************************************************
      if(ntpot.eq.11) then
c Parameters:
c      epsilon ---- linear parameter...
      r_0 = parab(1,itypea(ia),itypeb(ib)) 
      qa = param(1,itypea(ia))    
      qb = param(1,itypeb(ib))    
c
      numt = 0
      if(r_0.gt.0.d0) then
       numt = numt + 1
       r6 = (rij/r_0)**(-6)
       r12 = r6*r6
       values(numt) = (r12 - 2*r6)
      endif
      numt = numt + 1
      values(numt) = qa*qb/rij
      return
      endif
c ********
c Same as 11 but totally nonlienar....
c ********
      if(ntpot.eq.12) then
c Parameters:
      eps = parab(2,itypea(ia),itypeb(ib))
      r_0 = parab(1,itypea(ia),itypeb(ib))
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      values(1) = 0.d0
      if(r_0.gt.0.d0) then
       r6 = (r_0/rij)**6
       r12 = r6*r6
c       values(1) = values(1) + eps*(r12 - 2*r6)
       values(1) = values(1) + eps*(r12 - r6)
      endif
      values(1) = values(1) + qa*qb/rij
      return
      endif
c ********
c ntpot=19: LJ + elst, different formula...
c ********
      if(ntpot.eq.19) then
c Parameters:
      eps = parab(2,itypea(ia),itypeb(ib))
c convert eps from kJ/mol into kcal/mol (consistently with the charges.
      eps = 0.239006d0*eps
      r_0 = parab(1,itypea(ia),itypeb(ib))
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      values(1) = 0.d0
      if(r_0.gt.0.d0) then
       r6 = (r_0/rij)**6
       r12 = r6*r6
       values(1) = values(1) + 4.d0*eps*(r12 - r6)
      endif
      values(1) = values(1) + qa*qb/rij
      return
      endif
c ***********************************************************
c  ntpot = 20 : exponential + r^-6 + ... + r^-12 + elst, 
c  Swiss potential
c ***********************************************************
      if(ntpot.eq.20) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = parab(2,itypea(ia),itypeb(ib))
      c12 = parab(3,itypea(ia),itypeb(ib))
      c10 = parab(4,itypea(ia),itypeb(ib))
      c8 = parab(5,itypea(ia),itypeb(ib))
      c6 = parab(6,itypea(ia),itypeb(ib))
      c6 = c6*au2cal*a0**6
      c8 = c8*au2cal*a0**8
      c10 = c10*au2cal*a0**10
      c12 = c12*au2cal*a0**12
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      values(1) = 0.d0
      if(beta.gt.0.d0) then
       values(1) = values(1) + a*dexp(-beta*rij)
      endif
      r2 = rij**(-2)
      r6 = r2*r2*r2
      r8 = r6*r2
      r10 = r8*r2
      r12 = r10*r2
      dmp = 1.d0+dexp(-2.d0*(rij/a0 - 2.1d0))
      dmp = 1.d0/dmp**8
      asym = dmp*(c6*r6 + c8*r8 + c10*r10 + c12*r12 + qa*qb/rij)
      values(1) = values(1) + asym
      return
      endif
c ***********************************************************
c  ntpot = 17 : exponential + r^-6,7,8 + elst, different damping
c ***********************************************************
      if(ntpot.eq.17) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c7 = parab(4,itypea(ia),itypeb(ib))
      c8 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
      dmp6 = parab(7,itypea(ia),itypeb(ib))
      dmp7 = parab(8,itypea(ia),itypeb(ib))
      dmp8 = parab(9,itypea(ia),itypeb(ib))
c      dmp1 = parab(6,itypea(1),itypeb(1))
c      dmp6 = parab(7,itypea(1),itypeb(1))
c      dmp7 = parab(8,itypea(1),itypeb(1))
c      dmp8 = parab(9,itypea(1),itypeb(1))
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      values(1) = 0.d0
      if(beta.gt.0.d0) then
       values(1) = values(1) + a*dexp(-beta*rij)
      endif
c        x = 0.2d0*rij   
c        CALL ERROR(X, ERRORF)
       d1 = d(8,dmp1,rij)
       d6 = d(8,dmp6,rij)
       d7 = d(8,dmp7,rij)
       d8 = d(8,dmp8,rij)
       p_els = d1*qa*qb/rij
       p_6   = -d6*c6/(rij)**6
       p_7   = -d7*c7/rij**7
       p_8   = -d8*c8/rij**8
c       ee0 = -500.d0
       ee0 = 100.d0
       go to 222
c       if(p_els.lt.ee0) then
       if(d1*qa*qb.lt.ee0) then
c        write(6,*)'damelst',ia,ib
        inoutp = 7
       endif
c       if(p_6.lt.ee0) then
       if(-d6*c6.lt.ee0) then
c        write(6,*)'dam6',ia,ib
        inoutp = 7
       endif
c       if(p_7.lt.ee0) then
       if(-d7*c7.lt.ee0) then
c        write(6,*)'dam7',ia,ib
        inoutp = 7
       endif
c       if(p_8.lt.ee0) then
       if(-d8*c8.lt.ee0) then
c        write(6,*)'dam8',ia,ib
        inoutp = 7
       endif
c 222  continue
       r_1 = 1.d0
       d10 = d(1,dmp1,r_1)
       d60 = d(6,dmp6,r_1)
       d70 = d(7,dmp7,r_1)
       d80 = d(8,dmp8,r_1)
      chk1a = a*dexp(-beta) + d10*qa*qb - d60*c6 - d70*c7 - d80*c8
c      write(6,*)'chk1a',chk1a
      if(chk1a.lt.ee0) then
       inoutp = inoutp + 1
c       write(6,*)'Problems with pair',ia,ib,inoutp,chk1a    
      endif
 222  continue
      values(1) = values(1) + p_els + p_6 + p_7 + p_8
      return
      endif
c ***********************************************************
c  ntpot = 15 : exponential + r^-6,8,10 + elst, different damping
c ***********************************************************
      if(ntpot.eq.15) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(5,itypea(ia),itypeb(ib))
      c10 = parab(8,itypea(ia),itypeb(ib))
      dmp1 = parab(4,itypea(ia),itypeb(ib))
      dmp6 = parab(6,itypea(ia),itypeb(ib))
      dmp8 = parab(7,itypea(ia),itypeb(ib))
      dmp10 = parab(9,itypea(ia),itypeb(ib))
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      values(1) = 0.d0
      if(beta.gt.0.d0) then
       values(1) = values(1) + a*dexp(-beta*rij)
      endif
c        x = 0.2d0*rij   
c        CALL ERROR(X, ERRORF)
       d1 = d(1,dmp1,rij)
       d6 = d(6,dmp6,rij)
       d8 = d(8,dmp8,rij)
       d10 = d(10,dmp10,rij)
      values(1) = values(1) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
     1 d8*c8/rij**8 - d10*c10/rij**10
      go to 111
      write(6,*)'atoms',ia,ib,itypea(ia),itypeb(ib)
      write(6,*)'rij',rij
      write(6,*)'dmp',dmp1,dmp6,dmp8
      write(6,*)'d',d1,d6,d8
      write(6,*)'C  ',c6,c8
      write(6,*)'a,beta',a,beta
      write(6,*)'values(1)',values(1)
 111  continue
      return
      endif
c ***********************************************************
c  ntpot = 162 : exponential*(1/r+sum r^n) + r^-(6-8-10) + elst, 
c partially nonlinear..
c ***********************************************************
      if(ntpot.eq.162) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
c      dmp1 = parab(6,1,1)
      dmp6 = parab(7,itypea(ia),itypeb(ib))
c      dmp6 = parab(7,1,1)
      dmp8 = parab(8,itypea(ia),itypeb(ib)) 
c      dmp8 = parab(8,1,1)
      dmp10= parab(9,itypea(ia),itypeb(ib))
c      dmp10 = parab(9,1,1)
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      if(beta.gt.0.d0) numt = npowers + 1
      values(numt) = 0.d0
      if(beta.gt.0.d0) then
       values(numt) = values(numt) + a*dexp(-beta*rij)/rij
      endif
       d1 = d(1,dmp1,rij)
       d6 = d(6,dmp6,rij)
       d8 = d(8,dmp8,rij)
       d10 = d(10,dmp10,rij)
      values(numt) = values(numt) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
     1 d8*c8/rij**8 - d10*c10/rij**10
       if(numt.ne.1) then
        values(1) = a*dexp(-beta*rij)
        do kk=2,npowers
         values(kk) = values(kk-1)*rij
        end do
       endif
      return
      endif
c ***********************************************************
c  ntpot = 163 : exponential*(1/r^3+sum 1/r^n) + r^-(6-8-10) + elst, partially  
c  nonlinear..
c ***********************************************************
      if(ntpot.eq.163) then
       dlambda = 1.d0
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))/dlambda
      a = dexp(parab(2,itypea(ia),itypeb(ib)))*dlambda
      c6 = parab(3,itypea(ia),itypeb(ib))*dlambda**7
      c8 = parab(4,itypea(ia),itypeb(ib))*dlambda**9
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))*dlambda
c      dmp1 = parab(6,1,1)
      dmp6 = parab(7,itypea(ia),itypeb(ib))*dlambda
c      dmp6 = parab(7,1,1)
      dmp8 = parab(8,itypea(ia),itypeb(ib))*dlambda
c      dmp8 = parab(8,1,1)
      dmp10= parab(9,itypea(ia),itypeb(ib))
c make dmp10 safe for potentials with no c10 component
      if(dmp10.eq.0.d0) dmp10 = 1.d0
c      dmp10 = parab(9,1,1)
c      qa = param(1,itypea(ia))*dlambda*2
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      if(beta.gt.0.d0) numt = npowers + 1
      values(numt) = 0.d0
      if(beta.gt.0.d0) then
       values(numt) = values(numt) + a*dexp(-beta*rij)/rij**npowers
      endif
c      write(6,*)'kappa term:',values(numt)
       d1 = d(1,dmp1,rij)
       d6 = d(6,dmp6,rij)
       d8 = d(8,dmp8,rij)
       d10 = d(10,dmp10,rij)
c       write(6,*)'Damping functions:'
c       write(6,*)'d1:',d1
c       write(6,*)'d6:',d6
c       write(6,*)'d8:',d8
c       write(6,*)'d10:',d10
      values(numt) = values(numt) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
     1 d8*c8/rij**8 - d10*c10/rij**10
c       write(6,*)'Damped asymptotics:'
c       write(6,*)'Elst:',d1*qa*qb/rij
c       write(6,*)'R^-6:',- d6*c6/(rij)**6
c       write(6,*)'R^-8:',- d8*c8/(rij)**8
c       write(6,*)'R^-10:',- d10*c10/(rij)**10
       if(numt.ne.1) then
        values(1) = a*dexp(-beta*rij)
        do kk=2,npowers
         values(kk) = values(kk-1)/rij
        end do
       endif
      return
      endif
c ***********************************************************
c  ntpot = 164 : exponential*(sum r^n) + r^-(6-8-10) + elst, partially  
c  nonlinear..
c ***********************************************************
      if(ntpot.eq.164) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
c      dmp1 = parab(6,1,1)
      dmp6 = parab(7,itypea(ia),itypeb(ib))
c      dmp6 = parab(7,1,1)
      dmp8 = parab(8,itypea(ia),itypeb(ib)) 
c      dmp8 = parab(8,1,1)
      dmp10= parab(9,itypea(ia),itypeb(ib))
c      dmp10 = parab(9,1,1)
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      if(beta.gt.0.d0) numt = npowers + 1
      values(numt) = 0.d0
       d1 = d(1,dmp1,rij)
       d6 = d(6,dmp6,rij)
       d8 = d(8,dmp8,rij)
       d10 = d(10,dmp10,rij)
      values(numt) = values(numt) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
     1 d8*c8/rij**8 - d10*c10/rij**10
       if(numt.ne.1) then
        values(1) = a*dexp(-beta*rij)
        do kk=2,npowers
         values(kk) = rij*values(kk-1)
        end do
       endif
      return
      endif
c ***********************************************************
c  ntpot = 161 : exponential*(1+sum r^n) + r^-(6-8-10) + elst, partially  
c  nonlinear..
c ***********************************************************
      if(ntpot.eq.161) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
c      a = parab(2,itypea(ia),itypeb(ib))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
c      dmp1 = parab(6,7,1)
      dmp6 = parab(7,itypea(ia),itypeb(ib))
c      dmp6 = parab(7,7,1)
      dmp8 = parab(8,itypea(ia),itypeb(ib)) 
c      dmp8 = parab(7,7,1)
      dmp10= parab(9,itypea(ia),itypeb(ib))
c      dmp10 = parab(7,7,1)
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      if(beta.gt.0.d0) numt = npowers + 1
      values(numt) = 0.d0
      if(beta.gt.0.d0) then
       values(numt) = values(numt) + a*dexp(-beta*rij)
      endif
       d1 = d(1,dmp1,rij)
       d6 = d(6,dmp6,rij)
       d8 = d(8,dmp8,rij)
       d10 = d(10,dmp10,rij)
c        d6=1.0d0
c	d8=1.0d0
c	d10=1.0d0
      values(numt) = values(numt) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
     1 d8*c8/rij**8 - d10*c10/rij**10
c       write(6,*)'C6  term:',-418.400d0*d6*c6/(rij)**6
c       write(6,*)'C8  term:',-418.400d0*d8*c8/rij**8
c       write(6,*)'C10 term:',-418.400d0*d10*c10/rij**10
c       write(6,*)'Coul trm:',d1*qa*qb/rij
       if(numt.ne.1) then
        values(1) = rij*a*dexp(-beta*rij)
        do kk=2,npowers
         values(kk) = rij*values(kk-1)
        end do
       endif
c       go to 333
c Calculate the derivative of this a-b term over the distance
       der = -d1*qa*qb/(rij*rij) + 6.d0*d6*c6/(rij)**7 
       der = der + 8.d0*d8*c8/rij**9 + 10.d0*d10*c10/rij**11
       dd1 = dd(1,dmp1,rij)
       dd6 = dd(6,dmp6,rij)
       dd8 = dd(8,dmp8,rij)
       dd10 = dd(10,dmp10,rij)
       der = der + dd1*qa*qb/rij - dd6*c6/(rij)**6 -
     1             dd8*c8/rij**8 - dd10*c10/rij**10
       valder(numt) = der
       if(beta.ne.0.d0) then
        valder(numt) = valder(numt) - beta*a*dexp(-beta*rij)
       if(numt.ne.1) then
        do kk=1,npowers
         valder(kk) = (dfloat(kk)/rij - beta)*values(kk)
        end do
       endif
       endif
 333  continue
      return
      endif


c ***********************************************************
c  ntpot = 262 : exponential*(1+sum r^n) + r^-(6-8-10) + elst, partially  
c  nonlinear.. only one dampling factor
c ***********************************************************
      if(ntpot.eq.262.or.ntpot.eq.362) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
c      a = parab(2,itypea(ia),itypeb(ib))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
c      dmp1 = parab(6,7,1)
       dmp6 = parab(7,itypea(ia),itypeb(ib))
c      dmp6 = parab(7,7,1)
c      dmp8 = parab(8,itypea(ia),itypeb(ib)) 
c      dmp8 = parab(7,7,1)
c      dmp10= parab(9,itypea(ia),itypeb(ib))
c      dmp10 = parab(7,7,1)
      dmp8=dmp6
      dmp10=dmp6
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      if(beta.gt.0.d0) numt = npowers + 1
      values(numt) = 0.d0
      if(beta.gt.0.d0) then
       values(numt) = values(numt) + a*dexp(-beta*rij)
      endif
       d1 = d(1,dmp1,rij)
       d6 = d(6,dmp6,rij)
       d8 = d(8,dmp8,rij)
       d10 = d(10,dmp10,rij)
c        d6=1.0d0
c	d8=1.0d0
c	d10=1.0d0
      values(numt) = values(numt) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
     1 d8*c8/rij**8 - d10*c10/rij**10
c       write(6,*)'C6  term:',-418.400d0*d6*c6/(rij)**6
c       write(6,*)'C8  term:',-418.400d0*d8*c8/rij**8
c       write(6,*)'C10 term:',-418.400d0*d10*c10/rij**10
c       write(6,*)'Coul trm:',d1*qa*qb/rij
       if(numt.ne.1) then
        values(1) = rij*a*dexp(-beta*rij)
        do kk=2,npowers
         values(kk) = rij*values(kk-1)
        end do
       endif
      end if


c ***********************************************************
c  ntpot = 362 : exponential*(1+sum r^n) + r^-(6-8-10) + elst, partially  
c  nonlinear.. damping factors contrained
c ***********************************************************
c      if(ntpot.eq.362) then
cc Parameters:
c      beta = parab(1,itypea(ia),itypeb(ib))
c      a = dexp(parab(2,itypea(ia),itypeb(ib)))
cc      a = parab(2,itypea(ia),itypeb(ib))
c      c6 = parab(3,itypea(ia),itypeb(ib))
c      c8 = parab(4,itypea(ia),itypeb(ib))
c      c10 = parab(5,itypea(ia),itypeb(ib))
c      dmp1 = parab(6,itypea(ia),itypeb(ib))
      
c       dmp6 = parab(7,itypea(ia),itypeb(ib))

c      IF (dmp1.lt.0.5d0) dmp1=0.5d0
c      IF (dmp6.lt.0.5d0) dmp6=0.5d0
c      IF (beta.lt.0.5d0) beta=0.5d0

c      dmp8=dmp6
c      dmp10=dmp6
c      qa = param(1,itypea(ia))
c      qb = param(1,itypeb(ib))
c
      
      
      
c      numt = 1
c      if(beta.gt.0.d0) numt = npowers + 1
c      values(numt) = 0.d0
c      if(beta.gt.0.d0) then
c       values(numt) = values(numt) + a*dexp(-beta*rij)
c      endif
c       d1 = d(1,dmp1,rij)
c       d6 = d(6,dmp6,rij)
c       d8 = d(8,dmp8,rij)
c       d10 = d(10,dmp10,rij)
c      values(numt) = values(numt) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
c     1 d8*c8/rij**8 - d10*c10/rij**10
c       if(numt.ne.1) then
c        values(1) = rij*a*dexp(-beta*rij)
c        do kk=2,npowers
c         values(kk) = rij*values(kk-1)
c        end do
c       endif
c      end if



c ***********************************************************
c  ntpot = 2622 : exponential*(1+sum r^n) + r^-(6-8-10) + elst, partially  
c  nonlinear.. only one dampling factor for all
c ***********************************************************
      if(ntpot.eq.2622) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
c      a = parab(2,itypea(ia),itypeb(ib))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = param(2,1)
c      dmp1 = parab(6,7,1)
       dmp6 = parab(7,itypea(ia),itypeb(ib))
c      dmp6 = parab(7,7,1)
c      dmp8 = parab(8,itypea(ia),itypeb(ib)) 
c      dmp8 = parab(7,7,1)
c      dmp10= parab(9,itypea(ia),itypeb(ib))
c      dmp10 = parab(7,7,1)
      dmp8=dmp6
      dmp10=dmp6
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c    
c      write(*,*)"dmp",dmp1
      numt = 1
      if(beta.gt.0.d0) numt = npowers + 1
      values(numt) = 0.d0
      if(beta.gt.0.d0) then
       values(numt) = values(numt) + a*dexp(-beta*rij)
      endif
       d1 = d(1,dmp1,rij)
       d6 = d(6,dmp6,rij)
       d8 = d(8,dmp8,rij)
       d10 = d(10,dmp10,rij)
c        d6=1.0d0
c	d8=1.0d0
c	d10=1.0d0
      values(numt) = values(numt) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
     1 d8*c8/rij**8 - d10*c10/rij**10
c       write(6,*)'C6  term:',-418.400d0*d6*c6/(rij)**6
c       write(6,*)'C8  term:',-418.400d0*d8*c8/rij**8
c       write(6,*)'C10 term:',-418.400d0*d10*c10/rij**10
c       write(6,*)'Coul trm:',d1*qa*qb/rij
       if(numt.ne.1) then
        values(1) = rij*a*dexp(-beta*rij)
        do kk=2,npowers
         values(kk) = rij*values(kk-1)
        end do
       endif
      end if


c ***********************************************************
c  ntpot = 2626 : Only c6,c8,c10 non linear. Fit squareroot
c ********* **************************************************
      if(ntpot.eq.2626) then
c Parameters:
      c6sq = parab(3,itypea(ia),itypeb(ib))
      c8sq = parab(4,itypea(ia),itypeb(ib))
      c10sq = parab(5,itypea(ia),itypeb(ib))

c
      numt = 1
      values(1) = - c6sq*c6sq/(rij)**6 - 
     1 c8sq*c8sq/rij**8 - c10sq*c10sq/rij**10
      return
      endif


c ***********************************************************
c  ntpot = 26266 : Only c6 non linear. Combining rules
c ********* **************************************************
      if(ntpot.eq.26266) then
c Parameters:
      c6a = param(3,itypea(ia))
      c6b = param(3,itypea(ib))      
c
      c6=sqrt(c6a*c6b)
c       write(*,*) c6a,c6b,ia,ib,itypea(ia),itypeb(ia)
      numt = 1
      values(1) = - c6/(rij)**6 
      return
      endif


c ***********************************************************
c  ntpot = 177 : exponential*(1+sum r^n) + r^-(6,7,9-8-10) + elst, partially  
c  nonlinear..
c ***********************************************************
      if(ntpot.eq.177) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
c      a = parab(2,itypea(ia),itypeb(ib))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c7 = parab(4,itypea(ia),itypeb(ib))
      c8 = parab(5,itypea(ia),itypeb(ib))
      c9 = parab(6,itypea(ia),itypeb(ib))
      c10 = parab(7,itypea(ia),itypeb(ib))
      dmp1 = parab(8,itypea(ia),itypeb(ib))
c      dmp1 = parab(6,7,1)
      dmp6 = parab(9,itypea(ia),itypeb(ib))
c      dmp6 = parab(7,7,1)
      dmp7 = parab(10,itypea(ia),itypeb(ib)) 
c      dmp8 = parab(7,7,1)
      dmp8= parab(11,itypea(ia),itypeb(ib))
      dmp9= parab(12,itypea(ia),itypeb(ib))
      dmp10= parab(13,itypea(ia),itypeb(ib))
c      dmp10 = parab(7,7,1)
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      if(beta.gt.0.d0) numt = npowers + 1
      values(numt) = 0.d0
      if(beta.gt.0.d0) then
       values(numt) = values(numt) + a*dexp(-beta*rij)
      endif
       d1 = d(1,dmp1,rij)
       d6 = d(6,dmp6,rij)
       d7 = d(7,dmp7,rij)
       d8 = d(8,dmp8,rij)
       d9 = d(9,dmp9,rij)
       d10 = d(10,dmp10,rij)
c        d6=1.0d0
c	d8=1.0d0
c	d10=1.0d0
      values(numt) = values(numt) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
     1 d8*c8/rij**8 - d10*c10/rij**10 - d7*c7/rij**7 - d9*c9/rij**9
c       write(6,*)'C6  term:',-418.400d0*d6*c6/(rij)**6
c       write(6,*)'C8  term:',-418.400d0*d8*c8/rij**8
c       write(6,*)'C10 term:',-418.400d0*d10*c10/rij**10
c       write(6,*)'Coul trm:',d1*qa*qb/rij
       if(numt.ne.1) then
        values(1) = rij*a*dexp(-beta*rij)
        do kk=2,npowers
         values(kk) = rij*values(kk-1)
        end do
       endif
      return
      endif


      
      
c ***********************************************************
c  ntpot = 1161 : exponential*(1+sum r^n) + r^-(6-8-10) + elst, partially  
c  nonlinear..
c ***********************************************************
      if(ntpot.eq.1161) then
c Parameters:
c
      numt = 4
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      values(1)=0.0d0
      values(2)=0.0d0
      values(3)=0.0d0
      
      if (c6 .ne. 0.0d0) then
        values(1) = - 1.0d0/(rij)**6 
        values(2) = - 1.0d0/(rij)**8 
        values(3) = - 1.0d0/(rij)**10 
      end if
c     1 d8*c8/rij**8 - d10*c10/rij**10
c       write(6,*)'C6  term:',-418.400d0*d6*c6/(rij)**6
c       write(6,*)'C8  term:',-418.400d0*d8*c8/rij**8
c       write(6,*)'C10 term:',-418.400d0*d10*c10/rij**10
c       write(6,*)'Coul trm:',d1*qa*qb/rij
      return
      endif

c ***********************************************************
c  ntpot = 1177 : exponential*(1+sum r^n) + r^-(6,7,8,9,10) + elst, partially  
c  nonlinear..
c ***********************************************************
      if(ntpot.eq.1177) then
c Parameters:
c
      numt = 6
      c6 = parab(3,itypea(ia),itypeb(ib))
      c7 = parab(4,itypea(ia),itypeb(ib))      
      c8 = parab(5,itypea(ia),itypeb(ib))
      c9 = parab(6,itypea(ia),itypeb(ib))
      c10 = parab(7,itypea(ia),itypeb(ib))      
      values(1)=0.0d0
      values(2)=0.0d0
      values(3)=0.0d0
      values(4)=0.0d0      
      values(5)=0.0d0      
      
      if (c6 .ne. 0.0d0) then
        values(1) = - 1.0d0/(rij)**6 
        values(2) = - 1.0d0/(rij)**7
        values(3) = - 1.0d0/(rij)**8 
        values(4) = - 1.0d0/(rij)**9
        values(5) = - 1.0d0/(rij)**10
      end if
c     1 d8*c8/rij**8 - d10*c10/rij**10
c       write(6,*)'C6  term:',-418.400d0*d6*c6/(rij)**6
c       write(6,*)'C8  term:',-418.400d0*d8*c8/rij**8
c       write(6,*)'C10 term:',-418.400d0*d10*c10/rij**10
c       write(6,*)'Coul trm:',d1*qa*qb/rij
      return
      endif


      
      
c ***********************************************************
c  ntpot = 261 : exponential*(1+sum r^n) + r^-(6-8-10) + elst, partially  
c  nonlinear..
c ***********************************************************
      if(ntpot.eq.261) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
c      a = parab(2,itypea(ia),itypeb(ib))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
c      dmp1 = parab(6,7,1)
      dmp6 = parab(7,itypea(ia),itypeb(ib))
c      dmp6 = parab(7,7,1)
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      if(beta.gt.0.d0) numt = npowers + 1
      values(numt) = 0.d0
      if(beta.gt.0.d0) then
       values(numt) = values(numt) + a*dexp(-beta*rij)
      endif
       d1 = d(1,dmp1,rij)
       d6 = d(6,dmp6,rij)
c        d6=1.0d0
c	d8=1.0d0
c	d10=1.0d0
      values(numt) = values(numt) + d1*qa*qb/rij - d6*c6/(rij)**6
c       write(6,*)'C6  term:',-418.400d0*d6*c6/(rij)**6
c       write(6,*)'C8  term:',-418.400d0*d8*c8/rij**8
c       write(6,*)'C10 term:',-418.400d0*d10*c10/rij**10
c       write(6,*)'Coul trm:',d1*qa*qb/rij
       if(numt.ne.1) then
        values(1) = rij*a*dexp(-beta*rij)
        do kk=2,npowers
         values(kk) = rij*values(kk-1)
        end do
       endif
c       go to 333
c Calculate the derivative of this a-b term over the distance
       der = -d1*qa*qb/(rij*rij) + 6.d0*d6*c6/(rij)**7 
       der = der + 8.d0*d8*c8/rij**9 + 10.d0*d10*c10/rij**11
       dd1 = dd(1,dmp1,rij)
       dd6 = dd(6,dmp6,rij)
       dd8 = dd(8,dmp8,rij)
       dd10 = dd(10,dmp10,rij)
       der = der + dd1*qa*qb/rij - dd6*c6/(rij)**6 -
     1             dd8*c8/rij**8 - dd10*c10/rij**10
       valder(numt) = der
       if(beta.ne.0.d0) then
        valder(numt) = valder(numt) - beta*a*dexp(-beta*rij)
       if(numt.ne.1) then
        do kk=1,npowers
         valder(kk) = (dfloat(kk)/rij - beta)*values(kk)
        end do
       endif
       endif
 334  continue
      return
      endif
      
      

      
c ***********************************************************
c  ntpot = 21 : QPEN potential...
c ********* **************************************************
      if(ntpot.eq.21) then
c Parameters:
      ba = parab(1,itypea(ia),itypea(ia))
      bb = parab(1,itypeb(ib),itypeb(ib))
      beta = 0.5d0*(ba + bb)
      aa = parab(2,itypea(ia),itypea(ia))
      ab = parab(2,itypeb(ib),itypeb(ib))
      a = dsqrt(aa*ab)
      ca = parab(3,itypea(ia),itypea(ia))
      cb = parab(3,itypeb(ib),itypeb(ib))
      c6 = dsqrt(ca*cb)
      qa = param(1,itypea(ia))/18.22262373d0
      qb = param(1,itypeb(ib))/18.22262373d0
c
      numt = 1
      values(1) = 0.d0
      if(beta.gt.0.d0) then
       values(1) = values(1) + qa*qb*a*dexp(-beta*rij)
      endif
      values(1) = values(1)+ 332.0719d0*qa*qb/rij -qa*qb*c6/(rij)**6 
      return
      endif
c ***********************************************************
c  ntpot = 18 : exponential + r^-6 + elst, for the MP2 potential (CH3CN)
c ********* **************************************************
      if(ntpot.eq.18) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
c !!!!! a, c6, and dij give energy in kJ/mol ! transform to kcal/mol
      a = parab(2,itypea(ia),itypeb(ib))
      c6 = parab(3,itypea(ia),itypeb(ib))
      a = 0.239006d0*a
      c6 = 0.239006d0*c6
c      dij = parab(4,itypea(ia),itypeb(ib))
c      dij = 0.239006d0*dij
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
      dij = qa*qb
c
      numt = 1
      values(1) = 0.d0
      if(beta.gt.0.d0) then
       values(1) = values(1) + a*dexp(-beta*rij)
      endif
      values(1) = values(1) + dij/rij - c6/(rij)**6 
      return
      endif


c ***********************************************************
c  ntpot = 1888 : exponential + r^-6 + elst, combining rules
c ********* **************************************************
      if(ntpot.eq.1888) then
c Parameters:
      betaa = parab(1,itypea(ia),itypeb(ia))
      betab = parab(1,itypea(ib),itypeb(ib))
c !!!!! a, c6, and dij give energy in kJ/mol ! transform to kcal/mol
      aa = parab(2,itypea(ia),itypeb(ia))
      ab = parab(2,itypea(ib),itypeb(ib))
      c6a = parab(3,itypea(ia),itypeb(ia))
      c6b = parab(3,itypea(ib),itypeb(ib))      
      
      if (betaa.gt.0.0d0 .and. betab.gt.0.0d0) then
        beta=(betaa+betab)*0.5d0
      else
        beta=0.0d0
      end if
      a=sqrt(aa*ab)
      c6=sqrt(c6a*c6b)
      
c      write(*,*)itypea(ia),itypea(ib),beta,a,c6
      
      a = 0.239006d0*a
      c6 = 0.239006d0*c6
c      dij = parab(4,itypea(ia),itypeb(ib))
c      dij = 0.239006d0*dij
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
      dij = qa*qb
c
      numt = 1
      values(1) = 0.d0
      if(beta.gt.0.d0) then
       values(1) = values(1) + a*dexp(-beta*rij)
      endif
      values(1) = values(1) + dij/rij - c6/(rij)**6 
      
      return
      endif

c ***********************************************************
c  ntpot = 18881 : exponential + r^-6 + elst, combining rules
c ********* **************************************************
      if(ntpot.eq.18881) then
c Parameters:
      betaa = parab(1,itypea(ia),itypeb(ia))
      betab = parab(1,itypea(ib),itypeb(ib))
c !!!!! a, c6, and dij give energy in kJ/mol ! transform to kcal/mol
      aa = parab(2,itypea(ia),itypeb(ia))
      ab = parab(2,itypea(ib),itypeb(ib))
      c6a = parab(3,itypea(ia),itypeb(ia))
      c6b = parab(3,itypea(ib),itypeb(ib))      
      
      if (betaa.gt.0.0d0 .and. betab.gt.0.0d0) then
        beta=(betaa+betab)*0.5d0
      else
        beta=0.0d0
      end if
      a=sqrt(aa*ab)
      c6=sqrt(c6a*c6b)
      
c      write(*,*)itypea(ia),itypea(ib),beta,a,c6
      
      a = 0.239006d0*a
      c6 = 0.239006d0*c6
c      dij = parab(4,itypea(ia),itypeb(ib))
c      dij = 0.239006d0*dij
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
      dij = qa*qb
c
      numt = 1
      values(1) = 0.d0
c      if(beta.gt.0.d0) then
c       values(1) = values(1) + a*dexp(-beta*rij)
c      endif
      values(1) = values(1) + dij/rij - c6/(rij)**6 
      
      return
      endif



c ***********************************************************
c  ntpot = 1818 : exponential + r^-6 + elst, for the MP2 potential (CH3CN)
c ********* **************************************************
      if(ntpot.eq.1818) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
c !!!!! a, c6, and dij give energy in kJ/mol ! transform to kcal/mol
      a = parab(2,itypea(ia),itypeb(ib))
      c6 = parab(3,itypea(ia),itypeb(ib))
      a = 0.239006d0*a
      c6 = 0.239006d0*c6
c      dij = parab(4,itypea(ia),itypeb(ib))
c      dij = 0.239006d0*dij
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
      dij = qa*qb
c
      numt = 1
      values(1) = 0.d0
c      if(beta.gt.0.d0) then
c       values(1) = values(1) + a*dexp(-beta*rij)
c      endif
      values(1) = values(1)  - c6/(rij)**6 
      return
      endif



c ***********************************************************
c  ntpot = 188 : exponential + r^-6 + elst, partially linear
c ********* **************************************************
      if(ntpot.eq.188) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
c !!!!! a, c6, and dij give energy in kJ/mol ! transform to kcal/mol
c      a = parab(2,itypea(ia),itypeb(ib))
      c6 = parab(3,itypea(ia),itypeb(ib))
      a = 0.239006d0*a
      c6 = 0.239006d0*c6
c      dij = parab(4,itypea(ia),itypeb(ib))
c      dij = 0.239006d0*dij
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
      dij = qa*qb
c
      numt = 0
      if(beta.gt.0.d0) then
       numt = numt + 1
       values(numt) = dexp(-beta*rij)*0.239006d0
      endif
      numt = numt + 1
      values(numt) = dij/rij - c6/(rij)**6 
      return
      endif


c ***********************************************************
c  ntpot = 16 : exponential + r^-(6-8-10) + elst, totally nonlinear..
c ********* **************************************************
      if(ntpot.eq.16) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
c      dmp1 = parab(6,1,1)
      dmp6 = parab(7,itypea(ia),itypeb(ib))
c      dmp6 = parab(7,1,1)
      dmp8 = parab(8,itypea(ia),itypeb(ib)) 
c      dmp8 = parab(8,1,1)
      dmp10= parab(9,itypea(ia),itypeb(ib))
c      dmp10 = parab(9,1,1)
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      values(1) = 0.d0
      if(beta.gt.0.d0) then
       values(1) = values(1) + a*dexp(-beta*rij)
      endif
       d1 = d(1,dmp1,rij)
       d6 = d(6,dmp6,rij)
       d8 = d(8,dmp8,rij)
       d10 = d(10,dmp10,rij)
      values(1) = values(1) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
     1 d8*c8/rij**8 - d10*c10/rij**10
      return
      endif
c ***********************************************************
c  ntpot = 14 : exponential + r^-6 + elst, totally nonlinear..
c ***********************************************************
      if(ntpot.eq.14) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(5,itypea(ia),itypeb(ib))
      dmp = parab(4,itypea(ia),itypeb(ib))
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 1
      values(1) = 0.d0
      pom = 0.d0
      if(beta.gt.0.d0) then
       pom = a*dexp(-beta*rij)
       values(1) = values(1) + pom
      endif
c        x = 0.2d0*rij   
c        CALL ERROR(X, ERRORF)
       d1 = d(1,dmp,rij)
c       d3 = d(3,dmp,rij)
c       d4 = d(4,dmp,rij)       
       d6 = d(6,dmp,rij)
       d8 = d(8,dmp,rij)
c       d10 = d(10,dmp,rij)
c      errorf = 1.d0
c      values(1) = values(1) + errorf*(qa*qb/rij - c6/(rij)**6)
      values(1) = values(1) + d1*qa*qb/rij - d6*c6/(rij)**6 - 
     1 d8*c8/rij**8 - d10*c10/rij**10
c compute the derivative over r_ab
c General expression to be modified
c        der = -b*pom - dmp1*c1*r1*r1 - 6*dmp6*c6*r6*r1
c        der = der - 7*dmp7*c7*r7*r1 - 8*dmp8*c8*r8*r1
c        der = der - 9*dmp9*c9*r9*r1 - 10*dmp10*c10*r10*r1
c        der = der + c1*r1*dd(1,d1,r) + c6*r6*dd(6,d6,r)
c        der = der + c7*r7*dd(7,d7,r) + c8*r8*dd(8,d8,r)
c        der = der + c9*r9*dd(9,d9,r) + c10*r10*dd(10,d10,r)
c ntpot = 14 (DMNA_CH3CN) - specific stuff
        der = -beta*pom - d1*qa*qb*r1*r1 + 6*d6*c6/rij**7
        der = der + 8*d8*c8/rij**9
        der = der + qa*qb*dd(1,dmp,r)/rij - c6*dd(6,dmp,r)/rij**6
        der = der - c8*dd(8,dmp,r)/rij**8
      go to 112
      write(6,*)'atoms',ia,ib,itypea(ia),itypeb(ib)
      write(6,*)'rij',rij
      write(6,*)'dmp',dmp,dmp,dmp
      write(6,*)'d',d1,d6,d8
      write(6,*)'C  ',c6,c8
      write(6,*)'a,beta',a,beta
      write(6,*)'values(1)',values(1)
c      values(1)=values(1) + d3*c3/rij**3 + d4*c4/rij**4 - d6*c6/rij**6
 112  continue
      return
      endif
c ***********************************************************
c  ntpot = 13 : exponential + r^-6 + elst, partially linear..
c ***********************************************************
      if(ntpot.eq.13) then
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
      dswitch = parab(2,itypea(ia),itypeb(ib))
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
      numt = 0
      values(1) = 0.d0
      if(beta.gt.0.d0) then
       numt = numt+1
       values(numt) = dexp(-beta*rij)
      endif
      if(dswitch.ne.0.d0) then
       numt = numt + 1
       values(numt) =  -1.d0/(rij)**6
      endif
      x = 0.2d0*rij
      CALL ERROR(X, ERRORF)
      numt = numt + 1
      values(numt) = errorf*qa*qb/rij 
      return
      endif
c *************************************************************
c ntpot = -11: the charge density + asymptotics model. Totally
c nonlinear. Asymptotics is not of the site-site form, hence--
c it must be added in the calling procedure.
c ************************************************************* 
      if(ntpot.eq.-11) then
      pi = dacos(-1.d0)
      qa = param(1,itypea(ia))
      exa = param(2,itypea(ia))
      ca = param(3,itypea(ia))
      qb = param(1,itypeb(ib))
      exb = param(2,itypeb(ib))
      cb = param(3,itypeb(ib))
      numt = 1
      val = 0.d0
      if(ca*cb.ne.0.d0) then
c --- el-el part
       pom = pi/(exa+exb)
       pom = dsqrt(pom*pom*pom)
       pom = pom*dexp(-exa*exb*rij*rij/(exa+exb))
       val = val + ca*cb*pom
      endif
      if(qb*ca.ne.0.d0) then
c --- el_A -- q_B part
       val = val - qb*ca*dexp(-exa*rij*rij)
      endif
      if(qa*cb.ne.0.d0) then
c --- el_B -- q_A part
       val = val - qa*cb*dexp(-exb*rij*rij)
      endif
      itmp=4
      values(1) = param(itmp,1)*val
      return
      endif
      if(ntpot.eq.-12) then
c *************************************************************
c ntpot = -12: the exp-site-site + charge + disp-ind asymptotics model. 
c Totally nonlinear. Asymptotics is not of the site-site form, hence--
c it must be added in the calling procedure.
c *************************************************************
c Parameters:
      beta = parab(1,itypea(ia),itypeb(ib))
c      a = parab(2,itypea(ia),itypeb(ib))
      a = dexp(parab(2,itypea(ia),itypeb(ib)))
      a1 = parab(3,itypea(ia),itypeb(ib))
c      qa = param(1,itypea(ia))
c      qb = param(1,itypeb(ib))
       qa = 0.d0
       qb = 0.d0
c
      numt = 1
      values(1) = 0.d0
      if(beta.gt.0.d0) then
       values(1) = values(1) + (1.d0+a1*rij)*a*dexp(-beta*rij)
      endif
c        x = 0.2d0*rij
c        CALL ERROR(X, ERRORF)
c      d1 = d(1,beta,rij)
c      d6 = d(6,beta,rij)
      errorf = 0.d0
      values(1) = values(1) + errorf*qa*qb/rij 
c
      return
      endif
c *************************
c ntpot = -2 ---> the asymptotic dispersion (or elst) linear fitting
c *************************
      if(ntpot.eq.-2) then
       numt = 3
       values(1) = -1.d0/rij**6
       values(2) = -1.d0/rij**8
       values(3) = -1.d0/rij**10
c       values(1) = 1.d0/rij**3
c       values(2) = 1.d0/rij**6
      return
      endif
c *****************************************
c ntpot = -14 ----> the asymptotic electrostatics fitting...
c *****************************************
      if(ntpot.eq.-14) then
       numt = 1
       qa = param(1,itypea(ia))
       qb = param(1,itypeb(ib))
       values(1) = qa*qb/rij
      return
      endif
c *****************************************
c ntpot = 200 ----> the MCY potential
c *****************************************
      if(ntpot.eq.200) then
      b1 = parab(1,itypea(ia),itypeb(ib))/a0
      b2 = parab(2,itypea(ia),itypeb(ib))/a0
      a1 = parab(3,itypea(ia),itypeb(ib))*au2cal
      a2 = parab(4,itypea(ia),itypeb(ib))*au2cal
       qa = param(1,itypea(ia))
       qb = param(1,itypeb(ib))
       numt = 1
       pom1 = dexp(-b1*rij)
       pom2 = dexp(-b2*rij)
       values(1) = qa*qb/rij + a1*pom1 - a2*pom2
       valder(1) = -qa*qb/(rij*rij) - b1*a1*pom1 + b2*a2*pom2
      return
      endif
      end
c
C **********************************************************************
        SUBROUTINE ERROR(X, ERFC)
C **********************************************************************

C    *******************************************************************
C    ** APPROXIMATION TO THE COMPLEMENTARY ERROR FUNCTION             **
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,    **
C    **    NATIONAL BUREAU OF STANDARDS, FORMULA 7.1.26               **
C    *******************************************************************

        IMPLICIT DOUBLE PRECISION (A-H, O-Z)
        IMPLICIT INTEGER (I-N)
        DOUBLE PRECISION        A1, A2, A3, A4, A5, P

        PARAMETER ( A1 = 0.254829592, A2 = -0.284496736 )
        PARAMETER ( A3 = 1.421413741, A4 = -1.453152027 )
        PARAMETER ( A5 = 1.061405429, P  =  0.3275911   )

        DOUBLE PRECISION        T, X, XSQ, TP, TXP

        COMMON /TRANSFER/T, TXP, TP 

C    *******************************************************************

        T  = 1.0 / ( 1.0 + P * X )
        XSQ = X * X

        TP = T * ( A1 + T * ( A2 + T * ( A3 + T * ( A4 + T * A5 ) ) ) )

        TXP = EXP ( -XSQ )
        ERFC = TP * TXP

        RETURN
        END
c
c---- Compute the damping constant for a given geometry
c
        subroutine dampcomp(oa,ob,dmp)
        implicit real*8 (a-h,o-z)
        parameter (maxb=500,maxp=140000)
        parameter (mxl=100)
        dimension oa(3),ob(3),oc(2)
        common/damp/ ldum(5,mxl),cdum(maxb),ndum
        data oc /2*0.d0/
c
        iperm = 1   
        dmp = 0.d0
        do i=1,ndum
         la = ldum(1,i)
         ka = ldum(2,i)
         lb = ldum(3,i)
         kb = ldum(4,i)
         l  = ldum(5,i)
         glam = almre(la,ka,lb,kb,l,oa,ob,oc)
c -- symmetrize if molecules identical...
         if(iperm.eq.1) then
          iphase = (-1)**(la+lb)
          glam = glam + iphase*almre(lb,kb,la,ka,l,oa,ob,oc)
         endif
         dmp = dmp - cdum(i)*glam
        end do
c        if(dmp.lt.0.d0) write(6,*)'Wrong damping!',dmp
        dmp = dabs(dmp)
        return
        end
c
c   Calculate the total damped asymptotics from elst and dispind..
c   th1, th2, phi in radians...., R in Angstroms.
c
        subroutine asymp(rpt,value,valder)
        implicit real*8 (a-h,o-z)
       parameter (mxl=100)
       parameter (maxb=500,maxp=140000)
        dimension rpt(6),oa(3), ob(3), oc(2)
        dimension en(20), el(20), ei(20), ed(20)
      common/damp/ ldum(5,mxl),cdum(maxb),ndum
      common/fit/rrpt(6,maxp),dat0(maxp),sigma(maxp),iweight,ntpot,
     1 idonl,iopt
      data oa(1) /0.d0/, oc /0.d0,0.d0/

c
        iperm = 0
        r = rpt(1)
        oa(2) = rpt(2)
        oa(3) = rpt(3)
        ob(1) = rpt(4)
        ob(2) = rpt(5)
        ob(3) = rpt(6)
c
c---- call the asymptotics procedures 
c
        call dispind(rpt,ed,ei)
        call elst(rpt,el)
c        call shrtas(R,th1,th2,phi,srval)
c
c---- Compute the damping constant for a given geometry
c
        dump = 0.d0
        do i=1,ndum
         la = ldum(1,i)
         ka = ldum(2,i)
         lb = ldum(3,i)
         kb = ldum(4,i)
         l  = ldum(5,i)
c glamc changed into glam to test the real alm version...
         glam = almre(la,ka,lb,kb,l,oa,ob,oc)
c -- symmetrize iif molecules identical...
         if(iperm.eq.1) then
          iphase = (-1)**(la+lb)
          glam = glam + iphase*almre(lb,kb,la,ka,l,oa,ob,oc)
         endif
c         glam = real(glamc)
         dump = dump - cdum(i)*glam
        end do
c
c--- Damping factor ready, proceed to the asymptotics...
c
        do i=1,20
         en(i) = el(i) + ed(i) + ei(i)
        end do
        value = 0.d0
        valder = 0.d0
        do i=1,20
         ddd = d(i,dump,r)
         value = value + ddd*en(i)
       valder = valder + en(i)*(dd(i,dump,r)-i*ddd/r)
        end do
       return
       end
c
      function d(n,beta,r)
c
c     calculate the damping factor (small R correct)
c
      implicit real*8 (a-h,o-z)
      br=beta*r
c The following line added by RB, Sept. 18 1997
      if(br.eq.0.d0) then
       d = 0.d0
       return
      endif
      sum=1.0d0
      term=1.0d0
      ncn=n
      do i=1,ncn
        term=term*br/i
        sum=sum+term
      enddo
      d=1.0d0 - dexp(-br)*sum
c     in case of d --> 0 use
c     d=1.0d0 - dexp(-br)*sum = sum_m=ncn+1^\infty br^m/m!
      if(dabs(d).lt.1.0d-8) then
        d=0.0d0
        do i=ncn+1,1000
          term=term*br/i
          d=d+term
          if(term/d .lt. 1.0d-8) go to 111
        enddo
        write(6,*) 'No convergence in d'
  111 continue
      d=d*dexp(-br)
      endif
c     write(6,'(i4,2f10.5,e20.10)') n,beta,r,d
      return
      end
c
      function dd(n,b,r)
      implicit real*8 (a-h,o-z)
      common/factorial/ f(0:40)
      br = b*r
      dd = b*dexp(-br)*(br)**n/f(n)
      return
      end

c
c Subroutine for the calculation of the multipole
c part of the electrostatic energy for a given
c dimer conformation. The angles are expressed
c in radians. 
c
        subroutine elst(rpt,el)
        implicit real*8 (a-h,o-z)
       parameter (maxb=500,maxp=140000)
        dimension rrev(20), oa(3), ob(3), oc(2), rpt(6), el(20)
        common/factorial/ fac(0:40)
        common/moments/ qa(100,3), qb(100,3),la(maxb),ka(maxb),lb(maxb),
     1 kb(maxb),nmoma, nmomb
        data Pi /3.1415926535897932D0/
        data efact /627.51d0/
        data a0 /0.529177249d0/
        data oa(1) /0.d0/, oc /0.d0,0.d0/
        data sqt2 / 1.414213562373d0/
c
        r = rpt(1)
        oa(2) = rpt(2)
        oa(3) = rpt(3)
        ob(1) = rpt(4)
        ob(2) = rpt(5)
        ob(3) = rpt(6)
        rrev(1) = 1.d0/(R/a0)
        do i=2,20
         rrev(i) = rrev(i-1)*rrev(1)
        end do
c
         els = 0.d0
        do i=1,20
         el(i) = 0.d0
        end do
c
        do ia=1,nmoma
        do ib=1,nmomb
         phase = fac(2*la(ia)+2*lb(ib)+1)
         phase = phase/(fac(2*la(ia))*fac(2*lb(ib)))
         phase = (-1.d0)**la(ia) * dsqrt(phase)
         ll = la(ia) + lb(ib)
         rrr = rrev(ll+1)
         ipa = (-1)**ka(ia)
         ipb = (-1)**kb(ib)
         ikp = ka(ia)*kb(ib)
         iks = ka(ia) + kb(ib)
         aa  = ipa*almre(la(ia),ka(ia),lb(ib),kb(ib),ll,OA,OB,OC)
         if(ikp.ne.0)
     1 aa = aa+almre(la(ia),-ka(ia),lb(ib),kb(ib),ll,OA,OB,OC)
c calculate the product of multipole moments at the proper level
c (currently -- MP3-resp minus the SCF part) 
c         qpr = qa(ia,1)*qb(ib,1) + qa(ia,1)*(qb(ib,2)+qb(ib,3))
c         qpr =  qa(ia,1)*(qb(ib,2)+qb(ib,3))
c         qpr = qpr + (qa(ia,2)+qa(ia,3))*qb(ib,1)
         qpr=(qa(ia,1)+qa(ia,2)+qa(ia,3))*(qb(ib,1)+qb(ib,2)+qb(ib,3))
c         qpr = qa(ia,1)*qb(ib,1)
         if((ikp.eq.0).and.(iks.ne.0)) qpr = sqt2*qpr
c limitation to R^-15 in electrostatics.... (15 never reached anyway)
         if((ll+1).le.15)
c     1   el(ll+1) = el(ll+1) + phase*rrr*real(a)*qpr
     1   el(ll+1) = el(ll+1) + phase*rrr*aa*qpr*ipb
        end do
        end do
c
        do i=1,20
         el(i) = efact*el(i)
        end do
c
c That's all
c
        return
        end
c
c Read the multipole moment components
c
        subroutine rdmom
        implicit real*8 (a-h,o-z)
       parameter (maxb=500,maxp=140000)
        common/moments/ qa(100,3), qb(100,3),la(maxb),ka(maxb),lb(maxb),
     1 kb(maxb),nmoma, nmomb
c
        write(6,*)'RDmom entered'
c---- read the moments for molecule A
        open(unit=3,file='moments.a',form='formatted')
c
        do i=1,1000
         read(3,*,END=100) la(i),ka(i),(qa(i,j),j=1,3)
        end do
 100    continue
        nmoma = i - 1
        write(6,*)nmoma,' A moments read in:'
        write(6,*)' L            SCF             MP2             MP3'
        do i=1,nmoma
         write(6,'(2i3,3e16.7)') la(i),ka(i),(qa(i,j),j=1,3)
        end do
        close(3)
c---- read the moments for molecule B
        open(unit=3,file='moments.b',form='formatted')
c
        do i=1,1000
         read(3,*,END=101) lb(i),kb(i),(qb(i,j),j=1,3)
        end do
 101    continue
        nmomb = i - 1
        write(6,*)nmomb,' B moments read in:'
        write(6,*)' L            SCF             MP2             MP3'
        do i=1,nmomb
         write(6,'(2i3,3e16.7)') lb(i),kb(i),(qb(i,j),j=1,3)
        end do
        close(3)
c
        return
        end

c
c Subroutine for the calculation of the multipole
c part of the dispersion and induction energies for a given
c dimer conformation. 
c
        subroutine dispind(rpt,eld,eli)
        implicit real*8 (a-h,o-z)
        parameter (maxc=10000)
        dimension rr(20),oa(3),ob(3),oc(2), rpt(6), eld(20), eli(20)
c        complex*16 a, alm
c--- common/dind/ contains dispersion and induction coefs.
        common/dind/ cd(maxc),ci(maxc),ld(6,maxc),li(6,maxc),
     1 idisp,iind
        data Pi /3.1415926535897932D0/
        data a0 /0.529177249d0/
        data efact /627.51d0/        
        data oa(1) /0.d0/, oc /0.d0,0.d0/
c
        r = rpt(1)
        oa(2) = rpt(2)
        oa(3) = rpt(3)
        ob(1) = rpt(4)
        ob(2) = rpt(5)
        ob(3) = rpt(6)
        rrev = 1.d0/(R/a0)
        rr(1) = rrev
        do i=2,12
         rr(i) = rr(i-1)*rrev
        end do
        do i=1,20
         eld(i) = 0.d0
         eli(i) = 0.d0
        end do
c
        do i=1,idisp
         ka = ld(2,i)
         kb = ld(4,i)
         if(ka.lt.0) go to 14
         if((ka.eq.0).and.kb.lt.0) go to 14
         la = ld(1,i)
         lb = ld(3,i)
         l  = ld(5,i)
         n  = ld(6,i)
         aa = almre(la,ka,lb,kb,l,oa,ob,oc)
         if((ka.eq.0).and.(kb.eq.0)) aa = 0.5d0*aa
c         eld(n) = eld(n) + cd(i)*rr(n)*real(a)
         eld(n) = eld(n) + 2*cd(i)*rr(n)*aa
 14     continue
        end do
c
        do i=1,iind
         ka = li(2,i)
         kb = li(4,i)
         if(ka.lt.0) go to 12
         if((ka.eq.0).and.kb.lt.0) go to 12
         la = li(1,i)
         lb = li(3,i)
         l  = li(5,i)
         n  = li(6,i)
         aa = almre(la,ka,lb,kb,l,oa,ob,oc)
         if((ka.eq.0).and.(kb.eq.0)) aa = 0.5d0*aa
c         eli(n) = eli(n) + ci(i)*rr(n)*real(a)
         eli(n) = eli(n) + 2*ci(i)*rr(n)*aa
 12     continue
        end do
        do i=1,20
         eli(i) = -efact*eli(i)
         eld(i) = -efact*eld(i)
        end do
c
c That's all
c
        return
        end
c
c Read the dispersion ond induction coefficients.
c It is assumed that the induction coefficients are
c complete, i.e., that both A->B and B->A are included,
c with the proper phases after the induct calculation...
c
        subroutine rdcoef
        implicit real*8 (a-h,o-z)
        parameter (maxc=10000)
        common/dind/ cd(maxc),ci(maxc),ld(6,maxc),li(6,maxc),
     1 idisp,iind
c
        write(6,*)'RDcoef entered'
        open(unit=1,file='coefd.dat',form='formatted')
        open(unit=2,file='coefi.dat',form='formatted')
c
        do i=1,100000
         read(1,*,end=100)(ld(j,i),j=1,6), c0,c1,ctot
c !!!! update !!!!
         cd(i) = ctot
        end do
 100    continue
        idisp = i-1
        write(6,*)idisp,' nonzero dispersion coefficients read in'
c
        do i=1,100000
         read(2,*,END=200)(li(j,i),j=1,6),c0,c1,ctot
c !!!! update !!!!
         ci(i) = ctot
        end do
 200    continue
        iind = i-1
        write(6,*)iind,' nonzero induction coefficients read in'
c
        close(1)
        close(2)
        return
        end
c
      subroutine rddamp
      implicit real*8 (a-h,o-z)
       parameter (mxl=100)
       parameter (maxb=500,maxp=140000)
      common/damp/ ldum(5,mxl),cdum(maxb),ndum
c
      open(unit=4,file='dampbas',form='formatted')
c
        do i=1,1000
         read(4,*,END=100)(ldum(j,i),j=1,5), cdum(i)
        end do
 100    continue
        ndum = i-1
        write(6,*)'The damping factor basis and coefs:'
        do i=1,ndum
         write(6,10)(ldum(j,i),j=1,5), cdum(i)
        end do
        close(4)
 10     format(5i4,e20.10)
       return
       end
c
c compute the matrix of N!
c
        subroutine fct(nmax)
        implicit real*8 (a-h,o-z)
        common/factorial/ f(0:40)
c
        f(0) = 1.d0
        do i=1,nmax
         f(i) = f(i-1)*i
        end do
        return
        end
c
c The data procedure modified to read the Gaussian basis parameters.
c
      subroutine data1
      implicit real*8 (a-h,o-z)
c
      parameter (maxb=500,maxp=140000)
      parameter (nsitemax=42,ntypemax=25,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=3,maxpar2=13)
      parameter (mxl=100)
      parameter (maxn=6*nsitemax+maxpar1*ntypemax+maxpar2*ntypemax2)
c
      character*8 label(9)
      common/npar/param(maxpar1,ntypemax),
     1 parab(maxpar2,ntypemax,ntypemax),nparm,nparab
      common/sites/ sitea(3,nsitemax), siteb(3,nsitemax),nsitea,nsiteb
      common/types/ itypea(nsitemax),itypeb(nsitemax)
      common/optdrv/ iosita(2,3*nsitemax),iositb(2,3*nsitemax),
     1 iopar(2,ntypemax),ioparab(3,ntypemax2),nosa,nosb,nopar,noparab
c
      common/fit/rpt(6,maxp),dat0(maxp),sigma(maxp),iweight,ntpot,
     1 idonl, iopt
c
      common/leastsq/dat(maxp),a(maxp,maxb),a0(maxp),g(maxp),c(maxb),
     1 ata(maxb*maxb),b(maxb),works(2*maxb),iworks(maxb),np,nb
      COMMON/EOPT/TOLF,TOLR,ANOISE,TIMLIM,LINP,LOUT,ALFST,IPR,NPAR,
     &            SAFETL
      common/inou/ icon(80)
      common/grdfit/ rlow(maxn),rhigh(maxn),npts(maxn)
      common/readyasm/ asmval(15,maxp),iasdone
      common/damp/ ldum(5,mxl),cdum(maxb),ndum
      common/misc/ R_0,isyst,npowers
      data bohr2a /0.529177249d0/
c
      ntreg = 5
c
c Zero out the A and B site position optimization indicators...
c
      do i=1,3*nsitemax
      do j=1,2
       iosita(j,i) = 0
       iositb(j,i) = 0
      end do
      end do 
c
c Zero out the nonlinear parameters matrices
c
      do i=1,ntypemax
      do k=1,maxpar1
       param(k,i) = 0.d0
      end do
      end do
      do i=1,ntypemax
      do j=1,ntypemax
      do k=1,maxpar2
       parab(k,i,j) = 0.d0
      end do
      end do
      end do
c
c Zero out the nonl params optimization indicators...
c
      do i=1,ntypemax
      do k=1,2
       iopar(k,i) = 0
      end do
      end do
      do i=1,ntypemax2
      do k=1,3
       ioparab(k,i) = 0
      end do
      end do
c
      nopar = 0
	  open(5, file = "hemetfor", status = "old")
      read(5,*)iunit
c
c If iunit = 0, then lengths in bohrs, pseudofunction coefs normalized to 22
c (No of electrons for CO2), iunit = 1, then lengths in A, pseudofunctions
c normalized to unity... 
c
c read in the sites of monomer A
c
      nosa = 0
      read(5,*) nsitea,ngaussa
      write(6,*)'nsitea read in',nsitea
      do i=1,nsitea
       read(5,*) itypea(i),(sitea(j,i),j=1,3),k1,k2,k3
       if(iunit.eq.0) then
c transform to Angstroms...
       do j=1,3
        sitea(j,i) = bohr2a*sitea(j,i)
       end do
       endif  
       if(k1.ne.0) then
        nosa = nosa + 1
        iosita(1,nosa) = i
        iosita(2,nosa) = 1
       endif
       if(k2.ne.0) then
        nosa = nosa + 1
        iosita(1,nosa) = i
        iosita(2,nosa) = 2
       endif
       if(k3.ne.0) then
        nosa = nosa + 1
        iosita(1,nosa) = i
        iosita(2,nosa) = 3
       endif
      end do
c
c Read the Gaussian basis of monomer A
c assume the number of electrons for monomer A
      dnela = 22.d0
      do i=nsitea+1,nsitea+ngaussa
       read(5,*) itypea(i),(sitea(j,i),j=1,3),expon,k1,k2,k3,k4,coef
c       itypea(i) = ntreg+i-nsitea
c transform to Angstroms...
       if(iunit.eq.0) then
       do j=1,3
        sitea(j,i) = bohr2a*sitea(j,i)
       end do
       param(2,itypea(i)) = expon/bohr2a**2
       param(3,itypea(i)) = coef/(dsqrt(dnela)*dsqrt(bohr2a)**3)
       else
       param(2,itypea(i)) = expon
       param(3,itypea(i)) = coef
       endif
c
       if(k1.ne.0) then
        nosa = nosa + 1
        iosita(1,nosa) = i
        iosita(2,nosa) = 1
       endif
       if(k2.ne.0) then
        nosa = nosa + 1
        iosita(1,nosa) = i
        iosita(2,nosa) = 2
       endif
       if(k3.ne.0) then
        nosa = nosa + 1
        iosita(1,nosa) = i
        iosita(2,nosa) = 3
       endif
       if(k4.ne.0) then
        nopar = nopar + 1
        iopar(1,nopar) = itypea(i)
        iopar(2,nopar) = 2
        nopar = nopar + 1
        iopar(1,nopar) = itypea(i)
        iopar(2,nopar) = 3
       endif
      end do
      nsitea = nsitea + ngaussa
      write(6,*)'nsitea after',nsitea
c
c read in the sites of monomer B
c
      nosb = 0
      read(5,*) nsiteb,ngaussb
      do i=1,nsiteb
       read(5,*) itypeb(i),(siteb(j,i),j=1,3),k1,k2,k3
c transform to Angstroms...
       if(iunit.eq.0) then
       do j=1,3
        siteb(j,i) = bohr2a*siteb(j,i)
       end do
       endif
       if(k1.ne.0) then
        nosb = nosb + 1
        iositb(1,nosb) = i
        iositb(2,nosb) = 1
       endif
       if(k2.ne.0) then
        nosb = nosb + 1
        iositb(1,nosb) = i
        iositb(2,nosb) = 2
       endif
       if(k3.ne.0) then
        nosb = nosb + 1
        iositb(1,nosb) = i
        iositb(2,nosb) = 3
       endif
      end do
c
c Read the Gaussian basis of monomer B (for the moment assume
c the same center types)
c Assume number of electrons for monomer b
      dnelb = 22.d0
      do i=nsiteb+1,nsiteb+ngaussb
       read(5,*)itypeb(i),(siteb(j,i),j=1,3),expon,k1,k2,k3,k4,coef
c       itypeb(i) = ntreg+i-nsiteb
c transform to Angstroms...
       if(iunit.eq.0) then
       do j=1,3
        siteb(j,i) = bohr2a*siteb(j,i)
       end do
       param(2,itypeb(i)) = expon/bohr2a**2
       param(3,itypeb(i)) = coef/(dsqrt(dnelb)*dsqrt(bohr2a)**3)
       else
       param(2,itypeb(i)) = expon
       param(3,itypeb(i)) = coef
       endif
c
       if(k1.ne.0) then
        nosb = nosb + 1
        iositb(1,nosb) = i
        iositb(2,nosb) = 1
       endif
       if(k2.ne.0) then
        nosb = nosb + 1
        iositb(1,nosb) = i
        iositb(2,nosb) = 2
       endif
       if(k3.ne.0) then
        nosb = nosb + 1
        iositb(1,nosb) = i
        iositb(2,nosb) = 3
       endif
       if(k4.ne.0) then
        nopar = nopar + 1
        iopar(1,nopar) = itypeb(i)
        iopar(2,nopar) = 2
        nopar = nopar + 1
        iopar(1,nopar) = itypeb(i)
        iopar(2,nopar) = 3
       endif
      end do
      nsiteb = nsiteb + ngaussb
c
c Read in the one-site nonlinear parameters and opt indicators..
c 
c      nopar = 0
      read(5,*) nparm
      do i=1,nparm
       read(5,*) ityp, inumpar, val, indopt
       param(inumpar,ityp) = val
c Convert charges to the proper units 
c Switch off for the total charge density variant.... (everything
c in atomic units)
       if(inumpar.eq.1) param(1,ityp) = 18.22262373d0*param(1,ityp)
       if(indopt.ne.0) then
        nopar = nopar + 1
        iopar(1,nopar) = ityp
        iopar(2,nopar) = inumpar
       endif
      end do
      write(6,*)'One-site read-in'
c
c Read in the two-site nonlinear parameters and opt indicators..
c
      noparab = 0
      read(5,*) nparab
      do i=1,nparab
       read(5,*) ityp1,ityp2,inumpar, val, indopt
       write(6,*) ityp1,ityp2,inumpar, val, indopt
       parab(inumpar,ityp1,ityp2) = val  
       parab(inumpar,ityp2,ityp1) = val
       if(indopt.ne.0) then
        noparab = noparab + 1
        ioparab(1,noparab) = ityp1
        ioparab(2,noparab) = ityp2  
        ioparab(3,noparab) = inumpar
       endif
      end do
c
c Read the type of the potential form to be used, inear coeff. reading
c indicator, and the optimization parameters....
c
      read(5,*) ntpot, idonl, iopt, iweight, iasdone
      read(5,*) TOLF,TOLR,ANOISE,TIMLIM,LINP,LOUT,ALFST,IPR,
     1            SAFETL,icon(2)
c Read in the misc common
      read(5,*) R_0,isyst,npowers
      if(iopt.eq.2) then
c
c read in the grid upper boundaries and mesh...
c
       read(5,*) ndum
       do i=1,ndum
        read(5,*)rhigh(i),npts(i)
       end do
      endif
c      return
c
      entry output1
c
c Print the read-in parameters....
c
      write(6,*)'            Input parameters'
c
      write(6,*)'Monomer A has',nsitea,' sites:'
      write(6,*)'type         x         y          z'
      do i=1,nsitea
       write(6,10) itypea(i), (sitea(j,i),j=1,3)
      end do
 10   format(2x,i4,3e20.7)
      write(6,*)nosa,' opitmized coordinates for monomer A:'
      do i=1,nosa
       write(6,*)'for site',iosita(1,i),' coordinate',iosita(2,i)
      end do
c
      write(6,*)'Monomer B has',nsiteb,' sites:'
      write(6,*)'type         x         y          z'
      do i=1,nsiteb
       write(6,10) itypeb(i), (siteb(j,i),j=1,3)
      end do
      write(6,*)nosb,' opitmized coordinates for monomer B:'
      do i=1,nosb
       write(6,*)'for site',iositb(1,i),' coordinate',iositb(2,i)
      end do
c
      write(6,*)'One-site nonlinear parameters:'
      do k=1,maxpar1
       write(6,*)'Parameter',k
       do i=1,ntypemax
        pom = param(k,i)
c Convert charges back to atomic units...
        if(k.eq.1) pom = pom/18.22262373d0
        indopt = 0
        do l=1,nopar
         if((iopar(1,l).eq.i).and.(iopar(2,l).eq.k)) indopt = 1
        end do
        if(param(k,i).ne.0.d0.or.indopt.eq.1) then
        if(indopt.eq.0) then
         write(6,*)'Type',i,pom,' fixed'
        else
         write(6,*)'Type',i,pom,' optimized'
        endif
        endif
       end do
      end do
c
      write(6,*)'Two-site nonlinear parameters:'
      do k=1,maxpar2
       write(6,*)'Parameter',k
       do i=1,ntypemax
       do j=1,i
        indopt = 0
        do l=1,noparab
         if((ioparab(1,l).eq.i).and.(ioparab(2,l).eq.j).
     1 and.(ioparab(3,l).eq.k)) indopt = 1
         if((ioparab(2,l).eq.i).and.(ioparab(1,l).eq.j).
     1 and.(ioparab(3,l).eq.k)) indopt = 1 
        end do
        if(parab(k,i,j).ne.0.d0.or.indopt.eq.1) then
        if(indopt.eq.0) then
         write(6,*)'Types',i,j,parab(k,i,j),' fixed'
        else
         write(6,*)'Types',i,j,parab(k,i,j),' optimized'
        endif
        endif
       end do
       end do
       end do
c
       write(6,*)'Potential type used:',ntpot
       write(6,*)'damping parameters'
       do i=1,ndum
        write(6,13)(ldum(k,i),k=1,5),cdum(i)
       end do
       write(6,*)'R_0, isyst, npowers',R_0, isyst, npowers
c
      write(6,*)'Two-site nonlinear parameters - transferable format'
      do k=1,maxpar2
       do i=1,ntypemax
       do j=1,i
        indopt = 0
        do l=1,noparab
         if((ioparab(1,l).eq.i).and.(ioparab(2,l).eq.j).
     1 and.(ioparab(3,l).eq.k)) indopt = 1
         if((ioparab(2,l).eq.i).and.(ioparab(1,l).eq.j).
     1 and.(ioparab(3,l).eq.k)) indopt = 1 
        end do
        if(parab(k,i,j).ne.0.d0.or.indopt.eq.1) then
         write(6,*)i,j,k,parab(k,i,j),indopt
        endif
       end do
       end do
       end do
c
 13     format(5i4,e20.10)
       return
       end
c
c Add the ia,ib pair contribution to the torques on A due to B
c (torquea) and on B due to A (torqueb), and to the force
c on B due to A (forceab)
c
       subroutine caltorque(torquea,torqueb,forceab,der,ia,ib,rij,R)
       implicit real*8 (a-h,o-z)
      parameter (maxb=500,maxp=140000)
      parameter (nsitemax=42,ntypemax=25,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=3,maxpar2=13)
      parameter (mxl=100)
       dimension torquea(3),torqueb(3),rab(3),pomvec(3)
       dimension forceab(3)
       common/temp_sites/ siteat(3,nsitemax), sitebt(3,nsitemax)
c
c Calculate the rab vector as R_ab*der/rij
c
       do i=1,3
        rab(i) = der*(sitebt(i,ib) - siteat(i,ia))/rij
        forceab(i) = forceab(i) - rab(i)
       end do
c The increment to the A torque
       call vecp(siteat(1,ia),rab,pomvec)
       do i=1,3
        torquea(i) = torquea(i) + pomvec(i)
       end do
c The increment to the B torque -- shift the z coordinate back...
       sitebt(3,ib) = sitebt(3,ib) - R
       call vecp(rab,sitebt(1,ib),pomvec)
       do i=1,3
        torqueb(i) = torqueb(i) + pomvec(i)
       end do       
       sitebt(3,ib) = sitebt(3,ib) + R
c that's it, exit
       return
       end
c
c Calculate the vector product of vectors a and b
c
        subroutine vecp(a,b,c)
        implicit real*8 (a-h,o-z)
        dimension a(3), b(3), c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = -a(1)*b(3) + a(3)*b(1)
        c(3) = a(1)*b(2) - a(2)*b(1)
        return
        end
c
c The two-body fluctuating-charge induction model for H2O.
c Provided by Michael Smith on May 4, 2001. Modified by RB
c to connect to the fitting program ssfit.

C --------------------------------------------------------------------------

      subroutine POLDER(R,coor1,coor2,E,forceab,torquea,torqueb)
      implicit real*8 (a-h,o-z)
c
      parameter (maxb=500,maxp=140000)
      parameter (nsitemax=42,ntypemax=25,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=3,maxpar2=13)
      parameter (mxl=100)
c The parameters needed for damping....
      common/npar/param(maxpar1,ntypemax),
     1 parab(maxpar2,ntypemax,ntypemax),nparm,nparab
      common/types/ itypea(nsitemax),itypeb(nsitemax)
c
      DOUBLE PRECISION SUM,E,poma,pomb,R
      DOUBLE PRECISION rot(3,3),coor1(3,8), coor2(3,8),geom(3,5),
     +       rpt(8,8),w(5),w2(5),wp(5,3),w2p(5,3),forceab(3),
     +       torquea(3),torqueb(3)
      DOUBLE PRECISION A(5,5),Q(5),COORDS(7),AGEOM(5,3)
      DOUBLE PRECISION r_a(3),r_b(3),r_ab(3),vpoma(3),vpomb(3),
     +                 wta(5,3),w2ta(5,3),
     +                 wtb(5,3),w2tb(5,3)
      data ((A(I,J),I=1,5),J=1,5),Q,COORDS,((AGEOM(K,L),K=1,5),L=1,3)/
     *  75.92287318493169d0,   5.59411445413799d0,  5.59411445413799d0,      ! response matrix A
     * -43.55555104660383d0, -43.55555104660384d0,
     *   5.59411445413800d0,   3.22751022617862d0,  1.20983624816850d0,
     *  -5.01573046424256d0,  -5.01573046424256d0,
     *   5.59411445413800d0,   1.20983624816850d0,  3.22751022617862d0,
     *  -5.01573046424256d0,  -5.01573046424256d0,
     * -43.55555104660383d0,  -5.01573046424256d0, -5.01573046424256d0,
     *  74.48235664533756d0, -20.89534467024863d0,
     * -43.55555104660383d0,  -5.01573046424256d0, -5.01573046424256d0,
     * -20.89534467024863d0,  74.48235664533759d0,
     *   0.258504607653328389d0,  0.564050362569283648d0,     ! charges on sites O,H(1),H(2),D1(1),D1(2); Q
     *   0.564050362569283648d0, -0.693302666395947953d0,
     *  -0.693302666395947953d0,
     *   0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 15.D0,        ! omegaA, omegaB, R; COORDS: not needed in this version
     * 0.d0, -1.45365196228170d0, 1.45365196228170d0,  0.d0,       0.d0,     ! geometry of one monomer  AGEOM: not needed
     * 0.d0,  0.d0,              0.d0,        -0.2067213d0, 0.2067213d0,
     * 1.24631924913843d-1, -9.97055399310745d-1, -9.97055399310745d-1,
     *  -0.247160511d0, -0.247160511d0 /
      data nsite / 5 /
      data a0 / 0.529177249d0 /

c     compute distance matrix

      do i=1,nsite
         do j=1,nsite
            sum=0d0
            do k=1,3
               sum=sum+(coor1(k,i)-coor2(k,j))**2
            end do
            rpt(i,j)=dsqrt(sum)
c            write(*,*)rpt(i,j)
         end do
      end do

c     compute induction energy

      DO 10 I=1,5
      W(I)=0.D0
      W2(I)=0.D0
      do j=1,3
       wp(i,j) = 0.d0
       w2p(i,j) = 0.d0
       wta(i,j) = 0.d0
       w2ta(i,j) = 0.d0
       wtb(i,j) = 0.d0
       w2tb(i,j) = 0.d0
      end do
  10  CONTINUE

      DO 20 IA=1,5
      DO 20 IB=1,5
c RB implement the damping: the same as elst., distances in A !
c       dmp1 = parab(6,itypea(ia),itypeb(ib))
       dmp1 = parab(10,itypea(ia),itypeb(ib))
       d1 = d(1,dmp1,rpt(IA,IB))
c Damp only once in W*AW and W2*AW2
c       d1 = dsqrt(d1)
      W(IA)=W(IA)+a0*d1*Q(IB)/RPT(IA,IB)
      W2(IB)=W2(IB)+a0*d1*Q(IA)/RPT(IA,IB)
c
c RB calculate the force: careful with units!
c
      poma = -a0*a0*Q(IB)*(d1/(RPT(IA,IB)*RPT(IA,IB))-
     1       dmp1*dmp1*dexp(-dmp1*RPT(IA,IB)))/RPT(IA,IB)
      pomb = -a0*a0*Q(IA)*(d1/(RPT(IA,IB)*RPT(IA,IB))-
     1       dmp1*dmp1*dexp(-dmp1*RPT(IA,IB)))/RPT(IA,IB)
      do j=1,3
       r_ab(j) = (coor2(j,ib)-coor1(j,ia))
      end do
      coor2(3,ib) = coor2(3,ib) - R
      call vecp(coor1(1,ia),r_ab,vpoma)
      call vecp(r_ab,coor2(1,ib),vpomb)
      coor2(3,ib) = coor2(3,ib) + R
      do j=1,3
       wp(ia,j) = wp(ia,j) + poma*r_ab(j)
       w2p(ib,j) = w2p(ib,j) + pomb*r_ab(j)
       wta(ia,j) = wta(ia,j) + poma*vpoma(j)
       w2ta(ib,j) = w2ta(ib,j) + pomb*vpoma(j)
       wtb(ia,j) = wtb(ia,j) + poma*vpomb(j)
       w2tb(ib,j) = w2tb(ib,j) + pomb*vpomb(j)
      end do
c
  20  CONTINUE

      E=0.D0
      do j=1,3
       forceab(j) = 0.d0
       torquea(j) = 0.d0
       torqueb(j) = 0.d0
      end do 
      DO 30 I=1,5
      DO 30 J=1,5
      E=E-W(I)*A(J,I)*W(J)-W2(I)*A(J,I)*W2(J)
c Forces and torques now...
      do ia=1,3
       forceab(ia) = forceab(ia) + W(I)*A(J,I)*wp(j,ia) +
     1                            W2(I)*A(J,I)*w2p(J,ia)
       torquea(ia) = torquea(ia) + W(I)*A(J,I)*wta(j,ia) +
     1                            W2(I)*A(J,I)*w2ta(J,ia)
       torqueb(ia) = torqueb(ia) + W(I)*A(J,I)*wtb(j,ia) +
     1                            W2(I)*A(J,I)*w2tb(J,ia)
      end do
  30  CONTINUE

      E=0.5D0*627.510d0*E   ! induction energy: convert to kcal/mol
c Convert force and torque to (kcal/mol)/A
      do j=1,3
       forceab(j) = 627.510d0*forceab(j)/a0
       torquea(j) = -627.510d0*torquea(j)/a0
       torqueb(j) = -627.510d0*torqueb(j)/a0
      end do

c     E=-.16982183E-05   in this example
c      write(*,70)E
  70  format(g16.8)

      END
c
c The two-body fluctuating-charge induction model for H2O.
c Provided by Michael Smith on May 4, 2001. Modified by RB
c to connect to the fitting program ssfit. Iterative process added
c by RB on July 31, 2001.

C --------------------------------------------------------------------------

      subroutine POLDERit(coor1,coor2,E)
      implicit real*8 (a-h,o-z)
c
      parameter (maxb=500,maxp=140000)
      parameter (nsitemax=42,ntypemax=25,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=3,maxpar2=13)
      parameter (mxl=100)
c The parameters needed for damping....
      common/npar/param(maxpar1,ntypemax),
     1 parab(maxpar2,ntypemax,ntypemax),nparm,nparab
      common/types/ itypea(nsitemax),itypeb(nsitemax)
c
      DOUBLE PRECISION SUM,E
      DOUBLE PRECISION rot(3,3),coor1(3,8), coor2(3,8),geom(3,5),
     +       rpt(8,8),w(5),w2(5)
      DOUBLE PRECISION A(5,5),Q(5),COORDS(7),AGEOM(5,3)
      real*8 dqa(5), dqb(5), wp(5),w2p(5) ! induced charges and potentials
      data ((A(I,J),I=1,5),J=1,5),Q,COORDS,((AGEOM(K,L),K=1,5),L=1,3)/
     *  75.92287318493169d0,   5.59411445413799d0,  5.59411445413799d0,      ! response matrix A
     * -43.55555104660383d0, -43.55555104660384d0,
     *   5.59411445413800d0,   3.22751022617862d0,  1.20983624816850d0,
     *  -5.01573046424256d0,  -5.01573046424256d0,
     *   5.59411445413800d0,   1.20983624816850d0,  3.22751022617862d0,
     *  -5.01573046424256d0,  -5.01573046424256d0,
     * -43.55555104660383d0,  -5.01573046424256d0, -5.01573046424256d0,
     *  74.48235664533756d0, -20.89534467024863d0,
     * -43.55555104660383d0,  -5.01573046424256d0, -5.01573046424256d0,
     * -20.89534467024863d0,  74.48235664533759d0,
     *   0.258504607653328389d0,  0.564050362569283648d0,     ! charges on sites O,H(1),H(2),D1(1),D1(2); Q
     *   0.564050362569283648d0, -0.693302666395947953d0,
     *  -0.693302666395947953d0,
     *   0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 15.D0,        ! omegaA, omegaB, R; COORDS: not needed in this version
     * 0.d0, -1.45365196228170d0, 1.45365196228170d0,  0.d0,       0.d0,     ! geometry of one monomer  AGEOM: not needed
     * 0.d0,  0.d0,              0.d0,        -0.2067213d0, 0.2067213d0,
     * 1.24631924913843d-1, -9.97055399310745d-1, -9.97055399310745d-1,
     *  -0.247160511d0, -0.247160511d0 /
      data nsite / 5 /
      data a0 / 0.529177249d0 /, itmax /100/
      data dqmax / 1.d0 /
c
c RB test: print out the coordinates
c
c      do i=1,nsite
c       write(6,*)(coor1(k,i),k=1,3)
c      end do
c      do i=1,nsite
c       write(6,*)(coor2(k,i),k=1,3)
c      end do

c     compute distance matrix
c      write(6,*)'Inverse distance matrix:'
      dst_max = 1.d10
      do i=1,nsite
         do j=1,nsite
            sum=0d0
            do k=1,3
               sum=sum+(coor1(k,i)-coor2(k,j))**2
            end do
            rpt(i,j)=dsqrt(sum)
            if(rpt(i,j).lt.dst_max) dst_max = rpt(i,j) ! look for minimum
c RB implement the damping: different than elst.
c            dmp1 = parab(10,itypea(i),itypeb(j))
c             dmp1 = 10.d0
c            d1 = d(1,dmp1,rpt(i,j))
c RB damp uniformly using the O-O distance for now...
c             if(i.eq.1.and.j.eq.1) then
c               r_oo = rpt(i,j) - 1.5d0
c               d11 = d(1,dmp1,r_oo)
c             endif
            rpt(i,j)=a0/rpt(i,j)   ! convert from Angstroms
c          write(6,*)i,j,rpt(i,j)
         end do
      end do
c RB test: damp everything uniformly according to the smallest distance
      dmp1 = 4.d0
c      r_oo = dst_max - 0.8d0
      r_oo = dst_max
c      d11 = d(1,dmp1,r_oo)
      d11 = r_oo/1.5d0
      d11 = dmin1(d11,1.d0)
      do i=1,nsite
      do j=1,nsite
       rpt(i,j) = d11*rpt(i,j)
      end do
      end do
c
c     compute induction energy
c Compute permanent-charge W vectors

      DO 10 I=1,5
      W(I)=0.D0
      W2(I)=0.D0
      dqa(i) = 0.d0
      dqb(i) = 0.d0
  10  CONTINUE

      DO 20 IA=1,5
      DO 20 IB=1,5
c       dmp1 = parab(10,itypea(ia),itypeb(ib))
c       dmp1 = a0*dmp1/RPT(IA,IB)
c       d1 = d(1,dmp1,dmp1)
       d1 = 1.d0
       W(IA)=W(IA)+d1*Q(IB)*RPT(IA,IB)
       W2(IB)=W2(IB)+d1*Q(IA)*RPT(IA,IB)
  20  CONTINUE
c RB test
c      write(6,*)'Permanent W vectors:'
c      do i=1,5
c       write(6,*)i,w(i)
c       write(6,*)i,w2(i)
c      end do
c
c Iterations:
c
      IT = 0
      change = 1.d10
      eps = 1.d-16
      do while((change.gt.eps).and.(it.lt.itmax))
      IT = IT + 1
c Compute the induced W's
      do I=1,5
       WP(I)=0.D0
       W2P(I)=0.D0
      end do
      do IA=1,5
      do IB=1,5
       WP(IA)=  WP(IA)+(Q(IB)+dqb(IB))*RPT(IA,IB)
       W2P(IB)=W2P(IB)+(Q(IA)+dqa(IA))*RPT(IA,IB)
      end do
      end do
c Compute the induced charges and the induction energy
      E = 0.d0
      do I=1,5
       dqa(i) = 0.d0
       dqb(i) = 0.d0
      do J=1,5
       dqa(i) = dqa(i) - A(i,j)*WP(j)
       dqb(i) = dqb(i) - A(i,j)*W2P(j)
c      E=E-W(I)*A(J,I)*W(J)-W2(I)*A(J,I)*W2(J)
      end do
c       if(dqa(i).gt.dqmax) dqa(i) = dqmax
c       if(dqa(i).lt.-dqmax) dqa(i) = -dqmax
c       if(dqb(i).gt.dqmax) dqb(i) = dqmax
c       if(dqb(i).lt.-dqmax) dqb(i) = -dqmax
       E = E + 0.5d0*(dqa(i)*W(i) + dqb(i)*W2(i))
      end do
c      write(6,*)'Iteration',IT,627.510d0*E
c      write(6,'(5e15.7)') dqa
c      write(6,'(5e15.7)') dqb
c
      if(IT.gt.1) change = dabs(EP-E)
      EP = E
c
      end do  ! iterations loop

      E=627.510d0*E   ! induction energy: convert to kcal/mol
c      if(it.ge.itmax) write(6,*)'No convergence in POLDER!'
c
      return
      END
c
c A simple routine to calculate the dipole-dipole induction,
c for now - only in the rigid case....
c
      subroutine dipind(siteat,sitebt,oa,ob,iaa,ibb,energy)
      implicit real*8 (a-h,o-z)
cccccc Parameters taken from potparts
      parameter (maxb=500,maxp=140000)
      parameter (nsitemax=42,ntypemax=25,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=3,maxpar2=13)
      parameter (mxl=100)
c
      dimension values(1),valder(1)
      common/npar/param(maxpar1,ntypemax),
     1 parab(maxpar2,ntypemax,ntypemax),nparm,nparab
cccccc end of parameters taken from potparts
      logical first
      dimension siteat(3,1), sitebt(3,1),rot(3,3)
      dimension dma(3), dmb(3), u(3), oa(3),ob(3)
      dimension dmiut(3,0:14),polm(3,3,0:14),polis(0:14)
      data dmiu / 0.7498471537d0 /   ! dipole moment of monomer 0.
      data pol  / 10.31387108d0  /   ! dipole-dipole polarizability
      data a0 /0.529177249d0/
      data har2kcal /627.510d0/
      data first /.true./
      data dmiut /      ! monomer-spec. dipole moment vectors (MP3_resp level)
     1    0.000000000 ,  0.000000000  , -0.749847233 ,
     1    0.000000000 ,  0.000000000  , -0.726146839 ,
     1    0.000000000 ,  0.000000000  , -0.765793865,
     1    0.000000000 ,  0.000000000  , -0.893552380,
     1    0.000000000 ,  0.000000000  , -0.515479148,
     1    -0.277574780,   0.000000000 ,  -0.701285586,
     1    0.277574781 ,  0.000000000  , -0.701285586,
     1    0.000000000 ,  0.000000000  , -0.861261480,
     1    0.000000000 ,  0.000000000  , -0.504350022,
     1    0.000000000 ,  0.000000000  , -0.894120590,
     1    0.000000000 ,  0.000000000  , -0.608127530,
     1    -0.518485087,   0.000000000 ,  -0.649914180,
     1    0.518485097 ,  0.000000000  , -0.649914180,
     1    -0.158697718,   0.000000000 ,  -0.622826603,
     1    0.158697718 ,  0.000000000  , -0.622826603  /
      data polm /
c Monomer 0
     1 1.02918405E+01,   0.00000000E+00,   0.00000000E+00,
     1 0.00000000E+00,   8.30127788E+00,   0.00000000E+00,
     1 0.00000000E+00,   0.00000000E+00,   9.43156373E+00,
c Monomer 1
     1 8.20576434E+00,   0.00000000E+00,   0.00000000E+00,
     1 0.00000000E+00,   7.75092804E+00,   0.00000000E+00,
     1 0.00000000E+00,   0.00000000E+00,   8.01742765E+00,
c Monomer 2
     1 1.30170844E+01,   0.00000000E+00,   0.00000000E+00,
     1 0.00000000E+00,   8.87616704E+00,   0.00000000E+00,
     1 0.00000000E+00,   0.00000000E+00,   1.10363627E+01,
c Monomer 3
     1 9.53389234E+00,   0.00000000E+00,   0.00000000E+00,
     1 0.00000000E+00,   8.31779501E+00,   0.00000000E+00,
     1 0.00000000E+00,   0.00000000E+00,   1.01503973E+01,
c Monomer 4
     1 1.11678499E+01,   0.00000000E+00,   0.00000000E+00,
     1 0.00000000E+00,   8.42529924E+00,   0.00000000E+00,
     1 0.00000000E+00,   0.00000000E+00,   8.82460647E+00,
c Monomer 5
     1 1.10392780E+01,   0.00000000E+00,   7.59549908E-01,
     1 0.00000000E+00,   8.32436654E+00,   0.00000000E+00,
     1 7.59549908E-01,   0.00000000E+00,   8.95332446E+00,
c Monomer 6
     1 1.10392780E+01,   0.00000000E+00,   -7.59549908E-01,
     1 0.00000000E+00,   8.32436654E+00,   0.00000000E+00,
     1 -7.59549908E-01,   0.00000000E+00,   8.95332446E+00,
c Monomer 7
     1 8.74037586E+00,   0.00000000E+00,   0.00000000E+00,
     1 0.00000000E+00,   8.04628922E+00,   0.00000000E+00,
     1 0.00000000E+00,   0.00000000E+00,   9.17711640E+00,
c Monomer 8
     1 9.90796535E+00,   0.00000000E+00,   0.00000000E+00,
     1 0.00000000E+00,   8.14207837E+00,   0.00000000E+00,
     1 0.00000000E+00,   0.00000000E+00,   8.33321834E+00,
c Monomer 9
     1 1.07298657E+01,   0.00000000E+00,   0.00000000E+00,
     1 0.00000000E+00,   8.61120280E+00,   0.00000000E+00,
     1 0.00000000E+00,   0.00000000E+00,   1.11294426E+01,
c Monomer 10
     1 1.25379660E+01,   0.00000000E+00,   0.00000000E+00,
     1 0.00000000E+00,   8.70238593E+00,   0.00000000E+00,
     1 0.00000000E+00,   0.00000000E+00,   9.66615909E+00,
c Monomer 11
     1 1.10330132E+01,   0.00000000E+00,   3.17226349E-01,
     1 0.00000000E+00,   8.31415042E+00,   0.00000000E+00,
     1 3.17226349E-01,   0.00000000E+00,   8.84620107E+00,
c Monomer 12
     1 1.10330132E+01,   0.00000000E+00,   -3.17226349E-01,
     1 0.00000000E+00,   8.31415042E+00,   0.00000000E+00,
     1 -3.17226349E-01,   0.00000000E+00,   8.84620107E+00,
c Monomer 13
     1 1.11032782E+01,   0.00000000E+00,   6.04559803E-01,
     1 0.00000000E+00,   8.37154680E+00,   0.00000000E+00,
     1 6.04559803E-01,   0.00000000E+00,   8.93699552E+00,
c Monomer 14
     1 1.11032782E+01,   0.00000000E+00,   -6.04559803E-01,
     1 0.00000000E+00,   8.37154680E+00,   0.00000000E+00,
     1 -6.04559803E-01,   0.00000000E+00,   8.93699552E+00 /
      save
c
c All components of dipole-dipole polarizability are hardcoded,
c but for now we shall use only the isotropic part. In the
c beginning, this part will be calculated and stored in polis(*)
c
      if(first) then
      do i=0,14
       pom = 0.d0
       do j=1,3
        pom = pom + polm(j,j,i) 
       end do
       polis(i) = pom/3.d0
      end do
      first = .false.
      endif

c Let's place the polarizability on Oxygen, just
c for simplicity
c
c Calculate static dipole moments of A and B
c This fragment replaced by explicit rotation of
c "initial" dipole
c
      go to 10
      dlen = 0.d0
      do i=1,3     
       dma(i) = siteat(i,2)+siteat(i,3)-2.d0*siteat(i,1)
       dlen = dlen + dma(i)*dma(i)
      end do
      dlen = 1.d0/dsqrt(dlen)
      do i=1,3
       dma(i) = dmiu*dma(i)*dlen
      end do

      dlen = 0.d0
      do i=1,3
       dmb(i) = sitebt(i,2)+sitebt(i,3)-2.d0*sitebt(i,1)
       dlen = dlen + dmb(i)*dmb(i)
      end do 
      dlen = 1.d0/dsqrt(dlen)
      do i=1,3
       dmb(i) = dmiu*dmb(i)*dlen
      end do

 10   continue

      call rotmat(oa(1),oa(2),oa(3),rot)
      call matvec1(rot,dmiut(1,iaa),dma)
      call rotmat(ob(1),ob(2),ob(3),rot)
      call matvec1(rot,dmiut(1,ibb),dmb)
c
c Calculate the distance between the O's (pol centers)
c
      dlen = 0.d0
      do i=1,3
       pom = siteat(i,1) - sitebt(i,1)
       dlen = dlen + pom*pom
      end do
      dlen = dsqrt(dlen)
c Compute the damping factor for induction
      dmpind = d(6,parab(10,1,1),dlen)
cccccc
      dlen = dlen**(-3.d0)   ! inverse cube of the distance
c
c u(*): The field of A on B
c
      call Tprod(siteat(1,1),sitebt(1,1),dma,dlen,u)
c
c ...and the induction energy damped with dmpind
c
      energy = polis(iaa)*scalp(u,u)
      call Tprod(siteat(1,1),sitebt(1,1),dmb,dlen,u)
      energy = energy + polis(ibb)*scalp(u,u)
      energy = -0.5d0*(a0**6)*har2kcal*energy*dmpind
c
      return
      end
c
c Calculate the vector V resulting from the action
c of the dipole propagator tensor T_ij on U
c Ri, Rj - position vectors of molecules i and j, respectively
c rij = |Ri - Rj|^(-3)
c
        subroutine Tprod(Ri,Rj,u,rij,v)
        implicit real*8 (a-h,o-z)
        dimension Ri(3),Rj(3),u(3),v(3)
c
        ddd = rij**(0.66666666666666666d0)
        scal = 0.d0
        do i=1,3
         v(i) = Ri(i) - Rj(i)
         scal = scal + v(i)*u(i)
        end do
        do i=1,3
         v(i) = (3.d0*v(i)*scal*ddd - u(i))*rij
        end do
c
        return
        end
c
        function scalp(a,b)
        implicit real*8 (a-h,o-z)
        dimension a(3),b(3)
        scalp = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
        return
        end 
c
        subroutine rotmat(ad,bd,gd,u)
c
c Construct the rotation matrix corresponding to
c the Euler angles ad,bd,gd (in radians).
c
        implicit real*8 (a-h,o-z)
        dimension u(3,3)
        sad = dsin(ad)
        sbd = dsin(bd)
        sgd = dsin(gd)
        cad = dcos(ad)
        cbd = dcos(bd)
        cgd = dcos(gd)
c----- construct the transformation matrix
        u(1,1) = cad*cbd*cgd - sad*sgd
        u(1,2) = -cad*cbd*sgd - sad*cgd
        u(1,3) = cad*sbd
        u(2,1) = sad*cbd*cgd + cad*sgd
        u(2,2) = -sad*cbd*sgd + cad*cgd
        u(2,3) = sad*sbd
        u(3,1) = -sbd*cgd
        u(3,2) = sbd*sgd
        u(3,3) = cbd
        return
        end

c***********************************************************************
c
c Compute the iterative induction between two molecules.
c The electrostatic field supplied by routine efield.
c
      Subroutine ind_iter(siteat,nsitea,sitebt,nsiteb,R,energy)
c Assume that the dipole center is at NOT the center
c of mass, but in special point.
      implicit real*8 (a-h,o-z)
      parameter(maxsite=100)
      parameter (nsitemax=42,ntypemax=25,ntypemax2=ntypemax*ntypemax)
      Dimension G2(3,4)
      dimension ex0(2),ey0(2),ez0(2),epom(3)
      dimension siteat(3,maxsite), sitebt(3,maxsite)
      dimension Ri(3),Rj(3)
      data pol /108.9750659d0/         ! exprmt average polarizability
      data sig /0.367911875040999981d0/ ! position of polarizable center 
                                        ! same as in 3B fits
      data plen / 1.1216873242d0/
      data a0 /0.529177249d0/
      
      parameter (maxpar1=5,maxpar2=10)
      common/npar/param(maxpar1,ntypemax),
     1 parab(maxpar2,ntypemax,ntypemax),nparm,nparab
      common/types/ itypea(nsitemax),itypeb(nsitemax)
      dimension qtmp(maxsite)

      
c
c Since we are getting coordinates in bohr, we need to convert them to A first
c
c      goto 100
      do i=1,nsitea
       do j=1,3
        siteat(j,i) = siteat(j,i)/a0
       end do
      end do
      do i=1,nsiteb
       do j=1,3
        sitebt(j,i) = sitebt(j,i)/a0
       end do
      end do

100   continue

c      a03 = a0*a0*a0
       a03=1.0d0
      do i=1,2
       do j=1,3
        G2(j,i) = 0.d0
       enddo
      enddo
c Compute positions of polarizable sites and the distance between them
c      dist1G1= 1.d0/rij**3
c Ignore the rij, compute the M-M distance instead
c
c      dist1G1 = 0.d0
c      do ii=1,3
c       pom = 0.5d0*(siteat(ii,2) + siteat(ii,3))
c       Ri(ii) = siteat(ii,1) + sig*(pom-siteat(ii,1))/plen
c       pom = 0.5d0*(sitebt(ii,2) + sitebt(ii,3))
c       Rj(ii) = sitebt(ii,1) + sig*(pom-sitebt(ii,1))/plen
c       dist1G1 = dist1G1 + (Ri(ii)-Rj(ii))**2
c      end do
      dmpfct = 1.d0
c      dist1G1 = dist1G1**(-1.5d0)

       rij=R/a0

       Ri(1)=0.0d0
       Ri(2)=0.0d0
       Ri(3)=0.0d0
       Rj(1)=0.0d0
       Rj(2)=0.0d0
       Rj(3)=rij
       
c       write(*,*)"rij",rij
       
       dist1G1 = rij**(-3)
c
c Permanent charges field -- in pfld matrix
c

      do i=1,nsiteb
        qtmp(i) = param(1,itypeb(i))
c	write(*,*)sitebt(1,i),sitebt(2,i),sitebt(3,i)
        qtot=qtot+qtmp(i)	
c	write(*,*)"b",i,qtmp(i),qtot,itypeb(i)	
      end do

      call efield(Ri,qtmp,sitebt,nsiteb,epom) ! field of b on a
      ex0(1) = epom(1)
      ey0(1) = epom(2)
      ez0(1) = epom(3)
      
      
c      write(*,*)epom
      qtot=0.0d0
      do i=1,nsitea
        qtmp(i) = param(1,itypea(i))
c	write(*,*)siteat(1,i),siteat(2,i),siteat(3,i)
        qtot=qtot+qtmp(i)	
c	write(*,*)i,qtmp(i),qtot,itypea(i)	

      end do
      call efield(Rj,qtmp,siteat,nsitea,epom) ! field of a on b
      ex0(2) = epom(1)
      ey0(2) = epom(2)
      ez0(2) = epom(3)
c      write(*,*)epom

c Distance between COM's - rij

C   LETS DO IT (the iterations, that is...)

c      thr_iter = 1.d-10
      thr_iter = 1.d-30
      change = 10.d0

      isteps = 0

      do while(change.gt.thr_iter.and.isteps.lt.100)   

      energy= 0.0d0
      change= 0.0d0
c
      do ii = 1,2
       E1X = ex0(ii)
       E1Y = ey0(ii) 
       E1Z = ez0(ii) 
      do jj = 1,2
      if(ii.eq.jj) go to 111

C Calculate the total electric field at G1 of molecule ii

      call Tprod(Ri,Rj,g2(1,jj),dist1G1,epom)

      E1X = E1X + dmpfct*epom(1)*a03
      E1Y = E1Y + dmpfct*epom(2)*a03
      E1Z = E1Z + dmpfct*epom(3)*a03
      
c      write(*,*)E1X,E1Y,E1Z
c      write(*,*)epom      
c      write(*,*)g2      
c
 111  continue
      end do     ! jj loop

C Total field at ii ready. Place induced dipole....

      pole1x = pol*E1X
      pole1y = pol*E1Y
      pole1z = pol*E1Z
      change=(G2(1,ii)-polE1X)**2
     &      +(G2(2,ii)-polE1Y)**2
     &      +(G2(3,ii)-polE1Z)**2 + change

      G2(1,ii)=polE1X
      G2(2,ii)=polE1Y
      G2(3,ii)=polE1Z
   
      energy =-0.5d0*pol*(E1X*EX0(ii)
     &                 +E1Y*EY0(ii)
     &                 +E1Z*EZ0(ii)) + energy
     
c      write(*,*)ii,change,energy

      enddo     !   ii loop

       isteps = isteps +1

      enddo     !  end while

c Convert the energy in au to kcal/mol
c
      energy = 627.510d0*energy
c RB test
c      write(6,*)'isteps:',isteps
c end RB test

      return
      end
c     
c Calculate the contribution of molecule b to the field on a
c Use the CC-pol site charges.
c Units: input positions of A, 
c output fields in au
c
      subroutine efield(veci,q,sitebt,nsiteb,e)
      implicit real*8 (a-h,o-z)
      parameter (nsitemax=100)
      dimension veci(3),sitebt(3,nsiteb),e(3),chrg(8)
      dimension q(nsitemax)
      dimension sep(3,nsitemax),sepl(nsitemax)
      data a0 /0.529177249d0/
      data crgfct /18.22262373d0/

c
c Position of polarizable center of A is in veci(*)
      do isite=1,nsiteb
       if(q(isite).ne.0.d0) then   ! if charged site
       sepl(isite) = 0.d0
       do k=1,3
        sep(k,isite) = veci(k) - sitebt(k,isite)
        sepl(isite) = sepl(isite) + sep(k,isite)*sep(k,isite)
       end do
       sepl(isite) = sepl(isite)**(-1.5d0)
       endif                            ! if charged site
      end do

      do k=1,3
       e(k) = 0.d0
      end do

      totchrg=0.0d0
      do isite=1,nsiteb
       if(q(isite).ne.0.d0) then   ! if charged site
         totchrg=totchrg+q(isite)/crgfct
c	 write(*,*)isite,q(isite)/crgfct,totchrg
        do k=1,3
	  e(k) = e(k) + q(isite)*sep(k,isite)*sepl(isite)/crgfct
c         e(k) = e(k) + a0*a0*q(isite)*sep(k,isite)*sepl(isite)
        end do
       endif      ! if charged site
      end do
c      write(*,*)"totchrg",totchrg
c
      return
      end
