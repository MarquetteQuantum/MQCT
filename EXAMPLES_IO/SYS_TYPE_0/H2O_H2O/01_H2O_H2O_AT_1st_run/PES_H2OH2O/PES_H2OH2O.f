!----------!
! OPTION 1 !
!----------! 
!	  THIS IS THE DEFAULT OPTION IN THE CODE, NO KEYWORDS ARE REQUIRED 
      SUBROUTINE USER_DEFINED_PES(V,R,rvib,rvib2,alpha,beta,gamma,
     &										   aalpha,bbeta,ggamma)
!     INPUT:  R - distance betweenn COMs, rvib - vibrational coordinate
! 		      alpha,beta,gamma - Euler's angles of the first molecule
!   		  aalpha,bbeta, ggamma - Euler's angles of the second molecule
!     OUTPUT: V - value of the potential
      IMPLICIT NONE
      REAL*8 V,R,alpha,beta,gamma,aalpha,bbeta,ggamma,t,rvib,rvib2
	  CALL V_POT_H2O_H2O(V,R,alpha,beta,gamma,aalpha,bbeta,ggamma)
      RETURN	  
!     USER MUST INSERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE IS SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:
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
 
      MODULE POT_DATA_H2O_H2O   
      LOGICAL :: INI = .FALSE.
      INTEGER isurf,iemonomer,iembedang,ixz	
      REAL*8 omegA(3),omegB(3),RR_AB,vvalA(3),vvalB(3),Etot
      END MODULE POT_DATA_H2O_H2O

      SUBROUTINE V_POT_H2O_H2O(V,R,a1,b1,g1,a2,b2,g2) !!! ASYMETRC + ASYMETRIC
      USE POT_DATA_H2O_H2O
      IMPLICIT NONE	  
      REAL*8 V,a1,b1,g1,a2,b2,g2,R
      REAL*8 rHO,angHOH
      rHO = 0.9584d0/0.52917721067
      angHOH = 104.45d0/180d0*dacos(-1d0)	  
      IF(.not.INI) THEN
      isurf = 1
      iemonomer=0
      iembedang=1
      ixz=1
      vvalA(1) = rHO
      vvalB(1) = rHO
      vvalA(2) = rHO
      vvalB(2) = rHO
      vvalA(3) = angHOH
      vvalB(3) = angHOH	  
      CALL init_Potential(isurf,iemonomer,iembedang,ixz)
      INI = .TRUE.  
      ENDIF
      RR_AB = R 
      omegA(1) = a1
      omegA(2) = b1
      omegA(3) = g1
      omegB(1) = a2
      omegB(2) = b2
      omegB(3) = g2	  
      CALL Potential_ang(Etot,omegA,omegB,RR_AB,vvalA,vvalB)
      V = Etot
      END SUBROUTINE V_POT_H2O_H2O	  
	  
      subroutine init_Potential(isurf,iemon,iembedang,ixyz)
      implicit real*8 (a-h,o-z)
      dimension Oaa(3),Haa1(3),Haa2(3),Obb(3),Hbb1(3),Hbb2(3)
      dimension carta(3,3),cartb(3,3)
      character*50 nameofsaptdatafile
      character(len=11) :: path = "PES_H2OH2O/"

      common /ddaattaa/ param(18,6),parab(84,6,6),
     .       c(1000),cc(2000),params(1000),
     .       chrg(25),sites(3,25),ind_data(6250),
     .       nsitea,nsiteb,nparsall,iembed,ipotparts,icc,
     .       iemonomer,iembedinterf,ixz
	  open(6,file=path//'PES_H2O_H2O_Output.out',access='append')
      write(6,'(a)') "-----------------------------------------"
c one can choose from the following list of full-dimensional surfaces
c (the definitions are given in the paper):
c     isurf          short description:            
c       1            CCpol-8sf[2014] (Eckart embedding)
c       2            CCpol-8sfIR[2014] (Eckart embedding)
c       3            CCpol-8sfIR[2012] (Radau f=1 embedding)
c       4            CCpol-8sfIR[2012] (Eckart embedding)
c       5            SAPT-5s'f[2014]         ! 
c       6            SAPT-5s'fIR[2014]       !
c       7            SAPT-5s'fIR[2012]       !
c       8            SAPT-5s'f[2006]         !
c       9            SAPT-5s'fIR[2006]       !
c      10            CCpol-8sfIR[2014] (Radau f=1 embedding)

      if (isurf.eq.1) then
        write(6,*) 
     .      "generating CCpol-8sf[2014] (Eckart embedding) energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='data_SAPT5spf_2014'
        icc=1

      elseif (isurf.eq.2)then
        write(6,*) 
     .      "generating CCpol-8sfIR[2014] (Eckart embedding) energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='data_SAPT5spfIR_2014'
        icc=1

      elseif (isurf.eq.3)then
        write(6,*) 
     .    "generating CCpol-8sfIR[2012] (Radau f=1 embedding) energies"
        iembed=2
        ipotparts=1
        nameofsaptdatafile='data_SAPT5spfIR_2006'
        icc=1

      elseif (isurf.eq.4)then
        write(6,*) 
     .    "generating CCpol-8sfIR[2012] (Eckart embedding) energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='data_SAPT5spfIR_2006'
        icc=1

      elseif (isurf.eq.5)then
        write(6,*) "generating SAPT-5s'f[2014] energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='data_SAPT5spf_2014'
        icc=0

      elseif (isurf.eq.6)then
        write(6,*) "generating SAPT-5s'fIR[2014] energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='data_SAPT5spfIR_2014'
        icc=0

      elseif (isurf.eq.7)then
        write(6,*) "generating SAPT-5s'fIR[2012] energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='data_SAPT5spfIR_2006'
        icc=0

      elseif (isurf.eq.8)then
        write(6,*) "generating SAPT-5s'f[2006] energies"
        iembed=1
        ipotparts=0
        nameofsaptdatafile='data_SAPT5spf_2006'
        icc=0

      elseif (isurf.eq.9)then
        write(6,*) "generating SAPT-5s'fIR[2006] energies"
        iembed=1
        ipotparts=0
        nameofsaptdatafile='data_SAPT5spfIR_2006'
        icc=0

      elseif (isurf.eq.10)then
        write(6,*) 
     .     "generating CCpol-8sfIR[2014] (Radau f=1 embedding) energies"
        iembed=2
        ipotparts=1
        nameofsaptdatafile='data_SAPT5spfIR_2014'
        icc=1

      else
        write(6,*) "wrong value of isurf: ",isurf
        stop
      endif

c      iembed=2       ! 1 for Eckart     <--------!!!!!!!!!!!
c                     ! 2 for Radau f=1  <--------!!!!!!!!!!!

      pi=dacos(-1.d0)
      rad2deg=180.d0/pi

c..... initialize the interaction energy procedures .....
      open(55,file=path//nameofsaptdatafile,status='unknown')
      imode = -1
      call poten(r1,r2,theta1,r3,r4,theta2,
     .      rin,alphaAv,betaAv,gammaAv,alphaBv,betaBv,gammaBv,val,imode,
     .      icart,carta,cartb,ipotparts,
     .      param,parab,nsitea,nsiteb,c)
      close(55)

      imode = -1
      call ccpol8s_dimer(imode,Oaa,Haa1,Haa2,Obb,Hbb1,Hbb2,Erigid,
     .                   cc,params,nparsall,
     .                   ind_data,chrg,sites)
c........................................................

      iemonomer=iemon

      if (iembedang.ne.0) then
      iembedinterf=iembedang
      ixz=ixyz
      if (iembedang.eq.1) then
        write(6,"(a,a)") " Eckart embedding is applied to translate",
     .             " (angles + bonds) to the Cartesian coordinates"
      else
       if (iembedang.eq.2) then
        write(6,"(a,a)") " Radau f=1 embedding is applied to translate",
     .                 " (angles + bonds) to the Cartesian coordinates"
       else
          write(6,*) "wrong value of iembedang"
          stop
        endif
      endif

      if (ixyz.eq.1) then
        write(6,*) "initially the molecule is put in the xz plane"
      else
        if (ixyz.eq.0) then
          write(6,*) "initially the molecule is put in the yz plane"
        else
          write(6,*) "wrong value of ixz"
          stop
        endif
      endif
      endif ! <---- if (iembedang.eq.0) then

      if (iemon.eq.1) then
        write(6,*) 
     .    "energy of the monomers is added to the interaction energy"
      else
        if (iemon.ne.0) then
          write(6,*) "wrong value of iemonomer"
          stop
        endif
      endif
      write(6,'(a)') "-----------------------------------------"

      call test_parameters(isurf,iemon)

      return
      end
c
c------------------------------------------------------------------------------
      subroutine test_parameters(isurf,iemon)
      implicit real*8 (a-h,o-z)
      dimension val(10),valm(10)
      dimension Oa(3),Ha1(3),Ha2(3),Ob(3),Hb1(3),Hb2(3)
      data val / -1.45754, -1.36660, -1.34676, -1.34688, -0.76978,
     .           -0.67884, -0.65912, -0.75268, -0.65451, -1.36648/
      data valm/ 25.04582, 25.13676, 25.15660, 25.15648, 25.73358,
     .           25.82452, 25.84424, 25.75068, 25.84885, 25.13688/
      data Oa  /  0.6458557220D-01, 0.3399054992D-02,-0.1782922818D-01/
      data Ha1 / -0.5505396894D+00, 0.6383283738D-01,-0.6241648475D+00/
      data Ha2 / -0.4744802670D+00,-0.1177783097D+00, 0.9071276528D+00/
      data Ob  / -0.5658499752D-01,-0.1827211353D-03, 0.2523332441D+01/
      data Hb1 /  0.6334615998D+00, 0.2269642816D+00, 0.2055825014D+01/
      data Hb2 /  0.2645834235D+00,-0.2240643644D+00, 0.3293430563D+01/

      call Potential(Etot,Oa,Ha1,Ha2,Ob,Hb1,Hb2)

      write(6,'(a)') "to check if the parameters are read properly:"
      if (iemon.eq.0) then
        write(6,'(2(a,f10.5))') 
     .      "Calculated value: ",Etot,"  and should be:",val(isurf)
      else
        write(6,'(2(a,f10.5))') 
     .      "Calculated value: ",Etot,"  and should be:",valm(isurf)
      endif
      write(6,'(a)') "-----------------------------------------"

      return
      end
c
c------------------------------------------------------------------------------
c
c The input Cartesian coordinates should be given in Angstroms.
c
      subroutine Potential(Etot,Oa,Ha1,Ha2,Ob,HB1,Hb2)
      implicit real*8 (a-h,o-z)
      dimension Oa(3),Ha1(3),Ha2(3),Ob(3),Hb1(3),Hb2(3)

      common /ddaattaa/ param(18,6),parab(84,6,6),
     .       c(1000),cc(2000),params(1000),
     .       chrg(25),sites(3,25),ind_data(6250),
     .       nsitea,nsiteb,nparsall,iembed,ipotparts,icc,
     .       iemonomer,iembedinterf,ixz

      call CCpol_xyz
     .      (nnr,Oa,Ha1,Ha2,Ob,Hb1,Hb2,iembed,ipotparts,icc,Etot,
     .      param,parab,nsitea,nsiteb,c,
     .      cc,params,nparsall,
     .      ind_data,chrg,sites)

      if (iemonomer.eq.1) then
        rA1=0.0d0
        rA2=0.0d0
        rB1=0.0d0
        rB2=0.0d0
        ssA=0.0d0
        ssB=0.0d0
        do j=1,3
          rA1=rA1+(Ha1(j)-Oa(j))**2
          rA2=rA2+(Ha2(j)-Oa(j))**2
          rB1=rB1+(Hb1(j)-Ob(j))**2
          rB2=rB2+(Hb2(j)-Ob(j))**2
          ssA=ssA+(Ha1(j)-Oa(j))*(Ha2(j)-Oa(j))
          ssB=ssB+(Hb1(j)-Ob(j))*(Hb2(j)-Ob(j))
        enddo
        rA1=dsqrt(rA1)
        rA2=dsqrt(rA2)
        rB1=dsqrt(rB1)
        rB2=dsqrt(rB2)
        thA=dacos(ssA/(rA1*rA2))
        thB=dacos(ssB/(rB1*rB2))

c        pi=dacos(-1.d0)
c        rad2deg=180.d0/pi
        a0=0.529177249d0
        h2kcal=627.510d0

c        write(6,'(a,3f13.8)') "rA1, rA2, thA:  ",rA1, rA2, thA*rad2deg
c        write(6,'(a,3f13.8)') "rB1, rB2, thB:  ",rB1, rB2, thB*rad2deg
        rA1=rA1/a0
        rA2=rA2/a0
        rB1=rB1/a0
        rB2=rB2/a0
        call pots(vA,rA1,rA2,thA)
        call pots(vB,rB1,rB2,thB)
c        write(6,'(a,f13.8)') "vA= ",vA*h2kcal
c        write(6,'(a,f13.8)') "vB= ",vB*h2kcal

c        write(6,'(a,f13.8)') "Etot(dimer only): ",Etot
        Etot=Etot+(vA+vB)*h2kcal

      endif

      return
      end
c
c -----------------------------------------------------------------------------
      Subroutine CCpol_xyz(nnr,Oa,Ha1,Ha2,Ob,HB1,Hb2,
     .                 iembed,ipotparts,icc,Etot,
     .                 param,parab,nsitea,nsiteb,c,
     .                 cc,params,nparsall,
     .                 ind_data,chrg,sites)
      implicit real*8 (a-h,o-z)
      dimension Oa(3),Ha1(3),Ha2(3),COMa(3)
      dimension Ob(3),Hb1(3),Hb2(3),COMb(3)
      dimension Oaa(3),Haa1(3),Haa2(3)
      dimension Obb(3),Hbb1(3),Hbb2(3)
      character*1 c1
      dimension vecI(3),vecIa(3),vecIb(3),vecJa(3),vecJb(3)
      dimension carta(3,3),cartb(3,3)
      dimension cartaa(3,3),cartbb(3,3)
      character*50 nameofsaptdatafile
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=18,maxpar2=84)
c
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension c(1000)
      PARAMETER (maxlin=2000,mmaxpar=1000)
      dimension cc(maxlin),params(mmaxpar)

      dimension chrg(25)
      dimension sites(3,25)
      dimension ind_data(6250)
      integer omp_get_num_threads,omp_get_thread_num

c      ns=omp_get_num_threads()
c      id=omp_get_thread_num()

      a0=0.529177249d0
      vecI(1)=0.0d0
      vecI(2)=0.0d0
      vecI(3)=-1.0d0

      imode = 0

        call align_on_z_axis(Oa,Ha1,Ha2,Ob,Hb1,Hb2,Rcom)

        do jj=1,3      
          carta(1,jj)=Oa(jj)           
          carta(2,jj)=Ha1(jj)           
          carta(3,jj)=Ha2(jj)           
          cartb(1,jj)=Ob(jj)           
          cartb(2,jj)=Hb1(jj)           
          cartb(3,jj)=Hb2(jj)           
        enddo

c make embedding for molecule A
        if (iembed.eq.1) then
          call eck_rad_tst(Oa,Ha1,Ha2,vecIa,vecJa)
        endif
        if (iembed.eq.2) then
          call radau_f1_tst(Oa,Ha1,Ha2,vecIa,vecJa)
        endif

        call put_rigid(vecIa,vecJa,Oaa,Haa1,Haa2)

c make embedding for molecule B
        Ob(3)=Ob(3)-Rcom
        Hb1(3)=Hb1(3)-Rcom
        Hb2(3)=Hb2(3)-Rcom

        if (iembed.eq.1) then
          call eck_rad_tst(Ob,Hb1,Hb2,vecIb,vecJb)
        endif
        if (iembed.eq.2) then
          call radau_f1_tst(Ob,Hb1,Hb2,vecIb,vecJb)
        endif

        Ob(3)=Ob(3)+Rcom
        Hb1(3)=Hb1(3)+Rcom
        Hb2(3)=Hb2(3)+Rcom

        call put_rigid(vecIb,vecJb,Obb,Hbb1,Hbb2)
        Obb(3)=Obb(3)+Rcom
        Hbb1(3)=Hbb1(3)+Rcom
        Hbb2(3)=Hbb2(3)+Rcom

        do jj=1,3          
          cartaa(1,jj)=Oaa(jj)           
          cartaa(2,jj)=Haa1(jj)           
          cartaa(3,jj)=Haa2(jj)           
          cartbb(1,jj)=Obb(jj)           
          cartbb(2,jj)=Hbb1(jj)           
          cartbb(3,jj)=Hbb2(jj)           
        enddo

c        icc=1
c        icc=0
        if (icc.eq.1) then
         call driver_potss_sapt5sf(carta,cartb,val,
     .                     ipotparts,param,parab,nsitea,nsiteb,c)
         call driver_potss_sapt5sf(cartaa,cartbb,vall,
     .                     ipotparts,param,parab,nsitea,nsiteb,c)
         call ccpol8s_dimer(imode,Oaa,Haa1,Haa2,Obb,Hbb1,Hbb2,Erigid,
     .                     cc,params,nparsall,
     .                     ind_data,chrg,sites)
         Etot= Erigid + (val-vall)
        else 
         call driver_potss_sapt5sf(carta,cartb,val,
     .                     ipotparts,param,parab,nsitea,nsiteb,c)
         Etot= val
        endif

      return
      end
c------------------------------------------------------------------
c The rigid (reference) molecule is embedded into the plane defined by the original
c (distorted) molecule and by the axis I related to the embedding method. The bisection
c vector of the rigid molecule is put along the I vector and COM of the rigid molecule
c is put at the position of COM of the distorted one.
c
c Information about two normalized vectors is passed to the procedure.
c vi1 - the bisection vector of the reference water is put in this direction
c vi2 - it is perpendicular to vi1 and lies in the plane of the molecule
c
      subroutine put_rigid (vi1,vi2,O,H1,H2)
      implicit real*8 (a-h,o-z)
      dimension vi1(3),vi2(3)
      dimension O(3),H1(3),H2(3)
      dimension w1(3),w2(3),vshift(3)
      dimension Oposition(3),COM(3)
c ds=sin(th_ref_half), dc=cos(th_ref_half)
      ds=0.79170358110560535d0
      dc=0.61090542612139243d0
      rOHref=0.97162570027717354d0
      com_shift=0.66429466101803d-01
      a0=0.529177249d0

      pi=dacos(-1.d0)
      rad2deg=180.d0/pi

c Initially O is put in the origin, but next the whole rigid molecule should be shifted
c to put its COM in the origin. This shif is in the vi1 direction by com_shift.
      Oposition(1)=0.0d0
      Oposition(2)=0.0d0
      Oposition(3)=0.0d0

      do j=1,3
        w1(j)=dc*vi1(j)+ds*vi2(j)
        w2(j)=dc*vi1(j)-ds*vi2(j)
        vshift(j)=-com_shift*vi1(j)
      enddo
      do j=1,3
        w1(j)=rOHref*w1(j)
        w2(j)=rOHref*w2(j)
      enddo

      do j=1,3
        w1(j)=w1(j)+vshift(j)
        w2(j)=w2(j)+vshift(j)
        Oposition(j)=Oposition(j)+vshift(j)
      enddo
      do j=1,3
        O(j)=Oposition(j)
        H1(j)=w1(j)
        H2(j)=w2(j)
      enddo

      return
      end
c------------------------------------------------------------------
c This procedure is based on the Konrad's antigetgeo.py program.
c It transforms the cartesian coordinates of two water molecules given in the following
c order (OA, HA1, HA2, OB, HB1, HB2), into cartesian coordinates but in the system, where
c 1. COM of the first molecule is at the (0,0,0) position
c 2. COM of the second one is at (0,0,R), where R is the distance between COMs
c
      subroutine align_on_z_axis(O1A,H1A,H2A,O1B,H1B,H2B,Rcom)
      implicit real*8 (a-h,o-z)
      dimension xyzA(3,3), xyzB(3,3)
      dimension xyzAA(3,3), xyzBB(3,3)
      dimension xyzAAA(3,3), xyzBBB(3,3)
      dimension ama(3),comA(3),comB(3),s(3),s1(3)
      dimension O1A(3),H1A(3),H2A(3),O1B(3),H1B(3),H2B(3)
      character*1 c1
c      data ama/15.994915d0, 1.007825d0, 1.007825d0/
c      save ama

      bohr2a=0.529177d0
      ama(1)=15.994915d0
      ama(2)=1.007825d0 
      ama(3)=1.007825d0 

c read the cartesian coordiantes
      thr=1.0d-9

      do j=1,3
        xyzA(1,j)=O1A(j)
        xyzA(2,j)=H1A(j)
        xyzA(3,j)=H2A(j)
        xyzB(1,j)=O1B(j)
        xyzB(2,j)=H1B(j)
        xyzB(3,j)=H2B(j)
      enddo

c test a procedure
      call COMcalc3 (xyzA,comA)
      call COMcalc3 (xyzB,comB)

      sss=0.0d0
      do j=1,3
        sss=sss+(comB(j)-comA(j))**2
      enddo
      Rcom=dsqrt(sss)

c shift the whole system in such a way that comA coincides with the origin of the 
c system of coordinates
      do i=1,3
        do j=1,3
          xyzA(i,j)=xyzA(i,j)-comA(j)
          xyzB(i,j)=xyzB(i,j)-comA(j)
        enddo
        comB(i)=comB(i)-comA(i)
      enddo
      do i=1,3
        comA(i)=comA(i)-comA(i)
      enddo

      ss=dsqrt(comB(1)*comB(1)+comB(2)*comB(2))
      if (ss.lt.thr) then
c        write(66,*) "    both center of masses are on the Z axis, ss=",ss
c "Konrad problem"
c Even if both are on the Z axis, take care about position of comB with respect to comA: 
c the value of the z-coordinate of B should be larger than of A (means 0 in our case)
c Make an inversion when needed
        if (comB(3).lt.0.0d0) then
c          write(66,*) "I have to flip the complex"
          do i=1,3
            do j=1,3
              xyzA(i,j)=-xyzA(i,j)
              xyzB(i,j)=-xyzB(i,j)
            enddo
            comA(i)=-comA(i)
            comB(i)=-comB(i)
          enddo
        endif
        do i=1,3
          do j=1,3
            xyzAAA(i,j)=xyzA(i,j)
            xyzBBB(i,j)=xyzB(i,j)
          enddo
        enddo
      else
c        write(66,*) "not both center of masses are on Z axis, ss=",ss
c find a vector s (rotation axis) which is perpendicular both to the vector connecting
c comA and comB, and to z-axis.
c s = comB x (0,0,1) = (comB(2),-comB(1),0)
c + normalization
        xnorm=dsqrt(comB(1)*comB(1)+comB(2)*comB(2))
        s(1)=comB(2)/xnorm
        s(2)=-comB(1)/xnorm
        s(3)=0.0d0
c find an angle between comB and (0,0,1)
        rr=dsqrt(comB(1)*comB(1)+comB(2)*comB(2)+comB(3)*comB(3))
        ccos=comB(3)/rr
        ssin=dsqrt(1.0d0-ccos*ccos)
c now we have a new orthogonal basis {s,s',(0,0,01)}, where
c s'=(-comB(1)/sqrt(comB(1)**2+comB(2)**2),-comB(2)/sqrt(comB(1)**2+comB(2)**2),0)
        s1(1)=-comB(1)/xnorm
        s1(2)=-comB(2)/xnorm
        s1(3)=0.0d0

c let's express the positions of atoms in the new basis:
        do i=1,3
          xyzAA(i,1)=xyzA(i,1)*s(1)+xyzA(i,2)*s(2)
          xyzAA(i,2)=xyzA(i,1)*s1(1)+xyzA(i,2)*s1(2)
          xyzAA(i,3)=xyzA(i,3)
        enddo
        do i=1,3
          xyzBB(i,1)=xyzB(i,1)*s(1)+xyzB(i,2)*s(2)
          xyzBB(i,2)=xyzB(i,1)*s1(1)+xyzB(i,2)*s1(2)
          xyzBB(i,3)=xyzB(i,3)
        enddo

c rotate everything:
        do i=1,3
          xyzAAA(i,1)= xyzAA(i,1)
          xyzAAA(i,2)= xyzAA(i,2)*ccos+xyzAA(i,3)*ssin
          xyzAAA(i,3)=-xyzAA(i,2)*ssin+xyzAA(i,3)*ccos
        enddo
        do i=1,3
          xyzBBB(i,1)= xyzBB(i,1)
          xyzBBB(i,2)= xyzBB(i,2)*ccos+xyzBB(i,3)*ssin
          xyzBBB(i,3)=-xyzBB(i,2)*ssin+xyzBB(i,3)*ccos
        enddo
      endif

      do j=1,3
        O1A(j)=xyzAAA(1,j)
        H1A(j)=xyzAAA(2,j)
        H2A(j)=xyzAAA(3,j)
        O1B(j)=xyzBBB(1,j)
        H1B(j)=xyzBBB(2,j)
        H2B(j)=xyzBBB(3,j)
      enddo

      return
      end
c-------------------------------------------------------------------------------------
      subroutine COMcalc3 (xyz,COM)
      implicit real*8 (a-h,o-z)
      dimension xyz(3,3), xyzB(3,3)
      dimension O1(3),H1(3),H2(3),COM(3)
      real*8 mO,mH,M
      mO=15.994 914 6221d0
      mH= 1.007 825 032 1d0

      do j=1,3
        O1(j)=xyz(1,j)
        H1(j)=xyz(2,j)
        H2(j)=xyz(3,j)
      enddo
      M=mO+mH+mH
      do i=1,3
        COM(i)= (mO*O1(i) +mH*H1(i) +mH*H2(i))/M
      enddo 
      end
c-------------------------------------------------------------------------------------
c
c I assume that the molecule is already shifted to fit COM to (0,0,0) 
c 
      subroutine eck_rad_tst(r0,r1,r2,vecI,vecJ)
      implicit real*8 (a-h,o-z)
      dimension rot(3,3)
      dimension q1(3),q2(3),bv(3),bbv(3)
      dimension rcom(3),r0(3),r1(3),r2(3)
      dimension br0(3)
      dimension pom1(3),pom2(3)
      dimension temp1(3),temp2(3),vecI(3),vecJ(3)
      character*1 c1
c      data xmO /15.994915d0/, xmH /1.007825d0/
      xmO=15.9949146221d0
      xmH=1.0078250321d0
      xq1e=0.95111822d0
      xq2e=0.95111822d0
      theta_r_e_deg=107.952612d0
      theta_r_e=1.88412851d0
      a0=0.529177249d0

      deg2rad=dacos(-1.d0)/180.d0

c To perform the Eckart embedding we need info about the Radau coordinates of 
c the molecule in its equilibrium geometry. It was precalculated with the 
c data given above and hard-coded as xq1e, xq2e, and theta_r_e

      xm12=2*xmH
      xm=xm12+xmO
      xm0=xmO
      xm1=xmH
      xm2=xmH
      alpha=dsqrt(xm0/xm)               ! Eq. (26)    can be precalculated
      b=(alpha-alpha*alpha)*xm/xm12     ! Eq. (27)

c move the COM of the molecule to the origin of the coordinate system
      do j=1,3
        q1(j)=r1(j)-b*r0(j)    ! Eq. (28)
        q2(j)=r2(j)-b*r0(j)    ! Eq. (29)
        br0(j)=b*r0(j)
      enddo

c calculate length of q1 and q2 vectors, and the angle between them
        xq1=0.0d0
        xq2=0.0d0
        sss=0.0d0
        do j=1,3
          xq1=xq1+q1(j)*q1(j)
          xq2=xq2+q2(j)*q2(j)
          sss=sss+q1(j)*q2(j)
        enddo
        xq1=dsqrt(xq1)
        xq2=dsqrt(xq2)
        theta_r=dacos( sss/(xq1*xq2) )

        sss=0.0d0
        do j=1,3
          pom1(j)=q1(j)/xq1
          pom2(j)=q2(j)/xq2
          bv(j)=pom1(j)+pom2(j)
          sss=sss+bv(j)*bv(j)
        enddo
        sss=dsqrt(sss)
        do j=1,3
          bv(j)=bv(j)/sss
        enddo

c formula (33) from Wei & Carrington, JCP 107, 2813 (1997)
c masses are the same and can be skipped, in fact xq1e and xq2e could be also
        eta_e=0.5d0*theta_r_e
        ang=theta_r-theta_r_e+eta_e
        sss=(xq2e*xq2*dsin(ang) + xq1e*xq1*dsin(eta_e))/
     .      (xq2e*xq2*dcos(ang) + xq1e*xq1*dcos(eta_e))
        eta=datan(sss)

c To obtain the I vector which bisects the angle made by q1e and q2e I do not
c need an explicit information about these vectors, but I can get it by the
c rotation  of the q1 vector in direction of q2 by the angle equal to eta.
c Temporary axes temp1=norm(q1) and temp2=norm[q2-((q1*q2)/(q1*q1)) q1]
c                                or temp2=norm[q2-(temp1*q2) temp1]

        sss=0.0d0
        do j=1,3
          temp1(j)=q1(j)/xq1
          sss=sss+temp1(j)*q2(j)
        enddo
        ttt=0.0d0
        do j=1,3
          temp2(j)=q2(j)-sss*temp1(j)  
          ttt=ttt+temp2(j)**2
        enddo
        ttt=dsqrt(ttt)
        do j=1,3
          temp2(j)=temp2(j)/ttt
        enddo
        
c        write(67,666) (temp1(j),j=1,3)," temp1"
c        write(67,666) (temp2(j),j=1,3)," temp2"
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        do j=1,3
          s1=s1+temp1(j)*temp1(j)
          s2=s2+temp2(j)*temp2(j)
          s3=s3+temp1(j)*temp2(j)
        enddo
c        write(67,'(a,3f12.6)') "s1, s2, s3:", s1, s2, s3
c        write(67,*)
c vecI=cos(eta)*temp1 + sin(eta)*temp2
        do j=1,3
          vecI(j)= dcos(eta)*temp1(j)+dsin(eta)*temp2(j)
          vecJ(j)=-dsin(eta)*temp1(j)+dcos(eta)*temp2(j)
        enddo

      return

  666 format(3f12.6,3x,a)
  777 format(a,2f12.8,f12.6,f12.8)
  778 format(a,f12.6,f12.8)

      end
c-----------------------------------------------------------------------------------
c
c I assume that the molecule is already shifted to fit COM to (0,0,0) 
c 
      subroutine radau_f1_tst(r0,r1,r2,vecI,vecJ)
      implicit real*8 (a-h,o-z)
      dimension rot(3,3)
      dimension q1(3),q2(3),bv(3),bbv(3)
      dimension rcom(3),r0(3),r1(3),r2(3)
      dimension br0(3)
      dimension pom1(3),pom2(3)
      dimension temp1(3),temp2(3),vecI(3),vecJ(3)
      character*1 c1
c      data xmO /15.994915d0/, xmH /1.007825d0/
      xmO=15.9949146221d0
      xmH=1.0078250321d0
      xq1e=0.95111822d0
      xq2e=0.95111822d0
      theta_r_e_deg=107.952612d0
      theta_r_e=1.88412851d0
      a0=0.529177249d0

      deg2rad=dacos(-1.d0)/180.d0

c data given above and hard-coded as xq1e, xq2e, and theta_r_e


      xm12=2*xmH
      xm=xm12+xmO
      xm0=xmO
      xm1=xmH
      xm2=xmH
      alpha=dsqrt(xm0/xm)               ! Eq. (26)    can be precalculated
      b=(alpha-alpha*alpha)*xm/xm12     ! Eq. (27)

c move the COM of the molecule to the origin of the coordinate system
      do j=1,3
        q1(j)=r1(j)-b*r0(j)    ! Eq. (28)
        q2(j)=r2(j)-b*r0(j)    ! Eq. (29)
        br0(j)=b*r0(j)
      enddo

c calculate length of q1 and q2 vectors, and the angle between them
        xq1=0.0d0
        xq2=0.0d0
        sss=0.0d0
        do j=1,3
          xq1=xq1+q1(j)*q1(j)
          xq2=xq2+q2(j)*q2(j)
          sss=sss+q1(j)*q2(j)
        enddo
        xq1=dsqrt(xq1)
        xq2=dsqrt(xq2)
        theta_r=dacos( sss/(xq1*xq2) )

        sss=0.0d0
        do j=1,3
          pom1(j)=q1(j)/xq1
          pom2(j)=q2(j)/xq2
          bv(j)=pom1(j)+pom2(j)
          sss=sss+bv(j)*bv(j)
        enddo
        sss=dsqrt(sss)
        do j=1,3
          bv(j)=bv(j)/sss
          vecI(j)=bv(j)
        enddo

c I need also a vector which is perpendicular to vecI and lies in the plane
c of the molecule.
c I can obtain it from the following expression:
c temp2=norm[q2-(vecI*q2) vecI]     (q2 or q1 - does not matter)

        sss=0.0d0
        do j=1,3
          sss=sss+vecI(j)*q2(j)
        enddo
        ttt=0.0d0
        do j=1,3
          temp2(j)=q2(j)-sss*vecI(j)
          ttt=ttt+temp2(j)**2
        enddo
        ttt=dsqrt(ttt)
        do j=1,3
          temp2(j)=temp2(j)/ttt
c          vecJ(j)=temp2(j)          mozna wybrac zwrot
          vecJ(j)=-temp2(j)
        enddo

      return

  666 format(3f12.6,3x,a)
  777 format(a,2f12.8,f12.6,f12.8)
  778 format(a,f12.6,f12.8)

      end
c-----------------------------------------------------------------------------------
      subroutine vecprod(v1,v2,v3)
      implicit real*8 (a-h,o-z)
      dimension v1(3), v2(3), v3(3)

      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      return
      end
c-----------------------------------------------------------------------
      subroutine read_cc_data(ind_data,chrg,sites)
      implicit real*8 (a-h,o-z)
      character(len=11) :: path = "PES_H2OH2O/"
      character*1 cnull
      dimension ind_charge(25)
      dimension ind_beta(25,25)
      dimension ind_d1(25,25)
      dimension ind_d6(25,25)
      dimension ind_d8(25,25)
      dimension ind_d10(25,25)
      dimension ind_C6(25,25)
      dimension ind_C8(25,25)
      dimension ind_C10(25,25)
      dimension chrg(25)
      dimension sites(3,25)

      dimension ind_data(6250)
      dimension ind_data1(6250)
      equivalence (ind_data1(1),ind_charge)
      equivalence (ind_data1(26),ind_beta)
      equivalence (ind_data1(651),ind_d1)
      equivalence (ind_data1(1276),ind_d6)
      equivalence (ind_data1(1901),ind_d8)
      equivalence (ind_data1(2526),ind_d10)
      equivalence (ind_data1(3151),ind_c6)
      equivalence (ind_data1(3776),ind_c8)
      equivalence (ind_data1(4401),ind_c10)

      do i=1,6250
        ind_data1(i)=0
      enddo

      do i=1,25
       ind_charge(i)=0
       chrg(i)=0.0d0
       do j=1,25
        ind_beta(i,j)=0
        ind_d1(i,j)=0
        ind_d6(i,j)=0
        ind_d8(i,j)=0
        ind_d10(i,j)=0
        ind_C6(i,j)=0
        ind_C8(i,j)=0
        ind_C10(i,j)=0
       enddo
       do j=1,3
        sites(j,i)=0.0d0
       enddo
      enddo

      open (8,file=path//'data_ccdata')

      read(8,*) cnull
      do i=1,25
        read(8,*) (sites(j,i),j=1,3)
      enddo
      read(8,*) cnull
      read(8,*) (chrg(j),j=1,5)

      read(8,*) cnull
      read(8,*) (ind_charge(j),j=1,5)

      read(8,*) cnull
      do i=1,25
        read(8,*) (ind_beta(i,j),j=1,25)
      enddo
      read(8,*) cnull
      do i=1,5
        read(8,*) (ind_d1(i,j),j=1,5)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_d6(i,j),j=1,3)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_d8(i,j),j=1,3)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_d10(i,j),j=1,3)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_c6(i,j),j=1,3)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_c8(i,j),j=1,3)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_c10(i,j),j=1,3)
      enddo

      do i=1,6250
        ind_data(i)=ind_data1(i)
      enddo

      close(8)
      return
      end
c-------------------------------------------------------










c------------------------------------------------------------------------------
c
      subroutine Potential_ang(Etot,omegA,omegB,R_AB,valA,valB)
      implicit real*8 (a-h,o-z)
      dimension Oa(3),Ha1(3),Ha2(3),Ob(3),Hb1(3),Hb2(3)
      dimension omegA(3),omegB(3),valA(3),valB(3)

      common /ddaattaa/ param(18,6),parab(84,6,6),
     .       c(1000),cc(2000),params(1000),
     .       chrg(25),sites(3,25),ind_data(6250),
     .       nsitea,nsiteb,nparsall,iembed,ipotparts,icc,
     .       iemonomer,iembedinterf,ixz

      call interface(omegA,omegB,R_AB,valA,valB,
     .               Oa,Ha1,Ha2,Ob,Hb1,Hb2,
     .               iembedinterf,ixz)

      call CCpol_xyz
     .      (nnr,Oa,Ha1,Ha2,Ob,Hb1,Hb2,iembed,ipotparts,icc,Etot,
     .      param,parab,nsitea,nsiteb,c,
     .      cc,params,nparsall,
     .      ind_data,chrg,sites)

      if (iemonomer.eq.1) then

c        pi=dacos(-1.d0)
c        rad2deg=180.d0/pi
c        a0=0.529177249d0
        h2kcal=627.510d0

c        rA1=valA(1)*a0
c        rA2=valA(2)*a0
c        thA=valA(3)*rad2deg
c        rB1=valB(1)*a0
c        rB2=valB(2)*a0
c        thB=valB(3)*rad2deg

c        write(6,'(a,3f13.8)') "rA1, rA2, thA:  ",rA1, rA2, thA
c        write(6,'(a,3f13.8)') "rB1, rB2, thB:  ",rB1, rB2, thB

        call pots(vA,valA(1),valA(2),valA(3))
        call pots(vB,valB(1),valB(2),valB(3))
c        write(6,'(a,f13.8)') "vA= ",vA*h2kcal
c        write(6,'(a,f13.8)') "vB= ",vB*h2kcal

c        write(6,'(a,f10.5)') "Etot(dimer only): ",Etot
        Etot=Etot+(vA+vB)*h2kcal

      endif

      return
      end
c
c -----------------------------------------------------------------------------
c The Cartesian coordinates of the atoms in the complex are calculated from
c the geometry given using angles and distances.
c To perform this procedure one has to define the type of the embedding which
c is to be used.
c
      subroutine interface(omegA,omegB,R_AB,valA,valB,
     .                     OOOa,HHHa1,HHHa2,OOOb,HHHb1,HHHb2,
     .                     iembedang,ixz)
      implicit real*8 (a-h,o-z)

      dimension O1(3),H1(3),H2(3)

      dimension Rot(3,3),Rota(3,3),Rotb(3,3)
      dimension omegA(3),omegB(3),valA(3),valB(3)
      dimension Oa(3),Ha1(3),Ha2(3)
      dimension Ob(3),Hb1(3),Hb2(3)
      dimension OOOa(3),HHHa1(3),HHHa2(3)
      dimension OOOb(3),HHHb1(3),HHHb2(3)
      parameter (bohr2A=0.52917724924d0)

      A2bohr=1.0d0/bohr2A

      R_AB_ang=R_AB*bohr2A

c calculate the rotation matrix which puts the Oz axis at the position 
c indicated by the Euler angles
      phi  =omegA(1)
      theta=omegA(2)
      chi  =omegA(3)
      call eulerrot(Rota,phi,theta,chi,1)

c put the molecule A at the initial position in the xz (yz) surface (embedding)
      call put_flex(iembedang,valA,Oa,Ha1,Ha2,ixz)
c rotate the molecule
      do i=1,3
        OOOa(i)=0.0d0
        HHHa1(i)=0.0d0
        HHHa2(i)=0.0d0
        do k=1,3
         OOOa(i)=OOOa(i)+Rota(i,k)*Oa(k)
         HHHa1(i)=HHHa1(i)+Rota(i,k)*Ha1(k)
         HHHa2(i)=HHHa2(i)+Rota(i,k)*Ha2(k)
        enddo
      enddo

      phi  =omegB(1)
      theta=omegB(2)
      chi  =omegB(3)
      call eulerrot(Rotb,phi,theta,chi,1)
      call put_flex(iembedang,valB,Ob,Hb1,Hb2,ixz)
      do i=1,3
        OOOb(i)=0.0d0
        HHHb1(i)=0.0d0
        HHHb2(i)=0.0d0
        do k=1,3
         OOOb(i)=OOOb(i)+Rotb(i,k)*Ob(k)
         HHHb1(i)=HHHb1(i)+Rotb(i,k)*Hb1(k)
         HHHb2(i)=HHHb2(i)+Rotb(i,k)*Hb2(k)
        enddo
      enddo
c - for the molecule B shift the coordinates by R_AB
      OOOb(3)=OOOb(3)+R_AB_ang
      HHHb1(3)=HHHb1(3)+R_AB_ang
      HHHb2(3)=HHHb2(3)+R_AB_ang

c check the coordinates of the nuclei
c      write(6,*) "positions of the nuclei:"
c      write(6,'(3d18.10)') OOOa
c      write(6,'(3d18.10)') HHHa1
c      write(6,'(3d18.10)') HHHa2
c      write(6,'(3d18.10)') OOOb
c      write(6,'(3d18.10)') HHHb1
c      write(6,'(3d18.10)') HHHb2

      return
      end
c ---------------------------------------------------------

!     ===+=========+=========+=========+=========+=========+=========+==
      Subroutine Eulerrot(rot,alpha,beta,gamma,transpose)
!     ===+=========+=========+=========+=========+=========+=========+==
*     calculate rotation matrix given the euler angles (in radians)
*       alpha  corresponds to polar coordinate phi
*       beta   corresponds to polar coordinate theta
*       gamma  corresponds to the "spin" angle chi
*
      Implicit real*8 (a-h,o-z)
      double precision rot(3,3)
      integer transpose
c
      cosa=dcos(alpha)
      sina=dsin(alpha)
      cosb=dcos(beta)
      sinb=dsin(beta)
      cosc=dcos(gamma)
      sinc=dsin(gamma)
c
      if( transpose .eq. 0 ) then 
        rot(1,1) = cosb*cosa*cosc-sina*sinc
        rot(1,2) = cosb*sina*cosc+cosa*sinc
        rot(1,3) = -sinb*cosc
        rot(2,1) = -cosb*cosa*sinc-sina*cosc
        rot(2,2) = -cosb*sina*sinc+cosa*cosc
        rot(2,3) = sinb*sinc
        rot(3,1) = sinb*cosa
        rot(3,2) = sinb*sina
        rot(3,3) = cosb
      else
        rot(1,1) = cosb*cosa*cosc-sina*sinc
        rot(2,1) = cosb*sina*cosc+cosa*sinc
        rot(3,1) = -sinb*cosc
        rot(1,2) = -cosb*cosa*sinc-sina*cosc
        rot(2,2) = -cosb*sina*sinc+cosa*cosc
        rot(3,2) = sinb*sinc
        rot(1,3) = sinb*cosa
        rot(2,3) = sinb*sina
        rot(3,3) = cosb
      endif
      return
      End 
c ---------------------------------------------------------
c------------------------------------------------------------------
c - znieksztalcona czasteczke umiesc na plaszczyznie YZ (XZ) w taki sposob,
c   aby srodek masy znalazl sie w poczatku ukladu wspolrzednych,
c   a wektor dwusieczny skierowany byl wzdluz osi Z, 
c - zawolaj procedure eck_rad_tst lub radau_f1_tst ze wspolrzednymi kartezjanskimi
c   uzyskanymi w poprzednim punkcie
c - czasteczke obroc na plaszczyznie wokol poczatku ukladu wspolrzednych 
c   w taki sposob, aby wektor dwusieczny pokryl sie z vecIa 
!          valA(:)  : valence coordinates r1,r2,theta of monomer A

      subroutine put_flex(iembedang,val,O,H1,H2,ixz)
      implicit real*8 (a-h,o-z)
      dimension vi1(3),vi2(3)
      dimension O(3),H1(3),H2(3)
      dimension Or(3),H1r(3),H2r(3)
      dimension w1(3),w2(3),vshift(3)
      dimension Oposition(3),COM(3)

      dimension vecI(3),vecJ(3)
      dimension vz(3)
      dimension val(3),val_bohr(3)
      parameter (bohr2A=0.52917724924d0)
c ds=sin(th_ref_half), dc=cos(th_ref_half)
      a0=0.529177249d0

      pi=dacos(-1.d0)
      rad2deg=180.d0/pi
c distorted molecules put in the plane xz or yz in such a way that COM
c is in the origin and the bisection vector is oriented along the 0z axis.
      r1=val(1)*bohr2A
      r2=val(2)*bohr2A
      th=val(3)
      th_deg=th*rad2deg
      th2=th/2.0d0
      th2_deg=th_deg/2.0d0
      do i=1,3
        Or(i)=0.0d0
        H1r(i)=0.0d0
        H2r(i)=0.0d0
      enddo
      ds=dsin(th2)
      dc=dcos(th2)
      if (ixz.eq.1) then
       H1r(1)= r1*ds   
       H1r(3)=-r1*dc
       H2r(1)=-r2*ds
       H2r(3)=-r2*dc
      else
       H1r(2)=-r1*ds    ! Claude's convention
       H1r(3)=-r1*dc
       H2r(2)= r2*ds
       H2r(3)=-r2*dc
      endif

      call COMcalc(Or,H1r,H2r,COM)
      do i=1,3
        Or(i)=Or(i)-COM(i)
        H1r(i)=H1r(i)-COM(i)
        H2r(i)=H2r(i)-COM(i)
      enddo
c - call embedding procedure eck_rad_tst or radau_f1_tst 
      if (iembedang.eq.1) then
        call eck_rad_tst(Or,H1r,H2r,vecI,vecJ)
      endif
      if (iembedang.eq.2) then
        call radau_f1_tst(Or,H1r,H2r,vecI,vecJ)
      endif
c rotate the molecule around the origin to make the bisection vector
c aligned with the vector vecIa
      vz(1)=0.0d0
      vz(2)=0.0d0
      vz(3)=-1.0d0
c vecI is normalized
      ss=0.0d0
      vzn=0.0d0
      vecIn=0.0d0
      do i=1,3
        ss=ss+vz(i)*vecI(i)
        vzn=vzn+vz(i)*vz(i)
        vecIn=vecIn+vecI(i)*vecI(i)
      enddo
      vznn=dsqrt(vzn)
      vecInn=dsqrt(vecIn)
      angle=dacos(ss/(vznn*vecInn)) 

      call cross(w1,vz,vecI)
      zero=0.0d0
      if (ixz.eq.1) then
        if (w1(2).lt.zero) angle=-angle
      else
        if (w1(1).lt.zero) angle=-angle
      endif

      if (ixz.eq.1) then
c rotate the molecule in the XZ plane by the angle
       O(1)= dcos(angle)*Or(1)-dsin(angle)*Or(3)
       O(2)= 0.0d0    
       O(3)= dsin(angle)*Or(1)+dcos(angle)*Or(3)
       H1(1)= dcos(angle)*H1r(1)-dsin(angle)*H1r(3)
       H1(2)= 0.0d0
       H1(3)= dsin(angle)*H1r(1)+dcos(angle)*H1r(3)
       H2(1)= dcos(angle)*H2r(1)-dsin(angle)*H2r(3)
       H2(2)= 0.0d0
       H2(3)= dsin(angle)*H2r(1)+dcos(angle)*H2r(3)

c       O(1)= dcos(angle)*Or(1)+dsin(angle)*Or(3)
c       O(2)= 0.0d0    
c       O(3)=-dsin(angle)*Or(1)+dcos(angle)*Or(3)
c       H1(1)= dcos(angle)*H1r(1)+dsin(angle)*H1r(3)
c       H1(2)= 0.0d0
c       H1(3)=-dsin(angle)*H1r(1)+dcos(angle)*H1r(3)
c       H2(1)= dcos(angle)*H2r(1)+dsin(angle)*H2r(3)
c       H2(2)= 0.0d0
c       H2(3)=-dsin(angle)*H2r(1)+dcos(angle)*H2r(3)
      else
c rotate the molecule in the YZ plane by the angle
       O(1)= 0.0d0     ! for Claude's convention
       O(2)= dcos(angle)*Or(2)+dsin(angle)*Or(3)
       O(3)=-dsin(angle)*Or(2)+dcos(angle)*Or(3)
       H1(1)= 0.0d0
       H1(2)= dcos(angle)*H1r(2)+dsin(angle)*H1r(3)
       H1(3)=-dsin(angle)*H1r(2)+dcos(angle)*H1r(3)
       H2(1)= 0.0d0
       H2(2)= dcos(angle)*H2r(2)+dsin(angle)*H2r(3)
       H2(3)=-dsin(angle)*H2r(2)+dcos(angle)*H2r(3)
      endif
c sprawdzmy liczac vecI - powinien byc wzdluz osi OZ
c      if (iembedang.eq.1) then
c        call eck_rad_tst(O,H1,H2,vecI,vecJ)
c      endif
c      if (iembedang.eq.2) then
c        call radau_f1_tst(O,H1,H2,vecI,vecJ)
c      endif
c      write(37,*) ' check:',(vecI(i),i=1,3)
c sprawdzam geometrie czasteczki
c      rr1=0.0d0
c      rr2=0.0d0
c      sss=0.0d0
c      do j=1,3
c        rr1=rr1+(H1(j)-O(j))**2
c        rr2=rr2+(H2(j)-O(j))**2
c        sss=sss+(H1(j)-O(j))*(H2(j)-O(j))
c      enddo
c      rr1=dsqrt(rr1)
c      rr2=dsqrt(rr2)
c      ddd=sss/(rr1*rr2)
c      aaa=dacos(ddd)*rad2deg
c      write(6,'(a,3f13.8)') "rr1,rr2,angle:  ",rr1,rr2,aaa

      return
      end
c------------------------------------------------------------------







c----------------------------------------------------------
      subroutine ccpol8s_dimer(imode,Oa,Ha1,Ha2,Ob,Hb1,Hb2,Erigid,
     .                         c,params,nparsall,
     .                         ind_data,chrg,sites)
      implicit real*8 (a-h,o-z) 
      PARAMETER (nsite=25,maxdat=3000,maxlin=2000)
      PARAMETER (maxpar=1000)
      parameter (nsites=2*nsite*maxdat)
      dimension params(maxpar)
c      save params,nparsall

      dimension c(maxlin)
      dimension ElA(3),ElB(3)
      dimension RA(3),RB(3)
      dimension T(3,3),vec(3)
      dimension sites(3,nsite)     ! site coordinates for standard molecule
      dimension sitesAB(3,nsites) ! site coords for each of translated and rotated 2*ndat molecules
      dimension aj(maxlin)
      dimension rsites(3,nsite)     ! site coordinates for standard molecule - robocza
      dimension Oa(3),Ha1(3),Ha2(3)
      dimension Ob(3),Hb1(3),Hb2(3)
      dimension COMa(3),COMb(3)

      integer omp_get_num_threads,omp_get_thread_num
      character*1 cnull
      dimension ind_data(6250)
      dimension chrg_omp(25)
      dimension sites_omp(3,25)
      character(len=11) :: path = "PES_H2OH2O/"
c      save sites 

      zero=0.d0
      bohr2a=0.529177249d0
      h2kcal=627.510d0

      pi180=dacos(-1.d0)/180.d0

c        id=omp_get_thread_num()

      if (imode.eq.-1) then
      open (7,file=path//'data_CCpol8s')
      read (7,*) nparsall
      if (nparsall.gt.maxpar) stop 010
      do 10 i=1,nparsall
         read (7,*) ii,params(i)
         if (ii.ne.i) stop 020
 10   continue

      read (7,*) nlin0
      if (nlin0.gt.maxlin) stop 030
      do 20 i=1,nlin0
         read (7,*) ii,c(i)
         if (ii.ne.i) stop 040
 20   continue
      call read_cc_data(ind_data,chrg,sites)
      close(7)
      return
      endif

c      write (*,'(a)') 'dimer energies (hartree,kcal/mol)'
c      read (*,*) ndat
c      if (ndat.gt.maxdat) stop 111

        do j=1,3
          Oa(j)=Oa(j)/bohr2a
          Ha1(j)=Ha1(j)/bohr2a
          Ha2(j)=Ha2(j)/bohr2a
          Ob(j)=Ob(j)/bohr2a
          Hb1(j)=Hb1(j)/bohr2a
          Hb2(j)=Hb2(j)/bohr2a
        enddo

c If our program is used independently, use this subroutine to check whether our
c molecule is indeed the water in the reference geometry
c      call check_molecule(Oa,Ha1,Ha2)   
c      call check_molecule(Ob,Hb1,Hb2)   

        call fill_sites(Oa,Ha1,Ha2,sites,rsites)
        indA = 0
        do ns=1,nsite
          do jj=1,3
            sitesAB(jj,indA+ns) = rsites(jj,ns)
          enddo
        enddo

        call fill_sites(Ob,Hb1,Hb2,sites,rsites)
        indB = nsite
        do ns=1,nsite
          do jj=1,3
            sitesAB(jj,indB+ns) = rsites(jj,ns)
          enddo
        enddo

      itwo=2
      call indN_iter (itwo, sitesAB(1,indA+1), nsite, Eind, chrg)

      Eind = Eind           ! hartree
c      Eind = Eind*h2kcal    ! kcal/mol

      indA = 1
      indB = nsite +1
      call U0(nsite,sitesAB(1,indA),sitesAB(1,indB),dummy,a0,aj,nlin,
     .        params,nparsall,
     .        ind_data,chrg,sites)

      E=Eind
      do 210 nl=1,nlin
        E=E +c(nl)*aj(nl)
 210  continue
      E=E+a0

c      write (*,'(i5,2f20.8)') nd,E,E*h2kcal
      Erigid=E*h2kcal

      return
      end
c----------------------------------------------------------------------
      subroutine U0 (nsite,sitesA,sitesB,energy, a0,aj,nlin,
     .               params,nparsall,
     .               ind_data,chrg,sites)
      implicit real*8 (a-h,o-z)
      PARAMETER (maxpar=1000)

c      dimension sitesA(3,*),sitesB(3,*), aj(*)
      dimension sitesA(3,25),sitesB(3,25), aj(144)
      dimension params(maxpar)

      character*1 cnull
      dimension ind_charge(25)
      dimension ind_beta(25,25)
      dimension ind_d1(25,25)
      dimension ind_d6(25,25)
      dimension ind_d8(25,25)
      dimension ind_d10(25,25)
      dimension ind_C6(25,25)
      dimension ind_C8(25,25)
      dimension ind_C10(25,25)
      dimension chrg(25)
      dimension sites(3,25)
      dimension ind_data(6250)
      dimension ind_data1(6250)
      equivalence (ind_data1(1),ind_charge)
      equivalence (ind_data1(26),ind_beta)
      equivalence (ind_data1(651),ind_d1)
      equivalence (ind_data1(1276),ind_d6)
      equivalence (ind_data1(1901),ind_d8)
      equivalence (ind_data1(2526),ind_d10)
      equivalence (ind_data1(3151),ind_c6)
      equivalence (ind_data1(3776),ind_c8)
      equivalence (ind_data1(4401),ind_c10)
      integer omp_get_num_threads,omp_get_thread_num

      zero=0.d0
      h2kcal=627.510d0
      do i=1,6250
        ind_data1(i)=ind_data(i)
      enddo
      
      nlin=144

      do 5 i=1,nlin
      aj(i)=zero
 5    continue

      E_exp=zero
      E_ele=zero
      E_ind=zero

      do 10 nsA=1,nsite
      do 10 nsB=1,nsite

      call distan (sitesA(1,nsA),sitesB(1,nsB),R)

      if (ind_beta(nsA,nsB).ne.0) then     ! exponential
         beta = params(ind_beta(nsA,nsB))
         eks = dexp(-beta*R)
         indlin=ind_beta(nsA,nsB)-98
         if (indlin.lt.0) indlin=indlin+65
         ind0= indlin
         ind1= ind0+36
         ind2= ind1+36
         ind3= ind2+36
         aj(ind0) = aj(ind0) +eks
         aj(ind1) = aj(ind1) +eks*R
         aj(ind2) = aj(ind2) +eks*R*R
         aj(ind3) = aj(ind3) +eks*R*R*R
c         E_exp=...
      end if

      if (ind_charge(nsA)*ind_charge(nsB).ne.0) then     ! elst
         qA = params(ind_charge(nsA))
         qB = params(ind_charge(nsB))
         d1 = params(ind_d1(nsA,nsB))
         f1 = damp(1,d1,R)
         E_ele = E_ele+ f1* qA*qB/R
      end if

      if (ind_d6(nsA,nsB).ne.0) then     ! ind-disp
         d6 = params(ind_d6(nsA,nsB))
         d8 = params(ind_d8(nsA,nsB))
         d10 = params(ind_d10(nsA,nsB))
         C6 = params(ind_C6(nsA,nsB))
         C8 = params(ind_C8(nsA,nsB))
         C10 = params(ind_C10(nsA,nsB))
         f6 = damp(6,d6,R)
         f8 = damp(8,d8,R)
         f10 = damp(10,d10,R)
         R2 = R*R
         R6 = R2*R2*R2
         R8 = R6*R2
         R10 = R8*R2
         E_ind = E_ind -f6*C6/R6 -f8*C8/R8 -f10*C10/R10
      end if

 10   continue

c      energy = E_exp +E_ele +E_ind
      a0 = E_ele +E_ind

c!test1
c      write (8,*) 'E_exp=',E_exp*h2kcal
c      write (8,*) 'E_ele=',E_ele*h2kcal
c      write (8,*) 'E_ind=',E_ind*h2kcal
c      write (8,*) 'energy=',energy*h2kcal
c      sum=zero
c      do 200 i=1,40
c      sum=sum + aj(i)*params(i+43)
c 200  continue
c      write (8,*) 'test:',sum*h2kcal
c!test0      


      end
c---------------------------------------------------------------------
      subroutine indN_iter(N,sites,nsite,energy,chrg)
c Compute the iterative induction between N water molecules.
c The electrostatic field supplied by routine efield.
c
c Assume that the dipole center is at NOT the center
c of mass, but in special point.
      implicit real*8 (a-h,o-z)
      PARAMETER (Nmax=1024)
      PARAMETER (maxsite=30)
      parameter (maxNsite=Nmax*maxsite)
      parameter (maxit=200)

      dimension sites(3,maxNsite) ! site coords for all molecules
      dimension R(3,Nmax)         ! positions of polarizable sites
      dimension G2(3,Nmax)        ! induced dipoles
      dimension E0x(Nmax),E0y(Nmax),E0z(Nmax) ! fields from static charges
      dimension epom(3)
      dimension dist(Nmax,Nmax)   ! distances between pol.centers (**(-3))

      dimension chrg(25)
      integer omp_get_num_threads,omp_get_thread_num

c        ns=omp_get_num_threads()
c        id=omp_get_thread_num()
c        nn=30+id
c        write(nn,*) "in indN_iter  0  id=",id
c        call flush(nn)
      
      pol=9.922d0                ! exprmt average polarizability
      sig=0.367911875040999981d0 ! position of polarizable center 
                                 ! same as in 3B fits
      plen=1.1216873242d0
      dmpfct=1.d0
      zero=0.d0

      do j=1,N
       do i=1,3
        G2(i,j) = zero
       enddo
      enddo

c Compute positions of polarizable sites and the distance between them
c      dist(i,j)= 1.d0/rij**3

      do 10 i=1,N
      ind1= (i-1)*nsite +1
      ind2= (i-1)*nsite +2
      ind3= (i-1)*nsite +3
      do 11 ii=1,3
       pom = 0.5d0* (sites(ii,ind2) + sites(ii,ind3))
       R(ii,i) = sites(ii,ind1) + sig*(pom-sites(ii,ind1))/plen
 11   continue
 10   continue

      do 20 i=1,N
      do 21 j=i+1,N
      dist(i,j) = 0.d0
      do 25 ii=1,3
 25   dist(i,j) = dist(i,j) +(R(ii,i)-R(ii,j))**2.d0
      dist(i,j)=dist(i,j)**(-1.5d0)
      dist(j,i)=dist(i,j)
 21   continue
 20   continue

c
c Permanent charges fields
c
      do 30 i=1,N
      E0x(i) = zero
      E0y(i) = zero
      E0z(i) = zero
      do 35 j=1,N
      if (j.eq.i) goto 35
      index= (j-1)*nsite +1
      call efield_bohr(R(1,i),sites(1,index),nsite,epom,chrg) ! field of j on i
      E0x(i) = E0x(i) +epom(1)
      E0y(i) = E0y(i) +epom(2)
      E0z(i) = E0z(i) +epom(3)
 35   continue
 30   continue


C   LETS DO IT (the iterations, that is...)

      thr_iter = 1.d-20
      change = 10.d0
      isteps = 0

      DO WHILE (change.gt.thr_iter.and.isteps.lt.maxit)   
c     -----------------------------------------------

      energy= 0.0d0
      change= 0.0d0

      do 100 i=1,N
       E1X = E0x(i)
       E1Y = E0y(i) 
       E1Z = E0z(i) 
       do 200 j=1,N
          if(j.eq.i) go to 200
c calculate the total electric field at G1 of molecule i
         call TTTprod(R(1,i),R(1,j),G2(1,j),dist(i,j),epom)
          E1X = E1X + dmpfct*epom(1)
          E1Y = E1Y + dmpfct*epom(2)
          E1Z = E1Z + dmpfct*epom(3)
 200   continue

c total field at ii ready. Place induced dipole....
      polE1x = pol*E1X
      polE1y = pol*E1Y
      polE1z = pol*E1Z

      change=(G2(1,i)-polE1X)**2
     >      +(G2(2,i)-polE1Y)**2
     >      +(G2(3,i)-polE1Z)**2 + change
      G2(1,i)=polE1X
      G2(2,i)=polE1Y
      G2(3,i)=polE1Z
      energy =-0.5d0*pol*(E1X*E0X(i)
     >                   +E1Y*E0Y(i)
     >                   +E1Z*E0Z(i)) + energy

 100  continue

c      write (*,*) isteps,energy
      isteps = isteps +1
      END DO     !  end while
c     -----------------------

      if (isteps.ge.maxit) then
         write (*,*) 'No convergence in indN_iter'
         write (*,'(a,g12.3)') 'energy change=',change
         write (*,'(a,g12.3)') 'thr_iter=',thr_iter
         stop
      end if

c      energy = 627.510d0*energy   ! hartree--->kcal/mol
      end
c----------------------------------------------------------
c     
c Calculate the contribution of molecule b to the field on a
c Use the CC-pol site charges.
c Units: input positions of bohr, 
c output fields in au
c
      subroutine efield_bohr(veci,sitebt,nsiteb,e,chrg)
      implicit real*8 (a-h,o-z)
      PARAMETER (nsitemax=30)

      dimension veci(3),sitebt(3,1),e(3),chrg(25)
      dimension sep(3,nsitemax),sepl(nsitemax)

      integer omp_get_num_threads,omp_get_thread_num
      dimension sites(3,25)
      character*1 cnull
c       save chrg
c      data a0 /0.529177249d0/
      a0=1.d0
      crgfct=18.22262373d0

c
c Position of polarizable center of A is in veci(*)
      do isite=1,nsiteb
       if(chrg(isite).ne.0.d0) then   ! if charged site
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

      do isite=1,nsiteb
       if(chrg(isite).ne.0.d0) then   ! if charged site
        do k=1,3
         e(k) = e(k) + a0*a0*chrg(isite)*sep(k,isite)*sepl(isite)
        end do
       endif      ! if charged site
      end do
c
      return
      end
c----------------------------------------------------------
      subroutine Av (A,v,b)
c calculates A*v=b
      implicit real*8 (a-h,o-z)
      dimension A(3,3),v(3),b(3)
      zero=0.d0

      do 10 i=1,3
      sum=zero
      do 20 j=1,3
      sum=sum +a(i,j)*v(j)
 20   continue
      b(i)=sum
 10   continue

      end
c----------------------------------------------------------
      subroutine distan (R1,R2,d)
      implicit real*8 (a-h,o-z)
      dimension R1(3),R2(3)

      d=0.d0
      do 10 i=1,3
      R12i=R1(i)-R2(i)
      d=d +R12i*R12i
 10   continue
      d=sqrt(d)
      end
c----------------------------------------------------------
      function damp (n,beta,r)
c
c     calculate the damping factor (small R correct)
c
      implicit real*8 (a-h,o-z)
      br=beta*r
c The following line added by RB, Sept. 18 1997
      if(br.eq.0.d0) then
       damp = 0.d0
       return
      endif
      sum=1.0d0
      term=1.0d0
      ncn=n
      do i=1,ncn
        term=term*br/i
        sum=sum+term
      enddo
      damp=1.0d0 - dexp(-br)*sum
c     in case of d --> 0 use
c     d=1.0d0 - dexp(-br)*sum = sum_m=ncn+1^\infty br^m/m!
      if(dabs(damp).lt.1.0d-8) then
        damp=0.0d0
        do i=ncn+1,1000
          term=term*br/i
          damp=damp+term
          if(term/damp .lt. 1.0d-8) go to 111
        enddo
        write(6,*) 'No convergence in damp'
  111 continue
      damp=damp*dexp(-br)
      endif
c     write(6,'(i4,2f10.5,e20.10)') n,beta,r,d
      return
      end
c----------------------------------------------------------
      subroutine fill_sites(O,H1,H2,sites,rsites)
      implicit real*8 (a-h,o-z)
      PARAMETER (nsite=25,maxdat=3000,maxlin=2000)
      dimension sites(3,nsite)     ! site coordinates for standard molecule
      dimension rsites(3,nsite)  
      dimension O(3),H1(3),H2(3)
      dimension COM(3)
      dimension ex(3),ey(3),ez(3)
      dimension v1(3),v2(3),v3(3)

c      save sites 

      bohr2a=0.529177249d0

c ss1 - bond length
c      dv1pv2=2.0d0*ss1*dcos(theta/2.0d0)
c      dv1mv2=2.0d0*ss1*dsin(theta/2.0d0)
      dv1pv2=1.99230765895d0
      dv1mv2=2.907303924565d0

      pi180=dacos(-1.d0)/180.d0

c First prepare two vectors, which span the plane defined by the molecule.
c Let's call them ez and ex, and require that they are oriented with respect to 
c the molecule in a standard way.
c ez is a bisector vector  

      call COMcalc(O,H1,H2,COM)

      do j=1,3
        v1(j)=H1(j)-COM(j)
        v2(j)=H2(j)-COM(j)
      enddo
      do j=1,3
        ez(j)=-(v1(j)+v2(j))
        ex(j)=v2(j)-v1(j)
      enddo
      do j=1,3
        ez(j)=ez(j)/dv1pv2
        ex(j)=ex(j)/dv1mv2
      enddo
      call cross(ey,ez,ex)

      do kk=1,nsite
        rsites(1,kk)=
     .       ex(1)*sites(1,kk) + ey(1)*sites(2,kk) + ez(1)*sites(3,kk)
        rsites(2,kk)=
     .       ex(2)*sites(1,kk) + ey(2)*sites(2,kk) + ez(2)*sites(3,kk)
        rsites(3,kk)=
     .       ex(3)*sites(1,kk) + ey(3)*sites(2,kk) + ez(3)*sites(3,kk)
c        write(6,666) (rsites(j,kk),j=1,3)," rsites(*,kk)"
      enddo

      do kk=1,nsite
        do j=1,3
          rsites(j,kk)=rsites(j,kk)+COM(j)
        enddo
      enddo

  666 format(3f18.12,a) 
      return
      end
c----------------------------------------------------------
      subroutine cross (z,x,y)
c calculates z = vector cross poduct of x an y
      implicit real*8 (a-h,o-z)
      dimension x(3),y(3),z(3)

      z(1)=x(2)*y(3)-x(3)*y(2)
      z(2)=x(3)*y(1)-x(1)*y(3)
      z(3)=x(1)*y(2)-x(2)*y(1)
      end
c------------------------------------------------------------------
      subroutine COMcalc (O1,H1,H2,COM)
      implicit real*8 (a-h,o-z)
      dimension O1(3),H1(3),H2(3),COM(3)
      real*8 mO,mH,M
      mO=15.994 914 6221d0
      mH=1.007 825 032 1d0

      M=mO+mH+mH
      do 10 i=1,3
      COM(i)= (mO*O1(i) +mH*H1(i) +mH*H2(i))/M
 10   continue

      end
c------------------------------------------------------------------







       subroutine driver_potss_sapt5sf(carta,cartb,val,ipotparts,
     .                            param,parab,nsitea,nsiteb,c)
       implicit real*8 (a-h,o-z)
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
       dimension carta(3,3),cartb(3,3)
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension c(1000)

      a0=0.529177249d0

       icart=1  ! choose this flag if the cartesian coordinates are provided
                ! in the following format 
                ! (for (O1,H1,H2) interacting with (O2,H3,H4)):
                ! xO1  yO1  zO1
                ! xH1  yH1  zH1
                ! xH2  yH2  zH2
                ! xO2  yO2  zO2
                ! xH3  yH3  zH3
                ! xH4  yH4  zH4
                ! Angstrom is a default unit

        pi = dacos(-1.d0)
        d2rad = pi/180.d0
        rad2d = 1.d0/d2rad

        r0=0.9716257d0
        theta0=104.69d0

      imode = 0
      alphaAv=0.0d0

        do ii=1,3
          do jj=1,3
            carta(ii,jj)=carta(ii,jj)/a0
            cartb(ii,jj)=cartb(ii,jj)/a0
          enddo
        enddo

        call poten(r1,r2,theta1,r3,r4,theta2,
     .      rin,alphaAv,betaAv,gammaAv,alphaBv,betaBv,gammaBv,val,imode,
     .      icart,carta,cartb,ipotparts,
     .        param,parab,nsitea,nsiteb,c)

201   format(2(1x,3(1x,f9.4)),7(1x,f9.4),f15.8,f5.1)
202   format(1x,3f15.8/1x,3f15.8/1x,3f15.8/
     .       1x,3f15.8/1x,3f15.8/1x,3f15.8,f15.8)

      return
      end
c========================================================================
      subroutine poten(r1_ang,r2_ang,theta1_deg,
     .                 r3_ang,r4_ang,theta2_deg,
     .        R,alphaAv,betaAv,gammaAv,alphaBv,betaBv,gammaBv,val,mode,
     .        icart,carta,cartb,ipotparts,
     .        param,parab,nsitea,nsiteb,c)
      implicit real*8 (a-h,o-z)
c
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
      dimension oa(3), ob(3), values(100)  !     ???
c
      dimension siteat(3,nsitemax,-1:naamax),
     .          sitebt(3,nsitemax,-1:nbbmax)
      dimension sa(3,-1:naamax),sb(3,-1:nbbmax)
      dimension itypea(nsitemax),itypeb(nsitemax)
      dimension itypus(ntypemax,ntypemax,2)
      dimension numtm(2)   ! to be set in potparts....
      dimension sitex(3,nsitemax),sx(3)
      dimension carta(3,3),cartb(3,3)
c      save pi,d2rad,rad2d,iaa,ibb

      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension c(1000)
      integer omp_get_num_threads,omp_get_thread_num
      character*50 nameofsaptdatafile
c      save param,parab,nsitea,nsiteb,c

c        ns=omp_get_num_threads()
c        id=omp_get_thread_num()
c
c If mode = -1 -- read input
c
      if(mode.eq.-1) then
        call data1(param,parab,nsitea,nsiteb,c)
       return
      endif

        pi = dacos(-1.d0)
        d2rad = pi/180.d0
        rad2d = 1.d0/d2rad
        iaa=-1
        ibb=-1

       call set_sites(carta,sitex,sx,itypea,nsitea)

       do j=1,3
        sa(j,iaa)=sx(j)
       end do
       do ia=1,nsitea
        do j=1,3
         siteat(j,ia,iaa)=sitex(j,ia)
        end do
       end do

       call set_sites(cartb,sitex,sx,itypeb,nsiteb)

       do j=1,3
        sb(j,ibb)=sx(j)
       end do
       do ib=1,nsiteb
        do j=1,3
         sitebt(j,ib,ibb)=sitex(j,ib)
        end do
       end do

c
c Zero out the used-type matrix
c
      do k1=1,ntypemax
       do k2=1,ntypemax
        itypus(k1,k2,1) = 0
        itypus(k1,k2,2) = 0
       end do
      end do
c
c Loop over sites in A and B
c
      iii = 1
      val = 0.d0
c
      do ia=1,nsitea
       do ib=1,nsiteb
        valp = 0.d0
c
c Compute the distance between the two sites
c
        rij = dist(siteat(1,ia,iaa),sitebt(1,ib,ibb))
c        write(44,*) "r(",ia,",",ib,")=",rij
c
c Compute the "basis functions" values pertaining to the
c required form of the potential. numtm(1) is the number of terms 
c returned for the "symmetric" part, numtm(2) - number
c of antisymmetric terms, numt = total number of terms, including constant part 
c
       if (ipotparts.eq.1) then
        call potparts(iaa,ibb,rij,ia,ib,numt,numtm,values,
     .                c,sa,sb,param,parab,itypea,itypeb)
        elseif(ipotparts.eq.0) then
        call potparts_old(iaa,ibb,rij,ia,ib,numt,numtm,values,
     .                c,sa,sb,param,parab,itypea,itypeb)
       else
        write(6,*) "wrong value of ipotparts:",ipotparts
        stop
       endif
c
c If abs(ntpot)>10, the last values is the fixed part.
c
        ntpot=124161
        numt0 = numt
        if(abs(ntpot).gt.10) then
         numt0 = numt-1
         valp = valp + values(numt)
        endif
c
c decide which basis function it is (basis function is determined
c by the pair of site types and the position in values vector).
c iii is the first basis function for this pair of types. It is
c updated into the next available position for the next pair of types.
c
c RB !!!!!!! change the following if numtm matrix added to potparts !!!
c
c GM+RB: additional loop over permutational symmetry
c
        itsmax = 1
        if(itypea(ia).ne.itypeb(ib)) itsmax = 3
        do its=1,itsmax,2
c
         itsm = min(its,2)
         itu = itypus(itypea(ia),itypeb(ib),itsm)
         if(itu.eq.0) then
          itypus(itypea(ia),itypeb(ib),itsm) = iii
          itypus(itypeb(ib),itypea(ia),itsm) = iii
          itu = iii
          iii = iii + numtm(itsm)
         endif
c
         do i=1,numtm(itsm)
          itu1 = itu+i-1
          if(itsm.eq.1) then
           valp = valp + c(itu1)*values(i)
          else
           valp = valp + c(itu1)*values(i+numtm(1))
          endif
         end do
c
        end do   ! loop over its
        val = val + valp

       end do
      end do
c
c Don't forget to add asymptotics when the fitting part tested...
c
      call dipind(iaa,ibb,R,fcind,
     .     sa,sb,param,parab,siteat,sitebt,nsitea,nsiteb,itypea,itypeb)
      val = val + fcind

c         write(94,200)
c     .    iaa,r1_ang,r2_ang,theta1_deg,(sa(j,iaa),j=1,3),
c     .    ibb,r3_ang,r4_ang,theta2_deg,(sb(j,ibb),j=1,3),
c     .    id,R,(oa(j)*rad2d,j=2,3),(ob(j)*rad2d,j=1,3),val,enfitpj
200   format(2(1(1x,i4),6(1x,f9.4)),1(1x,i6),6(1x,f9.4),f15.8,f5.1)

      return
      end
c============================================================================
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
c============================================================================
c
      subroutine potparts(iaa,ibb,rij,ia,ib,numt,numtm,values,
     .                    c,sa,sb,param,parab,itypea,itypeb)
      implicit real*8 (a-h,o-z)
c
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
      dimension values(100),numtm(2)
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension itypea(nsitemax),itypeb(nsitemax)
      dimension c(1000)
      dimension sa(3,-1:naamax),sb(3,-1:nbbmax)
       
      a0=0.529177249d0
      au2cal=627.510d0
      
c
c ***********************************************************
c GM-modified: more linear basis functions inc squares (but now with couplings
c between monomers--differs from 24161 in that accidently neglected a3 
c coplings put in)
c RB-modified
c ntpot = 124161: exponential*(1+sum r^n, n=1...3) + r^-(6-8-10) + elst, totally
c non-linear      FLEXIBILIZED
c ***********************************************************
      ntpot=124161
      if(ntpot.eq.124161) then
c
c Previously fitted rigid contributions
c hen dealing with type=6, assume the parameters of H.
c
      ita = itypea(ia)
      itb = itypeb(ib)
c      if(ita.eq.6) ita = 2
c      if(itb.eq.6) itb = 2
      beta = parab(1,ita,itb)
      alpha = parab(2,ita,itb)
      a = dexp(alpha)
c Rigid C6,C8,C10 understood as symmetric in s_A, s_B
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
      dmp6 = parab(7,itypea(ia),itypeb(ib))
      dmp8 = parab(8,itypea(ia),itypeb(ib)) 
      dmp10= parab(9,itypea(ia),itypeb(ib))
c parameter index numbers updated
      a1 = parab(38,itypea(ia),itypeb(ib))
      a2 = parab(39,itypea(ia),itypeb(ib))
      a3 = parab(40,itypea(ia),itypeb(ib))
c ccccc
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
c and previously fitted flexible contributions
c
      s1=sa(1,iaa)
      s2=sa(2,iaa)
      s3=sa(3,iaa)
      s4=sb(1,ibb)
      s5=sb(2,ibb)
      s6=sb(3,ibb)
      signa=1.0d0
      signb=1.0d0
      if (ia .eq. 3) signa=-1.0d0
      if (ib .eq. 3) signb=-1.0d0
      s3 = signa*s3
      s6 = signb*s6
      qa
     & =qa
     & +param(2,itypea(ia))*s1
     & +param(3,itypea(ia))*s2
     & +param(4,itypea(ia))*s3
     & +param(5,itypea(ia))*s1*s2
     & +param(6,itypea(ia))*s2*s3
     & +param(7,itypea(ia))*s1*s1
     & +param(8,itypea(ia))*s2*s2
     & +param(9,itypea(ia))*s3*s3
      qb
     & =qb
     & +param(2,itypeb(ib))*s4
     & +param(3,itypeb(ib))*s5
     & +param(4,itypeb(ib))*s6
     & +param(5,itypeb(ib))*s4*s5
     & +param(6,itypeb(ib))*s5*s6
     & +param(7,itypeb(ib))*s4*s4
     & +param(8,itypeb(ib))*s5*s5
     & +param(9,itypeb(ib))*s6*s6
c Make basis functions even in s3, s6 if site is NOT H but not for 
c charges--thus charges have been moved above this--or in dipind 
c since charges and dipole in dipind were not fited with this definition
c      if(itypea(ia).ne.2) s3 = dabs(s3)
c      if(itypeb(ib).ne.2) s6 = dabs(s6)
      if(itypea(ia).ne.2) s3 = s3*s3
      if(itypeb(ib).ne.2) s6 = s6*s6
c Add the flexible part of beta here....
      if(itypea(ia).eq.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*(s3+s6)
       beta = beta + parab(46,itypea(ia),itypeb(ib))*(s3*s3+s6*s6)
      elseif(itypea(ia).lt.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*s3
       beta = beta + parab(42,itypea(ia),itypeb(ib))*s6
       beta = beta + parab(46,itypea(ia),itypeb(ib))*s3*s3
       beta = beta + parab(47,itypea(ia),itypeb(ib))*s6*s6
      elseif(itypea(ia).gt.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*s6
       beta = beta + parab(42,itypea(ia),itypeb(ib))*s3
       beta = beta + parab(47,itypea(ia),itypeb(ib))*s3*s3
       beta = beta + parab(46,itypea(ia),itypeb(ib))*s6*s6
      endif
      beta = dabs(beta)   ! make it positive
c end of flex beta part....
c Add the flexible part of alpha here....
      if(itypea(ia).eq.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*(s3+s6)
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*(s3*s3+s6*s6)
      elseif(itypea(ia).lt.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*s3
       alpha = alpha + parab(44,itypea(ia),itypeb(ib))*s6
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*s3*s3
       alpha = alpha + parab(49,itypea(ia),itypeb(ib))*s6*s6
      elseif(itypea(ia).gt.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*s6
       alpha = alpha + parab(44,itypea(ia),itypeb(ib))*s3
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*s6*s6
       alpha = alpha + parab(49,itypea(ia),itypeb(ib))*s3*s3
      endif
      a = dexp(alpha)   ! make it positive
c end of flex alpha part....
      d1 = d(1,dmp1,rij)
      d6 = d(6,dmp6,rij)
      d8 = d(8,dmp8,rij)
      d10 = d(10,dmp10,rij)
      c6
     & =c6
     & +parab(11,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(14,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(17,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(20,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(23,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(26,itypea(ia),itypeb(ib))*(s2*s5)
      c8
     & =c8
     & +parab(12,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(15,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(18,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(21,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(24,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(27,itypea(ia),itypeb(ib))*(s2*s5)
      c10
     & =c10
     & +parab(13,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(16,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(19,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(22,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(25,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(28,itypea(ia),itypeb(ib))*(s2*s5)
c Update c6, c8, c10 with asymptotic terms asymmetric in s_A, s_B:
      c6as = 0.d0
      c8as = 0.d0
      c10as = 0.d0
      if(itypea(ia).ne.itypeb(ib)) then
      c6as
     & =c6as
     & +parab(29,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(32,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(35,itypea(ia),itypeb(ib))*(s2-s5)
      c8as
     & =c8as
     & +parab(30,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(33,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(36,itypea(ia),itypeb(ib))*(s2-s5)
      c10as
     & =c10as
     & +parab(31,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(34,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(37,itypea(ia),itypeb(ib))*(s2-s5)
       if(itypea(ia).gt.itypeb(ib)) then
        c6as = -c6as
        c8as = -c8as
        c10as = -c10as
       endif
      endif
      c6 = c6 + c6as
      c8 = c8 + c8as
      c10 = c10 + c10as
c
      if (beta.gt.0.d0) then 
      if(itypea(ia).ne.2.and.itypeb(ib).ne.2) then ! no H involved
       if(itypea(ia).eq.itypeb(ib)) then
        numtm(1) = 40
        numtm(2) =  0 !28
       else   ! test: make it the same as symmetric
        numtm(1) = 40
        numtm(2) = 28
       endif
      else
       if(itypea(ia).eq.itypeb(ib)) then
        numtm(1) = 40
        numtm(2) =  0 !28  ! no s-expansion if none of sites is H
       else
        numtm(1) = 40
        numtm(2) = 28
       endif
      endif
      numt = numtm(1) + numtm(2) + 1
c
      val0=a*dexp(-beta*rij)
      val1=val0*rij
      val2=val1*rij
      val3=val2*rij
       values(numt) = val0 + a1*val1 + a2*val2 + a3*val3
     & + d1*qa*qb/rij
     & - d6*c6/rij**6 - d8*c8/rij**8 - d10*c10/rij**10
c
c Construct planning matrix by filling values vector
c s-Symmetric terms
      if(itypea(ia).ne.2.and.itypeb(ib).ne.2) then ! no H involved
c delta a0,1,2,3 contributions to planning matrix
      values(1)  = (s1+s4)*val0
      values(2)  = (s1+s4)*val1
      values(3)  = (s1+s4)*val2
      values(4)  = (s1+s4)*val3

      values(5)  = (s2+s5)*val0
      values(6)  = (s2+s5)*val1
      values(7)  = (s2+s5)*val2
      values(8)  = (s2+s5)*val3

      values(9)  = (s3+s6)*val0
      values(10) = (s3+s6)*val1
      values(11) = (s3+s6)*val2
      values(12) = (s3+s6)*val3


      values(13) = (s1*s2+s4*s5)*val0
      values(14) = (s1*s2+s4*s5)*val1
      values(15) = (s1*s2+s4*s5)*val2     !x
      values(16) = (s1*s2+s4*s5)*val3     !x

      values(17) = (s2*s3+s5*s6)*val0
      values(18) = (s2*s3+s5*s6)*val1
      values(19) = (s2*s3+s5*s6)*val2
      values(20) = (s2*s3+s5*s6)*val3

      values(21) = (s1*s1+s4*s4)*val0
      values(22) = (s1*s1+s4*s4)*val1
      values(23) = (s1*s1+s4*s4)*val2
      values(24) = (s1*s1+s4*s4)*val3

      values(25) = (s2*s2+s5*s5)*val0
      values(26) = (s2*s2+s5*s5)*val1
      values(27) = (s2*s2+s5*s5)*val2
      values(28) = (s2*s2+s5*s5)*val3


      values(29) = (s1*s4)*val0
      values(30) = (s1*s4)*val1
      values(31) = (s1*s4)*val2
      values(32) = (s1*s4)*val3

      values(33) = (s2*s5)*val0
      values(34) = (s2*s5)*val1
      values(35) = (s2*s5)*val2
      values(36) = (s2*s5)*val3

      values(37) = (s3*s6)*val0
      values(38) = (s3*s6)*val1
      values(39) = (s3*s6)*val2
      values(40) = (s3*s6)*val3
c For uneqal types, add s-antisymmetric terms
      if(itypea(ia).lt.itypeb(ib)) then
      values(41) = (s1-s4)*val0
      values(42) = (s1-s4)*val1
      values(43) = (s1-s4)*val2
      values(44) = (s1-s4)*val3

      values(45) = (s2-s5)*val0
      values(46) = (s2-s5)*val1
      values(47) = (s2-s5)*val2
      values(48) = (s2-s5)*val3

      values(49) = (s3-s6)*val0
      values(50) = (s3-s6)*val1
      values(51) = (s3-s6)*val2
      values(52) = (s3-s6)*val3


      values(53) = (s1*s2-s4*s5)*val0
      values(54) = (s1*s2-s4*s5)*val1
      values(55) = (s1*s2-s4*s5)*val2     !x
      values(56) = (s1*s2-s4*s5)*val3     !x

      values(57) = (s2*s3-s5*s6)*val0
      values(58) = (s2*s3-s5*s6)*val1
      values(59) = (s2*s3-s5*s6)*val2
      values(60) = (s2*s3-s5*s6)*val3

      values(61) = (s1*s1-s4*s4)*val0
      values(62) = (s1*s1-s4*s4)*val1
      values(63) = (s1*s1-s4*s4)*val2
      values(64) = (s1*s1-s4*s4)*val3

      values(65) = (s2*s2-s5*s5)*val0
      values(66) = (s2*s2-s5*s5)*val1
      values(67) = (s2*s2-s5*s5)*val2
      values(68) = (s2*s2-s5*s5)*val3
      endif
      if(itypea(ia).gt.itypeb(ib)) then
      values(41) =-(s1-s4)*val0
      values(42) =-(s1-s4)*val1
      values(43) =-(s1-s4)*val2
      values(44) =-(s1-s4)*val3

      values(45) =-(s2-s5)*val0
      values(46) =-(s2-s5)*val1
      values(47) =-(s2-s5)*val2
      values(48) =-(s2-s5)*val3

      values(49) =-(s3-s6)*val0
      values(50) =-(s3-s6)*val1
      values(51) =-(s3-s6)*val2
      values(52) =-(s3-s6)*val3


      values(53) =-(s1*s2-s4*s5)*val0
      values(54) =-(s1*s2-s4*s5)*val1
      values(55) =-(s1*s2-s4*s5)*val2     !x
      values(56) =-(s1*s2-s4*s5)*val3     !x

      values(57) =-(s2*s3-s5*s6)*val0
      values(58) =-(s2*s3-s5*s6)*val1
      values(59) =-(s2*s3-s5*s6)*val2
      values(60) =-(s2*s3-s5*s6)*val3

      values(61) =-(s1*s1-s4*s4)*val0
      values(62) =-(s1*s1-s4*s4)*val1
      values(63) =-(s1*s1-s4*s4)*val2
      values(64) =-(s1*s1-s4*s4)*val3

      values(65) =-(s2*s2-s5*s5)*val0
      values(66) =-(s2*s2-s5*s5)*val1
      values(67) =-(s2*s2-s5*s5)*val2
      values(68) =-(s2*s2-s5*s5)*val3
      endif
      else    ! at least 1 H involved
c delta a0,1,2,3 contributions to planning matrix
      values(1)  = (s1+s4)*val0
      values(2)  = (s1+s4)*val1
      values(3)  = (s1+s4)*val2
      values(4)  = (s1+s4)*val3

      values(5)  = (s2+s5)*val0
      values(6)  = (s2+s5)*val1
      values(7)  = (s2+s5)*val2
      values(8)  = (s2+s5)*val3

      values(9)  = (s3+s6)*val0
      values(10) = (s3+s6)*val1
      values(11) = (s3+s6)*val2
      values(12) = (s3+s6)*val3


      values(13) = (s1*s2+s4*s5)*val0
      values(14) = (s1*s2+s4*s5)*val1
      values(15) = (s1*s2+s4*s5)*val2     !x
      values(16) = (s1*s2+s4*s5)*val3     !x

      values(17) = (s2*s3+s5*s6)*val0
      values(18) = (s2*s3+s5*s6)*val1
      values(19) = (s2*s3+s5*s6)*val2
      values(20) = (s2*s3+s5*s6)*val3

      values(21) = (s1*s1+s4*s4)*val0
      values(22) = (s1*s1+s4*s4)*val1
      values(23) = (s1*s1+s4*s4)*val2
      values(24) = (s1*s1+s4*s4)*val3

      values(25) = (s2*s2+s5*s5)*val0
      values(26) = (s2*s2+s5*s5)*val1
      values(27) = (s2*s2+s5*s5)*val2
      values(28) = (s2*s2+s5*s5)*val3


      values(29) = (s1*s4)*val0
      values(30) = (s1*s4)*val1
      values(31) = (s1*s4)*val2
      values(32) = (s1*s4)*val3

      values(33) = (s2*s5)*val0
      values(34) = (s2*s5)*val1
      values(35) = (s2*s5)*val2
      values(36) = (s2*s5)*val3

      values(37) = (s3*s6)*val0
      values(38) = (s3*s6)*val1
      values(39) = (s3*s6)*val2
      values(40) = (s3*s6)*val3
c For uneqal types, add s-antisymmetric terms
      if(itypea(ia).lt.itypeb(ib)) then
      values(41) = (s1-s4)*val0
      values(42) = (s1-s4)*val1
      values(43) = (s1-s4)*val2
      values(44) = (s1-s4)*val3

      values(45) = (s2-s5)*val0
      values(46) = (s2-s5)*val1
      values(47) = (s2-s5)*val2
      values(48) = (s2-s5)*val3

      values(49) = (s3-s6)*val0
      values(50) = (s3-s6)*val1
      values(51) = (s3-s6)*val2
      values(52) = (s3-s6)*val3


      values(53) = (s1*s2-s4*s5)*val0
      values(54) = (s1*s2-s4*s5)*val1
      values(55) = (s1*s2-s4*s5)*val2     !x
      values(56) = (s1*s2-s4*s5)*val3     !x

      values(57) = (s2*s3-s5*s6)*val0
      values(58) = (s2*s3-s5*s6)*val1
      values(59) = (s2*s3-s5*s6)*val2
      values(60) = (s2*s3-s5*s6)*val3

      values(61) = (s1*s1-s4*s4)*val0
      values(62) = (s1*s1-s4*s4)*val1
      values(63) = (s1*s1-s4*s4)*val2
      values(64) = (s1*s1-s4*s4)*val3

      values(65) = (s2*s2-s5*s5)*val0
      values(66) = (s2*s2-s5*s5)*val1
      values(67) = (s2*s2-s5*s5)*val2
      values(68) = (s2*s2-s5*s5)*val3
      endif
      if(itypea(ia).gt.itypeb(ib)) then
      values(41) =-(s1-s4)*val0
      values(42) =-(s1-s4)*val1
      values(43) =-(s1-s4)*val2
      values(44) =-(s1-s4)*val3

      values(45) =-(s2-s5)*val0
      values(46) =-(s2-s5)*val1
      values(47) =-(s2-s5)*val2
      values(48) =-(s2-s5)*val3

      values(49) =-(s3-s6)*val0
      values(50) =-(s3-s6)*val1
      values(51) =-(s3-s6)*val2
      values(52) =-(s3-s6)*val3


      values(53) =-(s1*s2-s4*s5)*val0
      values(54) =-(s1*s2-s4*s5)*val1
      values(55) =-(s1*s2-s4*s5)*val2     !x
      values(56) =-(s1*s2-s4*s5)*val3     !x

      values(57) =-(s2*s3-s5*s6)*val0
      values(58) =-(s2*s3-s5*s6)*val1
      values(59) =-(s2*s3-s5*s6)*val2
      values(60) =-(s2*s3-s5*s6)*val3

      values(61) =-(s1*s1-s4*s4)*val0
      values(62) =-(s1*s1-s4*s4)*val1
      values(63) =-(s1*s1-s4*s4)*val2
      values(64) =-(s1*s1-s4*s4)*val3

      values(65) =-(s2*s2-s5*s5)*val0
      values(66) =-(s2*s2-s5*s5)*val1
      values(67) =-(s2*s2-s5*s5)*val2
      values(68) =-(s2*s2-s5*s5)*val3
      endif
c
      endif   ! if H involved
 
      else
c
      numt=1
      numtm(1) = 0
      numtm(2) = 0
      values(numt) =
     & + d1*qa*qb/rij
     & - d6*c6/rij**6 - d8*c8/rij**8 - d10*c10/rij**10
c
      end if 
c
      return
      endif
      end
c
c============================================================================
c
      subroutine potparts_old(iaa,ibb,rij,ia,ib,numt,numtm,values,
     .                        c,sa,sb,param,parab,itypea,itypeb)
      implicit real*8 (a-h,o-z)
c
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
      dimension values(100),numtm(2)
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension itypea(nsitemax),itypeb(nsitemax)
      dimension c(1000)
      dimension sa(3,-1:naamax),sb(3,-1:nbbmax)
       
      a0=0.529177249d0
      au2cal=627.510d0
      
c
c ***********************************************************
c GM-modified: more linear basis functions inc squares (but now with couplings
c between monomers--differs from 24161 in that accidently neglected a3 
c coplings put in)
c RB-modified
c ntpot = 124161: exponential*(1+sum r^n, n=1...3) + r^-(6-8-10) + elst, totally
c non-linear      FLEXIBILIZED
c ***********************************************************
      ntpot=124161
      if(ntpot.eq.124161) then
c
c Previously fitted rigid contributions
c hen dealing with type=6, assume the parameters of H.
c
      ita = itypea(ia)
      itb = itypeb(ib)
c      if(ita.eq.6) ita = 2
c      if(itb.eq.6) itb = 2
      beta = parab(1,ita,itb)
      alpha = parab(2,ita,itb)
      a = dexp(alpha)
c Rigid C6,C8,C10 understood as symmetric in s_A, s_B
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
      dmp6 = parab(7,itypea(ia),itypeb(ib))
      dmp8 = parab(8,itypea(ia),itypeb(ib)) 
      dmp10= parab(9,itypea(ia),itypeb(ib))
c parameter index numbers updated
      a1 = parab(38,itypea(ia),itypeb(ib))
      a2 = parab(39,itypea(ia),itypeb(ib))
      a3 = parab(40,itypea(ia),itypeb(ib))
c ccccc
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
c and previously fitted flexible contributions
c
      s1=sa(1,iaa)
      s2=sa(2,iaa)
      s3=sa(3,iaa)
      s4=sb(1,ibb)
      s5=sb(2,ibb)
      s6=sb(3,ibb)
      signa=1.0d0
      signb=1.0d0
      if (ia .eq. 3) signa=-1.0d0
      if (ib .eq. 3) signb=-1.0d0
      s3 = signa*s3
      s6 = signb*s6
      qa
     & =qa
     & +param(2,itypea(ia))*s1
     & +param(3,itypea(ia))*s2
     & +param(4,itypea(ia))*s3
     & +param(5,itypea(ia))*s1*s2
     & +param(6,itypea(ia))*s2*s3
     & +param(7,itypea(ia))*s1*s1
     & +param(8,itypea(ia))*s2*s2
     & +param(9,itypea(ia))*s3*s3
      qb
     & =qb
     & +param(2,itypeb(ib))*s4
     & +param(3,itypeb(ib))*s5
     & +param(4,itypeb(ib))*s6
     & +param(5,itypeb(ib))*s4*s5
     & +param(6,itypeb(ib))*s5*s6
     & +param(7,itypeb(ib))*s4*s4
     & +param(8,itypeb(ib))*s5*s5
     & +param(9,itypeb(ib))*s6*s6
c Make basis functions even in s3, s6 if site is NOT H but not for 
c charges--thus charges have been moved above this--or in dipind 
c since charges and dipole in dipind were not fited with this definition
c      if(itypea(ia).ne.2) s3 = dabs(s3)
c      if(itypeb(ib).ne.2) s6 = dabs(s6)
      if(itypea(ia).ne.2) s3 = s3*s3
      if(itypeb(ib).ne.2) s6 = s6*s6
c Add the flexible part of beta here....
      if(itypea(ia).eq.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*(s3+s6)
       beta = beta + parab(46,itypea(ia),itypeb(ib))*(s3*s3+s6*s6)
      elseif(itypea(ia).lt.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*s3
       beta = beta + parab(42,itypea(ia),itypeb(ib))*s6
       beta = beta + parab(46,itypea(ia),itypeb(ib))*s3*s3
       beta = beta + parab(47,itypea(ia),itypeb(ib))*s6*s6
      elseif(itypea(ia).gt.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*s6
       beta = beta + parab(42,itypea(ia),itypeb(ib))*s3
       beta = beta + parab(47,itypea(ia),itypeb(ib))*s3*s3
       beta = beta + parab(46,itypea(ia),itypeb(ib))*s6*s6
      endif
      beta = dabs(beta)   ! make it positive
c end of flex beta part....
c Add the flexible part of alpha here....
      if(itypea(ia).eq.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*(s3+s6)
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*(s3*s3+s6*s6)
      elseif(itypea(ia).lt.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*s3
       alpha = alpha + parab(44,itypea(ia),itypeb(ib))*s6
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*s3*s3
       alpha = alpha + parab(49,itypea(ia),itypeb(ib))*s6*s6
      elseif(itypea(ia).gt.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*s6
       alpha = alpha + parab(44,itypea(ia),itypeb(ib))*s3
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*s6*s6
       alpha = alpha + parab(49,itypea(ia),itypeb(ib))*s3*s3
      endif
      a = dexp(alpha)   ! make it positive
c end of flex alpha part....
      d1 = d(1,dmp1,rij)
      d6 = d(6,dmp6,rij)
      d8 = d(8,dmp8,rij)
      d10 = d(10,dmp10,rij)
      c6
     & =c6
     & +parab(11,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(14,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(17,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(20,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(23,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(26,itypea(ia),itypeb(ib))*(s2*s5)
      c8
     & =c8
     & +parab(12,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(15,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(18,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(21,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(24,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(27,itypea(ia),itypeb(ib))*(s2*s5)
      c10
     & =c10
     & +parab(13,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(16,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(19,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(22,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(25,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(28,itypea(ia),itypeb(ib))*(s2*s5)
c Update c6, c8, c10 with asymptotic terms asymmetric in s_A, s_B:
      c6as = 0.d0
      c8as = 0.d0
      c10as = 0.d0
      if(itypea(ia).ne.itypeb(ib)) then
      c6as
     & =c6as
     & +parab(29,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(32,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(35,itypea(ia),itypeb(ib))*(s2-s5)
      c8as
     & =c8as
     & +parab(30,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(33,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(36,itypea(ia),itypeb(ib))*(s2-s5)
      c10as
     & =c10as
     & +parab(31,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(34,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(37,itypea(ia),itypeb(ib))*(s2-s5)
       if(itypea(ia).gt.itypeb(ib)) then
        c6as = -c6as
        c8as = -c8as
        c10as = -c10as
       endif
      endif
      c6 = c6 + c6as
      c8 = c8 + c8as
      c10 = c10 + c10as
c
      if (beta.gt.0.d0) then 
      if(itypea(ia).ne.2.and.itypeb(ib).ne.2) then ! no H involved
       if(itypea(ia).eq.itypeb(ib)) then
        numtm(1) = 40
        numtm(2) =  0 !28
       else   ! test: make it the same as symmetric
        numtm(1) = 40
        numtm(2) = 28
       endif
      else
       if(itypea(ia).eq.itypeb(ib)) then
        numtm(1) = 40
        numtm(2) =  0 !28  ! no s-expansion if none of sites is H
       else
        numtm(1) = 40
        numtm(2) = 28
       endif
      endif
      numt = numtm(1) + numtm(2) + 1
c
      val0=a*dexp(-beta*rij)
      val1=val0*rij
      val2=val1*rij
      val3=val2*rij
       values(numt) = val0 + a1*val1 + a2*val2 + a3*val3
     & + d1*qa*qb/rij
     & - d6*c6/rij**6 - d8*c8/rij**8 - d10*c10/rij**10
c
c Construct planning matrix by filling values vector
c s-Symmetric terms
      if(itypea(ia).ne.2.and.itypeb(ib).ne.2) then ! no H involved
c delta a0,1,2,3 contributions to planning matrix
      values(1)  = (s1+s4)*val0
      values(2)  = (s1+s4)*val1
      values(3)  = (s1+s4)*val2
      values(4)  = (s1+s4)*val3

      values(5)  = (s2+s5)*val0
      values(6)  = (s2+s5)*val1
      values(7)  = (s2+s5)*val2
      values(8)  = (s2+s5)*val3

      values(9)  = (s3+s6)*val0
      values(10) = (s3+s6)*val1
      values(11) = (s3+s6)*val2
      values(12) = (s3+s6)*val3


      values(13) = (s1*s2+s4*s5)*val0
      values(14) = (s1*s2+s4*s5)*val1
      values(15) = (s1*s2+s4*s4)*val2
      values(16) = (s1*s2+s4*s4)*val3

      values(17) = (s2*s3+s5*s6)*val0
      values(18) = (s2*s3+s5*s6)*val1
      values(19) = (s2*s3+s5*s6)*val2
      values(20) = (s2*s3+s5*s6)*val3

      values(21) = (s1*s1+s4*s4)*val0
      values(22) = (s1*s1+s4*s4)*val1
      values(23) = (s1*s1+s4*s4)*val2
      values(24) = (s1*s1+s4*s4)*val3

      values(25) = (s2*s2+s5*s5)*val0
      values(26) = (s2*s2+s5*s5)*val1
      values(27) = (s2*s2+s5*s5)*val2
      values(28) = (s2*s2+s5*s5)*val3


      values(29) = (s1*s4)*val0
      values(30) = (s1*s4)*val1
      values(31) = (s1*s4)*val2
      values(32) = (s1*s4)*val3

      values(33) = (s2*s5)*val0
      values(34) = (s2*s5)*val1
      values(35) = (s2*s5)*val2
      values(36) = (s2*s5)*val3

      values(37) = (s3*s6)*val0
      values(38) = (s3*s6)*val1
      values(39) = (s3*s6)*val2
      values(40) = (s3*s6)*val3
c For uneqal types, add s-antisymmetric terms
      if(itypea(ia).lt.itypeb(ib)) then
      values(41) = (s1-s4)*val0
      values(42) = (s1-s4)*val1
      values(43) = (s1-s4)*val2
      values(44) = (s1-s4)*val3

      values(45) = (s2-s5)*val0
      values(46) = (s2-s5)*val1
      values(47) = (s2-s5)*val2
      values(48) = (s2-s5)*val3

      values(49) = (s3-s6)*val0
      values(50) = (s3-s6)*val1
      values(51) = (s3-s6)*val2
      values(52) = (s3-s6)*val3


      values(53) = (s1*s2-s4*s5)*val0
      values(54) = (s1*s2-s4*s5)*val1
      values(55) = (s1*s2-s4*s4)*val2
      values(56) = (s1*s2-s4*s4)*val3

      values(57) = (s2*s3-s5*s6)*val0
      values(58) = (s2*s3-s5*s6)*val1
      values(59) = (s2*s3-s5*s6)*val2
      values(60) = (s2*s3-s5*s6)*val3

      values(61) = (s1*s1-s4*s4)*val0
      values(62) = (s1*s1-s4*s4)*val1
      values(63) = (s1*s1-s4*s4)*val2
      values(64) = (s1*s1-s4*s4)*val3

      values(65) = (s2*s2-s5*s5)*val0
      values(66) = (s2*s2-s5*s5)*val1
      values(67) = (s2*s2-s5*s5)*val2
      values(68) = (s2*s2-s5*s5)*val3
      endif
      if(itypea(ia).gt.itypeb(ib)) then
      values(41) =-(s1-s4)*val0
      values(42) =-(s1-s4)*val1
      values(43) =-(s1-s4)*val2
      values(44) =-(s1-s4)*val3

      values(45) =-(s2-s5)*val0
      values(46) =-(s2-s5)*val1
      values(47) =-(s2-s5)*val2
      values(48) =-(s2-s5)*val3

      values(49) =-(s3-s6)*val0
      values(50) =-(s3-s6)*val1
      values(51) =-(s3-s6)*val2
      values(52) =-(s3-s6)*val3


      values(53) =-(s1*s2-s4*s5)*val0
      values(54) =-(s1*s2-s4*s5)*val1
      values(55) =-(s1*s2-s4*s4)*val2
      values(56) =-(s1*s2-s4*s4)*val3

      values(57) =-(s2*s3-s5*s6)*val0
      values(58) =-(s2*s3-s5*s6)*val1
      values(59) =-(s2*s3-s5*s6)*val2
      values(60) =-(s2*s3-s5*s6)*val3

      values(61) =-(s1*s1-s4*s4)*val0
      values(62) =-(s1*s1-s4*s4)*val1
      values(63) =-(s1*s1-s4*s4)*val2
      values(64) =-(s1*s1-s4*s4)*val3

      values(65) =-(s2*s2-s5*s5)*val0
      values(66) =-(s2*s2-s5*s5)*val1
      values(67) =-(s2*s2-s5*s5)*val2
      values(68) =-(s2*s2-s5*s5)*val3
      endif
      else    ! at least 1 H involved
c delta a0,1,2,3 contributions to planning matrix
      values(1)  = (s1+s4)*val0
      values(2)  = (s1+s4)*val1
      values(3)  = (s1+s4)*val2
      values(4)  = (s1+s4)*val3

      values(5)  = (s2+s5)*val0
      values(6)  = (s2+s5)*val1
      values(7)  = (s2+s5)*val2
      values(8)  = (s2+s5)*val3

      values(9)  = (s3+s6)*val0
      values(10) = (s3+s6)*val1
      values(11) = (s3+s6)*val2
      values(12) = (s3+s6)*val3


      values(13) = (s1*s2+s4*s5)*val0
      values(14) = (s1*s2+s4*s5)*val1
      values(15) = (s1*s2+s4*s4)*val2
      values(16) = (s1*s2+s4*s4)*val3

      values(17) = (s2*s3+s5*s6)*val0
      values(18) = (s2*s3+s5*s6)*val1
      values(19) = (s2*s3+s5*s6)*val2
      values(20) = (s2*s3+s5*s6)*val3

      values(21) = (s1*s1+s4*s4)*val0
      values(22) = (s1*s1+s4*s4)*val1
      values(23) = (s1*s1+s4*s4)*val2
      values(24) = (s1*s1+s4*s4)*val3

      values(25) = (s2*s2+s5*s5)*val0
      values(26) = (s2*s2+s5*s5)*val1
      values(27) = (s2*s2+s5*s5)*val2
      values(28) = (s2*s2+s5*s5)*val3


      values(29) = (s1*s4)*val0
      values(30) = (s1*s4)*val1
      values(31) = (s1*s4)*val2
      values(32) = (s1*s4)*val3

      values(33) = (s2*s5)*val0
      values(34) = (s2*s5)*val1
      values(35) = (s2*s5)*val2
      values(36) = (s2*s5)*val3

      values(37) = (s3*s6)*val0
      values(38) = (s3*s6)*val1
      values(39) = (s3*s6)*val2
      values(40) = (s3*s6)*val3
c For uneqal types, add s-antisymmetric terms
      if(itypea(ia).lt.itypeb(ib)) then
      values(41) = (s1-s4)*val0
      values(42) = (s1-s4)*val1
      values(43) = (s1-s4)*val2
      values(44) = (s1-s4)*val3

      values(45) = (s2-s5)*val0
      values(46) = (s2-s5)*val1
      values(47) = (s2-s5)*val2
      values(48) = (s2-s5)*val3

      values(49) = (s3-s6)*val0
      values(50) = (s3-s6)*val1
      values(51) = (s3-s6)*val2
      values(52) = (s3-s6)*val3


      values(53) = (s1*s2-s4*s5)*val0
      values(54) = (s1*s2-s4*s5)*val1
      values(55) = (s1*s2-s4*s4)*val2
      values(56) = (s1*s2-s4*s4)*val3

      values(57) = (s2*s3-s5*s6)*val0
      values(58) = (s2*s3-s5*s6)*val1
      values(59) = (s2*s3-s5*s6)*val2
      values(60) = (s2*s3-s5*s6)*val3

      values(61) = (s1*s1-s4*s4)*val0
      values(62) = (s1*s1-s4*s4)*val1
      values(63) = (s1*s1-s4*s4)*val2
      values(64) = (s1*s1-s4*s4)*val3

      values(65) = (s2*s2-s5*s5)*val0
      values(66) = (s2*s2-s5*s5)*val1
      values(67) = (s2*s2-s5*s5)*val2
      values(68) = (s2*s2-s5*s5)*val3
      endif
      if(itypea(ia).gt.itypeb(ib)) then
      values(41) =-(s1-s4)*val0
      values(42) =-(s1-s4)*val1
      values(43) =-(s1-s4)*val2
      values(44) =-(s1-s4)*val3

      values(45) =-(s2-s5)*val0
      values(46) =-(s2-s5)*val1
      values(47) =-(s2-s5)*val2
      values(48) =-(s2-s5)*val3

      values(49) =-(s3-s6)*val0
      values(50) =-(s3-s6)*val1
      values(51) =-(s3-s6)*val2
      values(52) =-(s3-s6)*val3


      values(53) =-(s1*s2-s4*s5)*val0
      values(54) =-(s1*s2-s4*s5)*val1
      values(55) =-(s1*s2-s4*s4)*val2
      values(56) =-(s1*s2-s4*s4)*val3

      values(57) =-(s2*s3-s5*s6)*val0
      values(58) =-(s2*s3-s5*s6)*val1
      values(59) =-(s2*s3-s5*s6)*val2
      values(60) =-(s2*s3-s5*s6)*val3

      values(61) =-(s1*s1-s4*s4)*val0
      values(62) =-(s1*s1-s4*s4)*val1
      values(63) =-(s1*s1-s4*s4)*val2
      values(64) =-(s1*s1-s4*s4)*val3

      values(65) =-(s2*s2-s5*s5)*val0
      values(66) =-(s2*s2-s5*s5)*val1
      values(67) =-(s2*s2-s5*s5)*val2
      values(68) =-(s2*s2-s5*s5)*val3
      endif
c
      endif   ! if H involved
 
      else
c
      numt=1
      numtm(1) = 0
      numtm(2) = 0
      values(numt) =
     & + d1*qa*qb/rij
     & - d6*c6/rij**6 - d8*c8/rij**8 - d10*c10/rij**10
c
      end if 
c
      return
      endif
      end
c
c============================================================================
c
c     calculate the damping factor (small R correct)
c
      function d(n,beta,r)
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
        write(66,*) 'No convergence in d'
  111 continue
      d=d*dexp(-br)
      endif
      return
      end
c============================================================================
c
c The parameters of the fit are read in this subroutine.
c
      subroutine data1(param,parab,nsitea,nsiteb,c)
      implicit real*8 (a-h,o-z)
c
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
      character*8 label(9)
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension c(1000)
      integer omp_get_num_threads,omp_get_thread_num

      a0=0.529177249d0

c      ns=omp_get_num_threads()
c      id=omp_get_thread_num()
      
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
      nsitea=8
      nsiteb=8

c
c Read in the one-site nonlinear parameters and opt indicators..
c 
c nparm used only here
      read(55,*) nparm
      do i=1,nparm
       read(55,*) ityp, inumpar, val, indopt
       param(inumpar,ityp) = val
c Convert charges to the proper units 
c and DELTA charge!!! GM 3/29/05
c Switch off for the total charge density variant.... (everything
c in atomic units)
       if(inumpar.le.9) then
        param(inumpar,ityp) = 18.22262373d0*param(inumpar,ityp)
       endif
      end do
c
c Read in the two-site nonlinear parameters and opt indicators..
c
c nparab used only here
      read(55,*) nparab
      do i=1,nparab
       read(55,*) ityp1,ityp2,inumpar, val, indopt
       parab(inumpar,ityp1,ityp2) = val
       parab(inumpar,ityp2,ityp1) = val
      end do
c
c Read the type of the potential form to be used, inear coeff. reading
c indicator, and the optimization parameters....
c
c w wersji ostatecznej ustawic:
c ntpot=124161
c idonl=1 (w procedurze wolajacej)
c wartosci iopt, iweight, iasdone niewazne - mozna pominac
c oraz dodac do commonu dataoffit
      read(55,*) ntpot, idonl, iopt, iweight, iasdone
c      ntpot=124161
c the following parameters are not necessary for evaluation of the interaction
c energy value, but are left to use the standard input - can be removed
c in the final "paper" version.
      read(55,*) TOLF,TOLR,ANOISE,TIMLIM,LINP,LOUT,ALFST,IPR,SAFETL
      read(55,*) R_0,isyst,npowers       ! to be removed in the final version
      read(55,*) RCOND                   ! to be removed in the final version
c      ntpot=124161
c      iasdone=1

c
        read(55,*)numlin
        do i=1,numlin
         read(55,*) c(i)
        end do

       return
       end
c============================================================================
c
c A simple routine to calculate the dipole-dipole induction,
c now in the flexible case
c units: charges and polarizability calculated in au; energy returned in kcal/mol
c
      subroutine dipind(iaa,ibb,R,energy,
     .    sa,sb,param,parab,siteat,sitebt,nsitea,nsiteb,itypea,itypeb)
      implicit real*8 (a-h,o-z)
cccccc Parameters taken from potparts
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension sa(3,-1:naamax),sb(3,-1:nbbmax)
      dimension siteat(3,nsitemax,-1:naamax),
     .          sitebt(3,nsitemax,-1:nbbmax)
      dimension itypea(nsitemax),itypeb(nsitemax)
cccccc end of parameters taken from potparts
c      dimension dma_unrot(3), dmb_unrot(3), dma(3), dmb(3), u(3)
      dimension dma(3), dmb(3), u(3)
      a0=0.529177249d0
      har2kcal=627.510d0
c
c Calculate static dipole moments of A and B in these orientations and
c Let's place the polarizability on Oxygen, just
c for simplicity
c
c monomer A
c
      s1=sa(1,iaa)
      s2=sa(2,iaa)
      s3=sa(3,iaa)
c      write(31,*) iaa,s1,s2,s3
      signa=1.0d0
      do i=1,3
c       dma_unrot(i)=0.0d0  !for checking
       dma(i)=0.0d0
      end do
      polisa=0.0d0
c
      do 11 ia=1,nsitea
      if (ia .eq. 3) signa=-1.0d0
      s3=signa*s3
c Make basis functions even in s6 if site is NOT H
c     if(itypea(ia).ne.2) s3=dabs(s3)
      qa
     & =param(1,itypea(ia))         !already in units of [kcal/mol & Ang]/Ang^n ; multiply by Ang^n to get bohr^3
     & +param(2,itypea(ia))*s1
     & +param(3,itypea(ia))*s2
     & +param(4,itypea(ia))*s3
     & +param(5,itypea(ia))*s1*s2
     & +param(6,itypea(ia))*s2*s3
     & +param(7,itypea(ia))*s1*s1
     & +param(8,itypea(ia))*s2*s2
     & +param(9,itypea(ia))*s3*s3
      qa=qa/18.22262373d0           !convert back to atomic units
c      write(31,101)
c     &iaa,ia,param(1,itypea(ia))/18.22262373d0,qa
      do i=1,3
c       dma_unrot(i)=dma_unrot(i)+qa*(sitea(i,ia,iaa))/a0
       dma(i)=dma(i)+qa*(siteat(i,ia,iaa))/a0
      end do
      if (ia.eq.1) then !Oxygen only
      polisa
     & =param(10,itypea(ia))
     & +param(11,itypea(ia))*s1
     & +param(12,itypea(ia))*s2
     & +param(13,itypea(ia))*s3
     & +param(14,itypea(ia))*s1*s2
     & +param(15,itypea(ia))*s2*s3
     & +param(16,itypea(ia))*s1*s1
     & +param(17,itypea(ia))*s2*s2
     & +param(18,itypea(ia))*s3*s3  !already in units of bohr^3/Ang^n ; multiply by Ang^n to get bohr^3
      end if
11    continue
c      write(31,201) iaa, (dma_unrot(i),i=1,3)
c      write(31,301) iaa, (dma(i),i=1,3)
c      write(31,*) 'iaa, polisa:', iaa, polisa
c
101   format('iaa,ia,qa_0,qa:',
     &2(1x,i6),2(1x,f12.7))
201   format('iaa,dma_unrot(i):',
     &1(1x,i6),3(1x,f12.7))
301   format('iaa,dma(i):',
     &1(1x,i6),3(1x,f12.7))
c
c monomer B
c
      s4=sb(1,ibb)
      s5=sb(2,ibb)
      s6=sb(3,ibb)
c      write(32,*) ibb,s4,s5,s6
      signb=1.0d0
      do i=1,3
c       dmb_unrot(i)=0.0d0  !for checking
       dmb(i)=0.0d0
      end do
      polisb=0.0d0
c
      do 12 ib=1,nsiteb
      if (ib .eq. 3) signb=-1.0d0
      s6=signb*s6
c Make basis functions even in s6 if site is NOT H
c     if(itypeb(ib).ne.2) s6=dabs(s6)
      qb
     & =param(1,itypeb(ib))         !already in units of [kcal/mol & Ang]/Ang^n ; multiply by Ang^n to get bohr^3
     & +param(2,itypeb(ib))*s4
     & +param(3,itypeb(ib))*s5
     & +param(4,itypeb(ib))*s6
     & +param(5,itypeb(ib))*s4*s5
     & +param(6,itypeb(ib))*s5*s6
     & +param(7,itypeb(ib))*s4*s4
     & +param(8,itypeb(ib))*s5*s5
     & +param(9,itypeb(ib))*s6*s6
      qb=qb/18.22262373d0           !convert back to atomic units
c      write(32,102)
c     &ibb,ib,param(1,itypeb(ib))/18.22262373d0,qb
      Rtemp=0.0
      do i=1,3
       if(i.eq.3)Rtemp=R
c       dmb_unrot(i)=dmb_unrot(i)+qb*(siteb(i,ib,ibb))/a0
       dmb(i)=dmb(i)+qb*(sitebt(i,ib,ibb)-Rtemp)/a0
      end do
      if (ib.eq.1) then !Oxygen only
      polisb
     & =param(10,itypeb(ib))
     & +param(11,itypeb(ib))*s4
     & +param(12,itypeb(ib))*s5
     & +param(13,itypeb(ib))*s6
     & +param(14,itypeb(ib))*s4*s5
     & +param(15,itypeb(ib))*s5*s6
     & +param(16,itypeb(ib))*s4*s4
     & +param(17,itypeb(ib))*s5*s5
     & +param(18,itypeb(ib))*s6*s6  !already in units of bohr^3/Ang^n ; multiply by Ang^n to get bohr^3
      end if
12    continue
c      write(32,202) ibb, (dmb_unrot(i),i=1,3)
c      write(32,302) ibb, (dmb(i),i=1,3)
c      write(32,*) 'ibb, polisb:', ibb, polisb
c
102   format('ibb,ib,qb_0,qb:',2(1x,i6),2(1x,f12.7))
202   format('ibb,dmb_unrot(i):',1(1x,i6),3(1x,f12.7))
302   format('ibb,dmb(i):',1(1x,i6),3(1x,f12.7))
c
c Calculate the distance between the O's (pol centers)
c
      dlen = 0.d0
      do i=1,3
       pom = sitebt(i,1,ibb) - siteat(i,1,iaa)
       dlen = dlen + pom*pom
      end do
      dlen = dsqrt(dlen)
c
c Compute the damping factor for induction
c
      dmpind = d(6,parab(10,1,1),dlen)
c
      dlen = dlen**(-3.d0)   ! inverse cube of the distance
c
c u(*): field of A on B and the induction energy
c
      call TTTprod(siteat(1,1,iaa),sitebt(1,1,ibb),dma,dlen,u)
      energy_a_on_b = polisa*scalp(u,u)
c
c u(*): field of B on A and the induction energy
c
      call TTTprod(siteat(1,1,iaa),sitebt(1,1,ibb),dmb,dlen,u)
      energy_b_on_a = polisb*scalp(u,u)
c
      energy = energy_a_on_b + energy_b_on_a
      energy = -0.5d0*(a0**6)*har2kcal*energy*dmpind
c
      return
      end
c=============================================================================
c
c Calculate the vector V resulting from the action
c of the dipole propagator tensor T_ij on U
c Ri, Rj - position vectors of molecules i and j, respectively
c rij = |Ri - Rj|^(-3)
c
        subroutine TTTprod(Ri,Rj,u,rij,v)
        implicit real*8 (a-h,o-z)
        dimension Ri(3),Rj(3),u(3),v(3)
c
c      goto 777
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
 777   continue
        return
        end
c=============================================================================
c
        function scalp(a,b)
        implicit real*8 (a-h,o-z)
        dimension a(3),b(3)
        scalp = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
        return
        end 
c=============================================================================
c In this procedure positions of sites are established on the base of 
c the cartesian coordinates of atoms. 
c Written on the base of data from the subroutine getgeo
c
c Also the sa and itypea matrices are filled in this routine

      subroutine set_sites(carta,sitea,sa,itypea,nsitea)
      implicit real*8 (a-h,o-z)

      parameter (lay_claude=0)
      parameter (nsitemax=8)
c
      dimension sitea(3,nsitemax),sa(3)
      dimension itypea(nsitemax)
      dimension carta(3,3),cartb(3,3)
      dimension v1(3),vn1(3),v2(3),vn2(3),v(3),vn(3),vv(3),vsm(3)
      dimension vb(3),vp(3),vd1a(3),vd1b(3),vd2a(3),vd2b(3)

      a0=0.529177249d0
      zero=0.d0
      r0_ang=0.9716257d0
      theta0_deg=104.69d0

      iaa=-1
      ibb=-1
c           
c Info regarding site positions for an undistorted monomer with approximation 
c that amO=16 amu, amH=1 amu 
c  sig   is the usual O site-COM site separation
c  sig2  is the O sep along the z coordinate of Bunny 1
c  sig3  is the y coordinate of Bunny 1
c  sig4  is the O sep along the z coordinate of Bunny 2
c  sig5  is the y coordinate of Bunny 2
      sig=  1.24631924913843d-01                 ! in bohr
      sig2= 0.371792435d0                        !
      sig3= 0.2067213d0                          !
      sig4= 0.125368076d0                        !
      sig5= 0.2d0                                !

      shift=9.01563628739252d-4  !accounts for shift to precise masses (in bohr)

      pi=dacos(-1.d0)
      rad2d=180.d0/pi

c Distances in carta are given in bohr, so we will use bohr as the default unit
c They should be changed to Angstroms at the end
c--- O1
      sitea(1,1)=carta(1,1)
      sitea(2,1)=carta(1,2)
      sitea(3,1)=carta(1,3)
c--- H1
      sitea(1,2)=carta(2,1)
      sitea(2,2)=carta(2,2)
      sitea(3,2)=carta(2,3)
c--- H2
      sitea(1,3)=carta(3,1)
      sitea(2,3)=carta(3,2)
      sitea(3,3)=carta(3,3)

c vn1 is a normalized vector pointing from O to H1
      v1(1)=carta(2,1)-carta(1,1)
      v1(2)=carta(2,2)-carta(1,2)
      v1(3)=carta(2,3)-carta(1,3)
      xnv1=dsqrt(v1(1)**2+v1(2)**2+v1(3)**2)
      vn1(1)=v1(1)/xnv1
      vn1(2)=v1(2)/xnv1
      vn1(3)=v1(3)/xnv1
      xnv1_ang=xnv1*a0
c      write(44,*) "xnv1_ang=",xnv1_ang," AA"

c vn2 is a normalized vector pointing from O to H2
      v2(1)=carta(3,1)-carta(1,1)
      v2(2)=carta(3,2)-carta(1,2)
      v2(3)=carta(3,3)-carta(1,3)
      xnv2=dsqrt(v2(1)**2+v2(2)**2+v2(3)**2)
      vn2(1)=v2(1)/xnv2
      vn2(2)=v2(2)/xnv2
      vn2(3)=v2(3)/xnv2
      xnv2_ang=xnv2*a0
c      write(44,*) "xnv2_ang=",xnv2_ang," AA"

c components of the bisector (vb) vector (normalized):
      v(1)=vn1(1)+vn2(1)
      v(2)=vn1(2)+vn2(2)
      v(3)=vn1(3)+vn2(3)
      xnv=dsqrt(v(1)**2+v(2)**2+v(3)**2)
      vb(1)=v(1)/xnv
      vb(2)=v(2)/xnv
      vb(3)=v(3)/xnv

c a normalized vector, perpendicular (vp) to the surface defined by 
c the v1 and v2 vectors (or vn1 and vn2)
      v(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v(3)=v1(1)*v2(2)-v1(2)*v2(1)
      xn=dsqrt(v(1)**2+v(2)**2+v(3)**2)
      vp(1)=v(1)/xn
      vp(2)=v(2)/xn
      vp(3)=v(3)/xn

c set bunny_ratio_A
      r0=r0_ang/a0
      theta0=theta0_deg/rad2d
      cta=dcos(0.5d0*theta0)
      prodv1vb=v1(1)*vb(1)+v1(2)*vb(2)+v1(3)*vb(3)
      prodv2vb=v2(1)*vb(1)+v2(2)*vb(2)+v2(3)*vb(3)
      IF (LAY_CLAUDE .EQ. 0) THEN
        bunny_ratio_A=(0.5d0*(prodv1vb+prodv2vb))/(r0*cta)
c        write(44,*) "bunny_ratio=",bunny_ratio_A
      END IF
      IF (LAY_CLAUDE .EQ. 1) THEN
        bunny_ratio_A=1.0d0
      END IF

c--- Bunny1 1 -- charged one
      vd1a(1)=sig3*vp(1)+sig2*vb(1)*bunny_ratio_A
      vd1a(2)=sig3*vp(2)+sig2*vb(2)*bunny_ratio_A
      vd1a(3)=sig3*vp(3)+sig2*vb(3)*bunny_ratio_A
      sitea(1,4)=carta(1,1)+vd1a(1)
      sitea(2,4)=carta(1,2)+vd1a(2)
      sitea(3,4)=carta(1,3)+vd1a(3)
c--- Bunny1 2 -- charged one
      vd1b(1)=-sig3*vp(1)+sig2*vb(1)*bunny_ratio_A
      vd1b(2)=-sig3*vp(2)+sig2*vb(2)*bunny_ratio_A
      vd1b(3)=-sig3*vp(3)+sig2*vb(3)*bunny_ratio_A
      sitea(1,5)=carta(1,1)+vd1b(1)
      sitea(2,5)=carta(1,2)+vd1b(2)
      sitea(3,5)=carta(1,3)+vd1b(3)
c--- Bunny2 1 -- exp-type one
      vd2a(1)=sig5*vp(1)-sig4*vb(1)*bunny_ratio_A
      vd2a(2)=sig5*vp(2)-sig4*vb(2)*bunny_ratio_A
      vd2a(3)=sig5*vp(3)-sig4*vb(3)*bunny_ratio_A
      sitea(1,6)=carta(1,1)+vd2a(1)
      sitea(2,6)=carta(1,2)+vd2a(2)
      sitea(3,6)=carta(1,3)+vd2a(3)
c--- Bunny2 2 -- exp-type one
      vd2b(1)=-sig5*vp(1)-sig4*vb(1)*bunny_ratio_A
      vd2b(2)=-sig5*vp(2)-sig4*vb(2)*bunny_ratio_A
      vd2b(3)=-sig5*vp(3)-sig4*vb(3)*bunny_ratio_A
      sitea(1,7)=carta(1,1)+vd2b(1)
      sitea(2,7)=carta(1,2)+vd2b(2)
      sitea(3,7)=carta(1,3)+vd2b(3)

c calculate COM position with Robert''s masses
      xm16=15.994915d0
      xm1=1.007825d0
      sm=xm16+2.0d0*xm1
      vsm(1)=(xm16*carta(1,1)+xm1*carta(2,1)+xm1*carta(3,1))/sm
      vsm(2)=(xm16*carta(1,2)+xm1*carta(2,2)+xm1*carta(3,2))/sm
      vsm(3)=(xm16*carta(1,3)+xm1*carta(2,3)+xm1*carta(3,3))/sm
c and then shift the position of the com-site to the position of COM calculated
c with the approximated masses.
c--- com
      sitea(1,8)= vsm(1) - shift*vb(1)
      sitea(2,8)= vsm(2) - shift*vb(2)
      sitea(3,8)= vsm(3) - shift*vb(3)
c      sss=sqrt(vsm(1)**2+vsm(2)**2+vsm(3)**2)
c      write(44,*) "|vsm|=",sss

c transform to Angstroms...
      do ia=1,nsitea
       do i=1,3
        sitea(i,ia) = sitea(i,ia)*a0
       end do
      end do

c As in the subroutine "monomer":
c distorts monomer A (or B) and writes the distorted symmetry coordinates sa
c (passed back as sa or sb) 

      sprod=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
      ccos=sprod/(xnv1*xnv2)
      theta1=dacos(ccos)
      theta1_deg=theta1*rad2d
c      write(44,*) "theta1_deg=",theta1_deg
      dsqrt2=dsqrt(2.0d0)
      sa(1)=( (xnv1_ang-r0_ang) + (xnv2_ang-r0_ang) )/dsqrt2
      sa(2)=dsqrt(xnv1_ang*xnv2_ang)*(theta1_deg-theta0_deg)/rad2d
      sa(3)=( (xnv1_ang-r0_ang) - (xnv2_ang-r0_ang) )/dsqrt2
     
      itypea(1)=1
      itypea(2)=2
      itypea(3)=2
      itypea(4)=3
      itypea(5)=3
      itypea(6)=4
      itypea(7)=4
      itypea(8)=5

      return
      end
c==============================================================================



      SUBROUTINE POTS(V,Q1,Q2,THETA)
C     Potential PJT2 due Polyansky, Jensen and Tennyson, 
C     J. Chem. Phys., 105, 6490-6497 (1996)
C     Update of Polyansky, Jensen and Tennyson, J Chem Phys 101, 7651 (1994))
C     Units: Hartree and Bohr
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     RZ = OH equilibrium value
C     RHO = equilibrium value of pi - bond angle(THETA)
      DATA TOANG/0.5291772/, CMTOAU/219474.624/
      DATA X1/1.0/
      DATA RHO1    /    75.50035308/
      DATA FA1     /     .00000000/
      DATA FA2     /18902.44193433/
      DATA FA3     /  1893.99788146/
      DATA FA4     /  4096.73443772/
      DATA FA5     /-1959.60113289/
      DATA FA6     /  4484.15893388/
      DATA FA7     /  4044.55388819/
      DATA FA8     / -4771.45043545/
      DATA FA9     /     0.00000000/
      DATA FA10    /     0.00000000/
      DATA RZ    /     .95792059/
      DATA A     /    2.22600000/
      DATA F1A1    /  -6152.40141181/
      DATA F2A1    / -2902.13912267/
      DATA F3A1    / -5732.68460689/
      DATA F4A1    /  953.88760833/
      DATA F11     / 42909.88869093/
      DATA F1A11   /  -2767.19197173/
      DATA F2A11   /  -3394.24705517/
      DATA F3A11   /     .00000000/
      DATA F13     /  -1031.93055205/
      DATA F1A13   /  6023.83435258/
      DATA F2A13   /     .00000000/
      DATA F3A13   /     .00000000/
      DATA F111    /     .00000000/
      DATA F1A111  /   124.23529382/
      DATA F2A111  /  -1282.50661226/
      DATA F113    /  -1146.49109522/
      DATA F1A113  /  9884.41685141/
      DATA F2A113  /  3040.34021836/ 
      DATA F1111   /  2040.96745268/
      DATA FA1111  /  .00000000/
      DATA F1113   /  -422.03394198/
      DATA FA1113  /-7238.09979404/
      DATA F1133   /     .00000000/
      DATA FA1133  /     .00000000/
      DATA F11111  / -4969.24544932/
      DATA f111111/  8108.49652354/
      DATA F71   /  90.00000000/
 
      data c1/50.0/,c2/10.0/,beta1/22.0/,beta2/13.5/,gammas/0.05/,
     *     gammaa/0.10/,delta/0.85/,rhh0/1.40/
                 RHO=RHO1*3.141592654/180.000000000
      fa11=0.0
      f1a3=f1a1
      f2a3=f2a1
      f3a3=f3a1
      f4a3=f4a1
      f33=f11
      f1a33=f1a11
      f2a33=f2a11
      f333=f111
      f1a333=f1a111
      f2a333=f2a111
      f133=f113
      f1a133=f1a113
      f2a133=f2a113
      f3333=f1111
      fa3333=fa1111
      f1333=f1113
      fa1333=fa1113
      f33333=f11111
      f333333 =f111111
      f73     =f71
 
C     Find value for DR and DS
      DR = TOANG*Q1 - RZ
      DS = TOANG*Q2 - RZ
 
C     Transform to Morse coordinates
      Y1 = X1 - EXP(-A * DR)
      Y3 = X1 - EXP(-A * DS)
 
C     transform to Jensens angular coordinate
      CORO = DCOS(THETA) + DCOS(RHO)
 
C     Now for the potential
      V0=(FA2+FA3*CORO+FA4*CORO**2+FA6*CORO**4+FA7*CORO**5)*CORO**2
      V0=V0+(FA8*CORO**6+FA5*CORO**3+FA9*CORO**7+FA10*CORO**8 )*CORO**2
      V0=V0+(                                    FA11*CORO**9 )*CORO**2
      FE1= F1A1*CORO+F2A1*CORO**2+F3A1*CORO**3+F4A1*CORO**4
      FE3= F1A3*CORO+F2A3*CORO**2+F3A3*CORO**3+F4A3*CORO**4
      FE11= F11+F1A11*CORO+F2A11*CORO**2
      FE33= F33+F1A33*CORO+F2A33*CORO**2
      FE13= F13+F1A13*CORO
      FE111= F111+F1A111*CORO+F2A111*CORO**2
      FE333= F333+F1A333*CORO+F2A333*CORO**2
      FE113= F113+F1A113*CORO+F2A113*CORO**2
      FE133= F133+F1A133*CORO+F2A133*CORO**2
      FE1111= F1111+FA1111*CORO
      FE3333= F3333+FA3333*CORO
      FE1113= F1113+FA1113*CORO
      FE1333= F1333+FA1333*CORO
      FE1133=       FA1133*CORO
      FE11111=F11111
      FE33333=F33333
      FE111111=F111111
      FE333333=F333333
      FE71    =F71
      FE73    =F73
      V   = V0 +  FE1*Y1+FE3*Y3
     1         +FE11*Y1**2+FE33*Y3**2+FE13*Y1*Y3
     2         +FE111*Y1**3+FE333*Y3**3+FE113*Y1**2*Y3
     3         +FE133*Y1*Y3**2
     4         +FE1111*Y1**4+FE3333*Y3**4+FE1113*Y1**3*Y3
     5         +FE1333*Y1*Y3**3+FE1133*Y1**2*Y3**2
     6         +FE11111*Y1**5+FE33333*Y3**5
     7         +FE111111*Y1**6+FE333333*Y3**6
     8         +FE71    *Y1**7+FE73    *Y3**7
C     modification by Choi & Light, J. Chem. Phys., 97, 7031 (1992).
      sqrt2=sqrt(2.0)
      xmup1=sqrt2/3.0+0.5
      xmum1=xmup1-x1
      term=2.0*xmum1*xmup1*q1*q2*cos(theta)
      r1=toang*sqrt((xmup1*q1)**2+(xmum1*q2)**2-term)
      r2=toang*sqrt((xmum1*q1)**2+(xmup1*q2)**2-term)
      rhh=sqrt(q1**2+q2**2-2.0*q1*q2*cos(theta))
      rbig=(r1+r2)/sqrt2
      rlit=(r1-r2)/sqrt2
 
      alpha=(x1-tanh(gammas*rbig**2))*(x1-tanh(gammaa*rlit**2))
      alpha1=beta1*alpha
      alpha2=beta2*alpha
      drhh=toang*(rhh-delta*rhh0)
      DOLEG=     (1.4500-THETA)
C     IF (THETA.LE.0.64  ) V=0.1E17
C     IF((DR.LE.-0.4).AND.(THETA.LE.1.1)) V=0.1E17
C     IF((DS.LE.-0.4).AND.(THETA.LE.1.1)) V=0.1E17
C     IF (DS.LE. 0.0  ) V=0.1E17
      v = v + c1*exp(-alpha1*drhh) + c2*exp(-alpha2*drhh)
 
C     Convert to Hartree
      V=V/CMTOAU
      RETURN
      END
	  	  