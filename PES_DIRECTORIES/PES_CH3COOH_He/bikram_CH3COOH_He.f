	    subroutine bk_Met_He(energ,rr,eua,eub)
! A simple program to calculate the total interaction energy
! from the fit.....
! Modified by KSz on 3/11/2010 to remove excessive printouts and print
! a file readable by ssfit. 

        implicit real*8 (a-h,o-z)
		real*8 energ, rr, eua(3), eub(3)
		logical :: init = .true.
        parameter (nsitemax=32)
        dimension rpt(6),oa(3),ob(3), valder(8)
        dimension forceab(3),torquea(3),torqueb(3)
        dimension coords(7)
	    dimension en(20),dsp(20),dind(20)
        common/temp_sites/ siteat(3,nsitemax), sitebt(3,nsitemax),
     1                   tranA(3,3), tranB(3,3)
! The following commons needed for the emp2 potential\
      common /unitc/ aerg,aev,acm1,atem,akm,amev,amu,ams,atim,aps,akjm
      common/param14/ipot,imod
      common/numbs/nmol
! End of the emp2 commons
        pi = dacos(-1.d0)
        d2rad = pi/180.d0

! read in the data
		if(init) then
        imode = -1
        call poten(r,oa,ob,tot,forceab,torquea,torqueb,imode)
		init = .false.
		print*,'Data Initilization Done'
		end if
!
!       call rdmom
!       call rdcoef
!       call fct(40)
!
!        do iq = 1,10000000


!       read distance (A) and Euler angles of monomers (in deg)

!        read(5,*,end=200) r, oa(2),oa(3),(ob(i),i=1,3)

!       r,  COM separation,  is in angstroms.

!       r=2.d0
!       do iir=1,3
!       r=r+0.5d0
!
!       oa(1)=0.0d0
!       ob(1)=0.0d0
!       ob(2)=0.0d0
!       ob(3)=0.0d0
!                              step of 10 degrees;
!       do ib=0,18
!       oa(2) = ib*10.0d0
!       do ig=0,18
!       oa(3) = ig*10.0d0

		r = rr
		do i=1,3
        oa(i) = eua(i)
        end do
		do i=1,3
        ob(i) = eub(i)
        end do

        do i=1,3
         oa(i) = d2rad*oa(i)
         ob(i) = d2rad*ob(i)
        end do
! Uncomment the following three lines to
! calculate the potential in radial loop
!        r0 = 2.00d0
!        do iir=300,0,-1
!        r = r0 + iir*0.02d0
!       
!       special case of zero elst
!       call elst(r,oa,ob,en)
!       ener=0.0d0
!       do i=1,20
!         ener=ener+en(i)
!         write(*,*)'1/R**',i,en(i)
!       end do
!       call dispind(r,oa,ob,dsp,dind)
!       write(*,*) dsp
!       write(*,*) dind

        imode = 4
        call poten(r,oa,ob,value1,forceab,torquea,torqueb,imode)
		energ = value1
!       disp=0.0d0
!       din=0.0d0
!       do i=1,20
!         disp=disp+dsp(i)
!         
!         din=din+dind(i)
!       end do
!       
!       write(*,*)"disp",disp
!       write(*,*)"din",din
!       
!       write(*,*)"disp",(i,dsp(i),i=1,20)
!       write(*,*)"ind",(i,dind(i),i=1,20)
!        do i=1,3
!         oa(i) = oa(i)/d2rad
!         ob(i) = ob(i)/d2rad
!        end do

!       write(6,*)'Geometry',r,oa(2),oa(3),ob
!        zero=0.0d0
!        ener=0.0d0
!        sum=ener+disp+din
!         write(6,'(A,6F10.3,4E18.10)') 'en',r, oa(2),oa(3),
!     *   ob(1),ob(2),ob(3),value1
!        write(6,'(A,6F10.3,4E18.10)') 'en',r, oa(2),oa(3),
!    *   ob(1),ob(2),ob(3),ener,disp,din,sum
!        write(1,'(6F8.2,E18.10,7f4.1)') 
!    *                         r, oa(2),oa(3),ob(1),ob(2),ob(3),
!    *                         sum,zero,zero,zero,zero,zero,zero,zero
!        write(6,'(A,6F14.8,5E18.10)') 'en 1/R3',r, oa(2),oa(3),
!    *   ob(1),ob(2),ob(3),en(3)
! Uncomment the following three lines to
! calculate the potential in radial loop
!         if(value1.ge.20.d0) go to 15
!        end do
! 15     continue
!
!       end do
!       end do
!        end do   ! loop over read-in coordinates
! 200    continue
		close(5)
		return
        end subroutine
