      SUBROUTINE USER_DEFINED_PES(V,R,rvib,rvib2,alpha,beta,gamma,
     &										   aalpha,bbeta,ggamma)
!     INPUT:  R - distance betweenn COMs, rvib - vibrational coordinate
! 		      alpha,beta,gamma - Euler's angles of the first molecule
!   		  aalpha,bbeta, ggamma - Euler's angles of the second molecule
!     OUTPUT: V - value of the potential
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 V,R,alpha,beta,gamma,aalpha,bbeta,ggamma,t,rvib,rvib2
      REAL*8 R1,COG, phi
	  real*8, parameter :: PI = 4.d0*datan(1.d0)
      DIMENSION rbohr(3)
	
	  CALL pot_init(i)
	  
	  !rvib = 2.1359d0
	  !rvib2 = 1.4011d0
	  phi = PI - gamma

	  CALL pot_h2co(R, rvib, rvib2, beta, bbeta, phi, V)
      
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


cc June 2, 2016
Cc before calling potential pot_h2co, call the initialization subroutine first to read
cc expansion coefficients
cc
	   subroutine pot_init(i)
	   implicit none
	   integer:: i 
	   !integer, intent(in)::i_pot 
	   integer,parameter::ncoff=502
	   integer, parameter :: wp=selected_real_kind(12,300)
	   real(kind=wp)::cof(ncoff)
	   common/pot_exp_coeff/cof
	   
	   open (11, file='fit_coefs.dat')

       do i=1, ncoff
        read(11, *) cof(i)
       end do

        close (11)

        return
         end


!! --------------------------------------------------------------------
!! distance in bohr,  angles in radian
!! potential Vpot is returned in hartree 
!! 
       subroutine pot_h2co(R, R_CO, R_HH, theta_1, theta_2, phi, V_int)
       implicit none
       real, parameter:: conv = 219474.631      ! hartree to 1/cm.
       real, parameter:: mass_C = 12.00000, mass_O = 15.99940
       integer, parameter :: wp=selected_real_kind(12,300)
       integer::i, j, k
       real(kind=wp):: cc(4,3), Vfit 
       real(kind=wp):: theta_1, theta_2, phi, R_HH, R_CO, R, Rc, Ro, Vpot
       real(kind=wp):: COPot, H2Pot, V_int

        Rc = R_CO * mass_O/(mass_C + mass_O)
        Ro = R_CO * mass_C/(mass_C + mass_O)
!!
!! Cartesian Coordinates for C, O, H, H atoms
!!
        cc(1,1) = 0.0d0
        cc(1,2) = -Rc * dsin(theta_1)
        cc(1,3) = -Rc*dcos(theta_1)

        cc(2,1) = 0.0d0
        cc(2,2) =  Ro * dsin(theta_1)
        cc(2,3) =  Ro*dcos(theta_1)
                   
        cc(3,1) = -R_HH * dsin(theta_2) * 0.5 * dsin(phi)
        cc(3,2) = -R_HH * dsin(theta_2) * 0.5 * dcos(phi)
        cc(3,3) = R - R_HH * dcos(theta_2) * 0.5

        cc(4,1) = R_HH * dsin(theta_2) * 0.5 * dsin(phi) 
        cc(4,2) = R_HH * dsin(theta_2) * 0.5 * dcos(phi)
        cc(4,3) = R + R_HH * dcos(theta_2) * 0.5

         call fith2co(cc, Vfit)

          V_int = Vfit   !!! in hartree (a.u.)
 
       return
      END

!!!!-----------------------------------------------------------
!!!!-----------------------------------------------------------

       subroutine fith2co(cood, vfit)
       implicit none
       real,parameter::power=6
       integer,parameter::ncoff=502
       integer, parameter :: wp=selected_real_kind(12,300)
       real,parameter::alpha=0.5
       real(kind=wp)::cof(ncoff),cood(4,3),y(6),bas(ncoff),getpot, vfit
       integer::i
       character(len=2)::symb
       common/pot_exp_coeff/cof

!       open(112,status='old',file='fit_coef.dat')
!       read(112,*)
!       do i=1,ncoff
!       read(112,*) cof(i)
!       end do
  
       call edis(y(1),cood(3,:),cood(4,:))
       call edis(y(2),cood(2,:),cood(3,:))
       call edis(y(3),cood(1,:),cood(3,:))
       call edis(y(4),cood(1,:),cood(4,:))
       call edis(y(5),cood(2,:),cood(4,:))
       call edis(y(6),cood(1,:),cood(2,:))

       call basis(bas,y)
       getpot=0.0d0
       do i=1,ncoff
        getpot=getpot+cof(i)*bas(i)
       end do
       vfit = getpot
       close(112)
         return
          end

       subroutine edis(edist,atomi,atomj)
       real,parameter::alpha=0.5
       real(kind=8)::atomi(3),atomj(3)
       real(kind=8)::dis
       real(kind=8),intent(out)::edist
       dis=sqrt((atomi(1)-atomj(1))**2+(atomi(2)-atomj(2))**2
     & +(atomi(3)-atomj(3))**2)
       edist=exp(-alpha*dis)
       return
       end subroutine edis

        subroutine basis(bas,y)
        integer,parameter::ncoff=502, power=6
        integer m,n,p,q,i,j,k,x,z,s,dex(919,6)
        real(kind=8)::bas(ncoff)
        real(kind=8)::y(6)
         k=1
          do i=0,power
           do j=0,power-i
            do m=0,power-i-j
             do n=0,power-i-j-m
              do p=0,power-i-j-m-n
               do q=0,power-i-j-m-n-p
                 z=0
                 dex(k,1)=i
                 dex(k,2)=j
                 dex(k,3)=m
                 dex(k,4)=n
                 dex(k,5)=p
                 dex(k,6)=q
                  do s=1,k-1
                  if(dex(k,1)==dex(s,1).and.dex(k,6)==dex(s,6)
     $ .and.dex(k,5)==dex(s,2).and.dex(k,3)==dex(s,4)
     $ .and.dex(k,4)==dex(s,3).and.dex(k,5)==dex(s,2)) then
                         z=1
                         exit
                       end if
                    end do
                    if(z==1) cycle
                    bas(k)=(y(1)**i)*(y(6)**q)*((y(2)**j)
     $       *(y(3)**m)*(y(4)**n)*(y(5)**p)+(y(2)**p)
     $   *(y(3)**n)*(y(4)**m)*(y(5)**j))
                    k=k+1
                 end do  
               end do
              end do
            end do
           end do
          end do
          return
          end subroutine
