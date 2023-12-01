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
!----------! 

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
!     USER MUST INSERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE IS SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:
!      STOP "ERROR: USER_DEFINED_PES IS NOT SUPPLIED"
	  
! Bikram Connecting PES to MQCT May 15, 2021:
	  real*8 pi, rad2deg, x
	  real*8 theta,phi,thetap,phip
	  real*8 thx,phx,thxp,phxp,energy
	  
      pi=dacos(-1d0)
      rad2deg=180d0/pi
! call first to initialize
      call nd3h2pot(5d0, 1d0, 1d0, 1d0, 1d0, x)
	  
! convert angles to radians
!      alpha1 = alpha!0d0
!      beta1 = beta!/rad2deg
!      gamma1 = gamma!/rad2deg
!      beta2 = bbeta!/rad2deg
!      alpha2 = aalpha!/rad2deg
!  The Euler angles used here differ from those used in the potential
!  Coordinate transformation
      call convert(alpha,beta,gamma,aalpha,bbeta,
     &              theta,phi,thetap,phip)
! Angles defined in Maret et al., MNRAS, 399, 425–431 (2009)
! convert to degrees
      thx = theta*rad2deg
      phx = phi*rad2deg
      thxp = thetap*rad2deg
      phxp = phip*rad2deg
! call with angles in degrees
      call nd3h2pot(R, thx, phx, thxp, phxp, energy)
!      call nd3h2pot(R, thx, 0.d0, 0.d0, 0.d0, energy)
	  V = energy
! Bikram End.
	  
      END SUBROUTINE USER_DEFINED_PES
	  
	  include "nd3h2.f"
	  include "convert.f"
	  include "cubspl.f"
	  include "t.f"
