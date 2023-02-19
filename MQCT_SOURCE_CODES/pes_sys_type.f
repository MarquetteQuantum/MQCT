!!!! SYS_TYPE 1	  !! INPUT R - distance, beta - angle, OUTPUT : V - potential 
	  SUBROUTINE V_POT_DIATOM_ATOM(V,R,beta)
      USE VARIABLES
      IMPLICIT NONE
      REAL*8 V,R,r_vib,alpha,beta,gamma,aalpha,bbeta,ggamma
      REAL*8 ra,rb,r_vib_1,r_vib_2
      EXTERNAL USER_DEFINED_PES	  
      CALL USER_DEFINED_PES (V,R,0d0,0d0,0d0,beta,0d0,0d0,0d0,0d0)
      END SUBROUTINE V_POT_DIATOM_ATOM
	  
!!!! SYS_TYPE 2 !! INPUT: r_vib - vibrational distance or distances ra and rb from each atom of the molecule to collider
      SUBROUTINE V_POT_VIB_DIATOM_ATOM(V,R,r_vib,beta)
      USE VARIABLES
      IMPLICIT NONE
      REAL*8 V,R,r_vib,alpha,beta,gamma,aalpha,bbeta,ggamma
      REAL*8 ra,rb,r_vib_1,r_vib_2
      EXTERNAL USER_DEFINED_PES	  
      IF(atom_coord_dist)	THEN
      r_vib_1=
     & r_vib*atomic_masses(1)/(atomic_masses(1)+atomic_masses(2))
      r_vib_2=r_vib
     & * atomic_masses(2)/(atomic_masses(1)+atomic_masses(2))	  
      ra = dsqrt(R**2+r_vib_1**2-2*R*r_vib_1*dcos(beta))
      rb = dsqrt(R**2+r_vib_2**2+2*R*r_vib_2*dcos(beta))	  
      CALL USER_DEFINED_PES(V,r_vib,ra,rb,0d0,0d0,0d0,0d0,0d0,0d0)	
      ELSE	  
      CALL USER_DEFINED_PES(V,R,r_vib,0d0,0d0,beta,0d0,0d0,0d0,0d0)
      ENDIF	  
      END SUBROUTINE V_POT_VIB_DIATOM_ATOM

!!!! SYS_TYPE 3,4 INPUT Euler's angles beta,gamma in the body fixed reference frame ( -beta=theta and -gamma=phi  where angles theta and phi defined in the molecule fixed referneced frame as defined by Green, JCP, 1976)  
      SUBROUTINE V_POT_TOP_ATOM(V,R,beta,gamma)
      IMPLICIT NONE
      REAL*8 V,R,beta,gamma
      EXTERNAL USER_DEFINED_PES	  
      CALL USER_DEFINED_PES(V,R,0d0,0d0,0d0,-beta,-gamma,0d0,0d0,0d0)	  
      END SUBROUTINE V_POT_TOP_ATOM
	  
!!!! SYS_TYPE 5	  
      SUBROUTINE V_POT_DIAT_DIAT(V,R,beta,bbeta,aalpha) !!INPUT: beta,bbeta - theta 1,2 , gamma - phi 
      IMPLICIT NONE
      REAL*8 V,R,beta,bbeta,aalpha
      EXTERNAL USER_DEFINED_PES	  
      CALL USER_DEFINED_PES(V,R,0d0,0d0,0,beta,aalpha,0d0,bbeta,0d0)	  
      END SUBROUTINE V_POT_DIAT_DIAT	  
	  
!!!! SYS_TYPE 6	!! INPUT: SAME AS SYS 5, BUT r_vib1 and r_vib2 - vibrational distances
      SUBROUTINE V_POT_VIB_DIAT_DIAT(V,R,r_vib1,r_vib2,
     & beta,bbeta,gamma)
      REAL*8 V,R,r_vib1,r_vib2,beta,bbeta,gamma
      EXTERNAL USER_DEFINED_PES	  
      CALL 
     & USER_DEFINED_PES(V,R,r_vib1,r_vib2,0d0,beta,gamma,0d0,bbeta,0d0)	 
      END SUBROUTINE V_POT_VIB_DIAT_DIAT

	  
!!!! SYS_TYPE 7 	!! INPUT: Euler' angles of #1 molecule: beta,gamma, #2 - aalpha,bbeta with respect to R axis of the body-fixed refence frame 
      SUBROUTINE V_POT_SYM_DIATOM(V,R,alpha,beta,gamma,
     & aalpha,bbeta)
      REAL*8 V,R,alpha,beta,gamma,
     & aalpha,bbeta,ggamma
      EXTERNAL USER_DEFINED_PES	 
      CALL USER_DEFINED_PES(V,R,0d0,0d0,0d0,beta,gamma,aalpha,bbeta,0) 
      END SUBROUTINE V_POT_SYM_DIATOM


!!!! SYS_TYPE 8		!! INPUT: Euler' angles of #1 molecule: beta,gamma, #2 - aalpha,bbeta with respect to R axis of the body-fixed refence frame   
      SUBROUTINE V_POT_ASYM_DIATOM(V,R,alpha,beta,gamma,
     & aalpha,bbeta)
      IMPLICIT NONE	 
      REAL*8 V,R,alpha,beta,gamma,
     & aalpha,bbeta,ggamma
      EXTERNAL USER_DEFINED_PES	 
      CALL USER_DEFINED_PES(V,R,0d0,0d0,0d0,beta,gamma,aalpha,bbeta,0)	 
      END SUBROUTINE V_POT_ASYM_DIATOM	  

!!!! SYS_TYPE 9	 !!!! alpha,beta,gamme - Eulers' angles #1 molecule, aalpha,bbeta,ggama for #2  
      SUBROUTINE V_POT_ASYM_SYM(V,R,alpha,beta,gamma,
     & aalpha,bbeta,ggamma)
      IMPLICIT NONE	 
      REAL*8 V,R,alpha,beta,gamma,
     & aalpha,bbeta,ggamma
      EXTERNAL USER_DEFINED_PES	 
      CALL
     & USER_DEFINED_PES(V,R,0,0,alpha,beta,gamma,aalpha,bbeta,ggamma)
      END SUBROUTINE V_POT_ASYM_SYM	 


!!!! SYS_TYPE 0		!!!! alpha,beta,gamme - Eulers' angles #1 molecule, aalpha,bbeta,ggama for #2  
      SUBROUTINE V_POT_ASYM_ASYM(V,R,alpha,beta,gamma,
     & aalpha,bbeta,ggamma)
      IMPLICIT NONE	 
      REAL*8 V,R,alpha,beta,gamma,
     & aalpha,bbeta,ggamma
      EXTERNAL USER_DEFINED_PES	 
      CALL
     & USER_DEFINED_PES(V,R,0,0,alpha,beta,gamma,aalpha,bbeta,ggamma)
      END SUBROUTINE V_POT_ASYM_ASYM	  