      SUBROUTINE EXPAND_PES(n_R, R, n_exp_terms, exp_terms)

!------------------------------------------------------------------------------------------------------------!
! This subroutine is to create the PES_EXPAN_TERMS.DAT file, with appropriate formating, 
! from expansion coefficents computed externally. The PES_EXPAN_TERMS.DAT file then 
! can be used with the MQCT code to compute state-to-state transition matrix elements.
!
! INPUTS: n_R, R, n_exp_terms, exp_terms
! n_R - the number of grid points along the molecule-molecule distance R
! R - 1D array, the values of molecule-molecule distance R (in Bohr) at grid points
! n_exp_terms - number of PES expansion terms
! exp_terms - 2D array (n_R x n_exp_terms) containing expansion coefficents (in 1/cm) at each grid point 
! 1st index is R number, and 2nd index is the expansion coefficents.
! For example, exp_terms(:,1) stores the values of coefficent for 1st expansion term at all values of R-grid
!
! OUTPUT: User ready PES_EXPAN_TERMS.DAT file.
!------------------------------------------------------------------------------------------------------------!

      implicit none
      integer i, j, n_R, n_exp_terms
      real*8 :: R(n_R)
      real*8 :: exp_terms(n_R, n_exp_terms)
      logical db_exst
      	   
      OPEN(1, FILE = "PES_EXPAN_TERMS.DAT", ACTION = "WRITE")                                    ! Creating the file
      write(1, '(a)') "Expansion terms calculated externally"                                    ! Header line
      
      DO j = 1, n_R                                                                              ! Loop over the grid-points
      DO i = 1, n_exp_terms                                                                      ! Loop over the expansion terms
      IF(i.eq.1) THEN	  
      WRITE(1,'(f6.2,1x,e19.12,3x)',ADVANCE="NO") R(j),exp_terms(j,i)                            ! Printing 1st and 2nd columns
      ELSE
      WRITE(1,'(e19.12,3x)',ADVANCE="NO") exp_terms(j,i)                                         ! Printing the rest
      ENDIF	  
      ENDDO
      WRITE(1,*)                                                                                 ! Required to move to a new line
      ENDDO 
      close(1)
      	  
      END SUBROUTINE EXPAND_PES