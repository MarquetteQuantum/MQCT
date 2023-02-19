      SUBROUTINE DEFINE_BASIS(n_ch, v1, j1, v2, j2, E,
     & n_r_vib1, r_vib1, wg1, wavefunction1, 
     & n_r_vib2, r_vib2, wg2, wavefunction2)
	
!-------------------------------------------------------------------------------------------------------------! 
! This subroutine can be used to generate USER_DEFINED_BASIS.DAT file for SYS_TYPE = 2 and 6.
!
! INPUTS:  n_ch, v1, j1, v2, j2, E, n_r_vib1, r_vib1, wg1, wavefunction1, n_r_vib2, r_vib2, wg2, wavefunction2     
! n_ch - number of channels in the basis     
! v1 - 1D array, vibrational quantum numbers for each channel of the 1st diatom 
! j1 - 1D array, rotational  quantum numbers for each channel of the 1st diatom
! v2 - 1D array, vibrational quantum numbers for each channel of the 2nd diatom 
! j2 - 1D array, rotational  quantum numbers for each channel of the 2nd diatom 
! E  - 1D array, energies of rotational-vibrational channels in (1/cm)
!
! MOLECULE#1:
! n_r_vib1 - the number of grid points along the vibrational coordinate of 1st diatom
! r_vib1 - 1D array, the values of vibrational coordinate (Bohr) at the grid points for the 1st diatom  
! wg1 - 1D array, the weights of the vibrational grid points (jacobian) for the 1st diatom
! wavefunction1 - 2D array (n_r_vib1 x n_ch), the values ro-vibrational wavefunctions on the grid
!
! MOLECULE#2:
! n_r_vib2 - the number of grid points along the vibrational coordinate of 2nd diatom
! r_vib2 - 1D array, the values of vibrational coordinate (Bohr) at the grid points for the 2nd diatom 
! wg2 - 1D array, the weights of the vibrational grid points (jacobian) for the 2nd diatom
! wavefunction2 - 2D array (n_r_vib2 x n_ch), the values ro-vibrational wavefunctions on the grid 
!
! OUTPUT: user ready USER_DEFINED_BASIS.DAT file.
!-------------------------------------------------------------------------------------------------------------!

      implicit none
      integer i, k, n_r_vib1, n_r_vib2, n_ch
      integer :: v1(n_ch), j1(n_ch), v2(n_ch), j2(n_ch)
      real*8 :: E(n_ch)
      real*8 :: r_vib1(n_r_vib1), wg1(n_r_vib1)
      real*8 :: r_vib2(n_r_vib2), wg2(n_r_vib2)
      real*8 :: wavefunction1(n_r_vib1, n_ch)
      real*8 :: wavefunction2(n_r_vib2, n_ch)
	  
      OPEN(1, FILE="USER_DEFINED_BASIS.DAT",ACTION="WRITE")                  ! Creating the file

      if(n_r_vib2.eq.0) then                                                 ! For SYS_TYPE = 2
      write(1, '(9x,a3,3(11x,a1))') '#N=', 'J', 'V', 'E'                     ! Header for the list of channels
      DO k = 1, n_ch                                                         ! Loop over the channels
      write(1, *) k,  j1(k), v1(k), E(k)                                     ! Channel number, ro-vib quantum numbers and energy
      ENDDO                                                                  
      else                                                                   ! For SYS_TYPE = 6
      write(1, '(9x,a3,5(11x,a2))') '#N=', 'J1', 'V1', 'J2', 'V2', 'E'       ! Header for the list of channels
      DO k = 1, n_ch                                                         ! Loop over the channels
      write(1, *) k, j1(k), v1(k), j2(k), v2(k), E(k)                        ! Channel number,  ro-vib quantum numbers and energy
      ENDDO                                                                  
      endif                                                                  

!     PRINTING THE GRID	                                                     
      WRITE(1, *) "r_vib1(i),          J1(i)"                                ! A header
      DO i=1, n_r_vib1                                                       ! Loop over the vibrational grid points	  
      WRITE(1, *) r_vib1(i), wg1(i)                                          ! Vibrational distances and weights for 1st diatom
      ENDDO                                                                  
      
      if(n_r_vib2.gt.0) then                                                 ! If SYS_TYPE = 6
      WRITE(1, *) "r_vib2(i),          J2(i)"                                ! A header
      DO i=1, n_r_vib2                                                       ! Loop over the vibrational grid points	  
      WRITE(1, *) r_vib2(i), wg2(i)                                          ! Vibrational distances and weights for 2nd diatom
      ENDDO                                                                  
      endif                                                                  
      
      WRITE(1, *) "WAVEFUNCTIONS ARE LISTED BELOW"	                         
      DO k=1, n_ch                                                           ! Loop over the channels
!     PRINTING WAVEFUNCTIONS	                                             
      WRITE(1, *) "MOLECULE#1, CHANNEL=", k                                  ! A header
      DO i=1, n_r_vib1	                                                     
      WRITE(1, *) wavefunction1(i, k)                                        ! The value of wavefunction for k-channel at i-point of the grid for 1st molecule
      ENDDO                                                                  
      WRITE(1, *)                                                            ! Leave it blank	  
      ENDDO                                                                  

      if(n_r_vib2.gt.0) then                                                 ! If SYS_TYPE = 6
      DO k=1, n_ch                                                           ! Loop over the channels
!     PRINTING WAVEFUNCTIONS	                                             
      WRITE(1, *) "MOLECULE#2, CHANNEL=", k                                  ! A header
      DO i=1, n_r_vib2	                                                     
      WRITE(1, *) wavefunction2(i, k)                                        ! The value of wavefunction for k-channel at i-point of the grid for 2nd molecule
      ENDDO                                                                  
      WRITE(1, *)                                                            ! LEAVE IT BLANK	  
      ENDDO
      endif
	  
      CLOSE(1)
      
      end SUBROUTINE