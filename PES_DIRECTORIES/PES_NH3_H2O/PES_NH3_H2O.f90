module PES_H2O_NH3
  implicit none
  real*8 water_cords(3,3), water_cords_norm(3,3)
  real*8 ammon_cords(3,4), ammon_cords_norm(3,4), dm1, dm2
  real*8, parameter :: bohr2ang = 0.52917721067d0
  real*8, parameter :: water_masses(3) = [ 1.00782503210d0, 15.99491461956d0, 1.00782503210d0]
  real*8, parameter :: ammon_masses(4) = [14.00307400520d0,  1.00782503210d0, 1.00782503210d0, 1.00782503210d0]
  real*8, parameter :: k_Boltz = 1.38065050d-23, J2eV = 6.241509470d18, eV2wn = 8065.5448d0
  real*8, parameter :: cm2bohr = 1.0d8/bohr2ang
  real*8, parameter :: eps2eps = k_Boltz*J2eV*eV2wn
  real*8, parameter :: sigma_H = 2.870d0/bohr2ang, eps_H =  29.20d0*eps2eps
  real*8, parameter :: sigma_N = 3.700d0/bohr2ang, eps_N =  95.00d0*eps2eps
  real*8, parameter :: sigma_O = 3.460d0/bohr2ang, eps_O = 118.00d0*eps2eps
  real*8, parameter :: n_atom1 = 3, n_atom2 = 4

contains
  ! This is the actual subroutine to generate PES of NH3-H2O
  subroutine PES(V, R, eu_ang1, eu_ang2)
    use isotopologue_coordinates_conversion_mod
    implicit none
    integer i, j, k
    real*8 V, R, eu_ang1(3), eu_ang2(3)
    real*8 water_euler_rot(3,3), water_cords_rot(3,3)
    real*8 ammon_euler_rot(3,3), ammon_cords_rot(3,4)
    real*16 sig1, sig2, eps1, eps2, sig, eps, calc_R, calc_V
    logical :: first = .false.
    
    ! Initilization
    if(.not.first) then
    
      !water coordinates as H-O-H
      water_cords(1,1) = -0.774880d0
      water_cords(2,1) =  0.0d0
      water_cords(3,1) = -0.53470d0
      water_cords(1,2) =  0.0d0
      water_cords(2,2) =  0.0d0
      water_cords(3,2) =  0.06530d0
      water_cords(1,3) =  0.774880d0
      water_cords(2,3) =  0.0d0
      water_cords(3,3) = -0.53470d0
      
      !ammonia coordinates as N-H-H-H
      ! ammon_cords(1,1) =  0.0d0
      ! ammon_cords(2,1) =  0.0d0
      ! ammon_cords(3,1) =  0.1127504580141590d0
      ! ammon_cords(1,2) =  1.7959038575956300d0
      ! ammon_cords(2,2) =  0.0d0
      ! ammon_cords(3,2) = -0.5221981188548880d0
      ! ammon_cords(1,3) = -0.8979532108530360d0
      ! ammon_cords(2,3) =  1.5552991036272100d0
      ! ammon_cords(3,3) = -0.5221981188548880d0
      ! ammon_cords(1,4) = -0.8979532108530360d0
      ! ammon_cords(2,4) = -1.5552991036272100d0
      ! ammon_cords(3,4) = -0.5221981188548880d0
	  ammon_cords(1,1) =  0.0d0
      ammon_cords(2,1) =  0.0d0
      ammon_cords(3,1) =  0.198900d0
      ammon_cords(1,2) =  1.772061181223720d0
      ammon_cords(2,2) =  0.0d0
      ammon_cords(3,2) = -0.522198d0
      ammon_cords(1,3) = -0.8860305906118590d0
      ammon_cords(2,3) =  1.534650d0
      ammon_cords(3,3) = -0.522198d0
      ammon_cords(1,4) = -0.8860305906118590d0
      ammon_cords(2,4) = -1.534650d0
      ammon_cords(3,4) = -0.522198d0
      
      water_cords_norm = normalize_cartesian_matrix(water_cords, water_masses, [2,1])
      ammon_cords_norm = normalize_cartesian_matrix(ammon_cords, ammon_masses, [1,2])
      
      print*, "Water normalized-cartesian-coordinates"
      do i = 1, 3
        write(*,*) water_cords_norm(1,i), water_cords_norm(2,i), water_cords_norm(3,i)
      end do
      print*,""
      print*, "Ammonia normalized-cartesian-coordinates"
      do i = 1, 4
        write(*,*) ammon_cords_norm(1,i), ammon_cords_norm(2,i), ammon_cords_norm(3,i)
      end do
      print*,""
      print*, "Water atom-atom bond distances (Ang)"
      do i = 1, 2
        do j = i+1, 3
          write(*,'(a,i0,a,i0,2x,e19.12)') "#atom_",i, "_",j, sqrt((water_cords_norm(1,i) - water_cords_norm(1,j))**2 + (water_cords_norm(2,i) - water_cords_norm(2,j))**2 + (water_cords_norm(3,i) - water_cords_norm(3,j))**2)*bohr2ang
        end do
      end do
      print*,""
      print*, "Ammonia atom-atom bond distances (Ang)"
      do i = 1, 3
        do j = i+1, 4
          write(*,'(a,i0,a,i0,2x,e19.12)') "#atom_",i, "_",j, sqrt((ammon_cords_norm(1,i) - ammon_cords_norm(1,j))**2 + (ammon_cords_norm(2,i) - ammon_cords_norm(2,j))**2 + (ammon_cords_norm(3,i) - ammon_cords_norm(3,j))**2)*bohr2ang
        end do
      end do
      print*, "Initialization of PES Done."
      first = .true.
    end if
    
    ! performing Euler angles rotation (intrinsic) to obtain new rotation matrix, water and ammonia, repsectively
    water_euler_rot = euler_to_rotation_matrix(eu_ang1)
    ammon_euler_rot = euler_to_rotation_matrix(eu_ang2)
    
    ! getting new coordinates after Euler rotation
    water_cords_rot = matmul(water_euler_rot, water_cords_norm)
    ammon_cords_rot = matmul(ammon_euler_rot, ammon_cords_norm)
    
    ! shiffting ammonia to distance of R (Bohr) along Z-axis
    ammon_cords_rot(3,:) = ammon_cords_rot(3,:) + R
    
    if(.not.first) then
      print*, "Water normalized-cartesian-coordinates"
      do i = 1, 3
        write(*,*) water_cords_rot(1,i), water_cords_rot(2,i), water_cords_rot(3,i)
      end do
      print*,""
      print*, "Ammonia normalized-cartesian-coordinates"
      do i = 1, 4
        write(*,*) ammon_cords_rot(1,i), ammon_cords_rot(2,i), ammon_cords_rot(3,i)
      end do
      print*,""
      print*, "Water atom-atom bond distances (Ang)"
      do i = 1, 2
        do j = i+1, 3
          write(*,'(a,i0,a,i0,2x,e19.12)') "#atom_",i, "_",j, sqrt((water_cords_rot(1,i) - water_cords_rot(1,j))**2 + (water_cords_rot(2,i) - water_cords_rot(2,j))**2 + (water_cords_rot(3,i) - water_cords_rot(3,j))**2)*bohr2ang
        end do
      end do
      print*,""
      print*, "Ammonia atom-atom bond distances (Ang)"
      do i = 1, 3
        do j = i+1, 4
          write(*,'(a,i0,a,i0,2x,e19.12)') "#atom_",i, "_",j, sqrt((ammon_cords_rot(1,i) - ammon_cords_rot(1,j))**2 + (ammon_cords_rot(2,i) - ammon_cords_rot(2,j))**2 + (ammon_cords_rot(3,i) - ammon_cords_rot(3,j))**2)*bohr2ang
        end do
      end do
    end if
        
    ! calculating sigma, eps, and R for L-J potential of only intermolecular atoms
    V = 0.0d0
    do i = 1, n_atom1
      do j = 1, n_atom2
        if(i == 2) then
          sig1 = sigma_O
          eps1 = eps_O
        else
          sig1 = sigma_H
          eps1 = eps_H
        end if
        if(j == 1) then
          sig2 = sigma_N
          eps2 = eps_N
        else
          sig2 = sigma_H
          eps2 = eps_H
        end if
        sig = (sig1 + sig2)/2d0
        eps = sqrt(eps1 * eps2)
        calc_R = get_norm(water_cords_rot(:,i) - ammon_cords_rot(:,j))
        calc_V = 4d0*eps*( (sig/calc_R)**12d0 - (sig/calc_R)**6d0 )
        V = V + calc_V
      end do
    end do
  end subroutine
end module

! program test_PES
  ! use PES_H2O_NH3
  ! implicit none
  ! integer i
  ! real*8 V, R, ang1(3), ang2(3)
  ! real*8, parameter :: pi = dacos(-1.0d0)
  
  ! do i = 1, 1
    ! R = 7.1150d0 !5.0d0 + (i-1)*0.0450d0
    ! ! generating grid of angles in Degree
	! ang1(1) = 0.0d0
	! ang1(2) = 0.0d0 + (i-1)*1.0d0
	! ang1(3) = 0.0d0
    ! ang2(1) = 0.0d0
    ! ang2(2) = 0.0d0
    ! ang2(3) = 0.0d0
	! ! converting angles to radian
	! ang1 = ang1*pi/180d0
	! ang2 = ang2*pi/180d0
    ! call PES(V, R, ang1, ang2)
	! print*, ang1(1)*180d0/pi, V
  ! end do

! end program
