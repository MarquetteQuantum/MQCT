module isotopologue_coordinates_conversion_mod
  use iso_fortran_env, only: real64
  use lapack_interface_mod
  implicit none

  ! Any number less than this is considered equal to 0
  real(real64), private :: zero_tolerance = 1d-10
  real(real64), private :: pi = acos(-1d0)

  private :: operator (-)
  interface operator (-)
    module procedure :: matrix_minus_vector
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Subtracts a given vector from all columns of a given matrix.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function matrix_minus_vector(matrix, vector) result(res)
    real(real64), intent(in) :: matrix(:, :)
    real(real64), intent(in) :: vector(:)
    real(real64), allocatable :: res(:, :)
    integer :: j

    res = matrix
    do j = 1, size(res, 2)
      res(:, j) = res(:, j) - vector
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates center of mass for a given matrix of cartesian coordinates (3xN, order: x, y, z) and corresponding particle masses (N).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mass_center(coords_cartesian, masses) result(mass_center)
    real(real64), intent(in) :: coords_cartesian(:, :)
    real(real64), intent(in) :: masses(:)
    real(real64), allocatable :: mass_center(:)
    real(real64), allocatable :: weights(:)

    weights = masses / sum(masses)
    mass_center = matmul(coords_cartesian, weights)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates tensor of inertia.
! coords_cartesian: matrix of cartesian coordinates, 3xN, each point is written in column in order: x, y, z.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_inertia_tensor(coords_cartesian, masses) result(inertia_tensor)
    real(real64), intent(in) :: coords_cartesian(:, :)
    real(real64), intent(in) :: masses(:)
    real(real64) :: inertia_tensor(3, 3)

    inertia_tensor(1, 1) = dot_product(masses, coords_cartesian(2, :) ** 2 + coords_cartesian(3, :) ** 2)
    inertia_tensor(2, 2) = dot_product(masses, coords_cartesian(1, :) ** 2 + coords_cartesian(3, :) ** 2)
    inertia_tensor(3, 3) = dot_product(masses, coords_cartesian(1, :) ** 2 + coords_cartesian(2, :) ** 2)
    inertia_tensor(1, 2) = -dot_product(masses, coords_cartesian(1, :) * coords_cartesian(2, :))
    inertia_tensor(1, 3) = -dot_product(masses, coords_cartesian(1, :) * coords_cartesian(3, :))
    inertia_tensor(2, 3) = -dot_product(masses, coords_cartesian(2, :) * coords_cartesian(3, :))
    inertia_tensor(2, 1) = inertia_tensor(1, 2)
    inertia_tensor(3, 1) = inertia_tensor(1, 3)
    inertia_tensor(3, 2) = inertia_tensor(2, 3)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes vector L2-norm.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_norm(v) result(norm)
    real(real64), intent(in) :: v(3)
    real(real64) :: norm
    norm = sqrt(sum(v ** 2))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes normalized cross product between 3D vectors v1 and v2.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function cross_product(v1, v2) result(v3)
    real(real64), intent(in) :: v1(3), v2(3)
    real(real64) :: v3(3)
    
    v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
    v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
    v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
    v3 = v3 / get_norm(v3)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Rearranges order and direction of principal axes of inertia to right-handed zxy according to the incides of atoms that have to have positive coordinates along z- and x-axes.
! If there is no atom in the direction of an axis, specify 0.
! Assumes zx_atom_index(1) != 0.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function rearrange_inertia_tensor(inertia_tensor, cartesian, zx_atom_index) result(principal_axes)
    real(real64), intent(in) :: inertia_tensor(3, 3)
    real(real64), intent(in) :: cartesian(:, :)
    integer, intent(in) :: zx_atom_index(2)
    real(real64) :: principal_axes(3, 3)
    integer :: axes_to_find, i, j
    integer :: max_axis_ind(2)
    real(real64) :: axis_projection(3)
    integer, parameter :: axes_inds(2) = [3, 1]
    integer, parameter :: dependent_axis_ind = 2

    if (zx_atom_index(2) == 0) then
      axes_to_find = 1
    else
      axes_to_find = 2
    end if

    principal_axes = 0
    do i = 1, axes_to_find
      do j = 1, 3
        axis_projection(j) = dot_product(cartesian(:, zx_atom_index(i)), inertia_tensor(:, j))
      end do

      max_axis_ind(i) = maxloc(abs(axis_projection), 1)
      principal_axes(:, axes_inds(i)) = inertia_tensor(:, max_axis_ind(i))
      if (axis_projection(max_axis_ind(i)) < 0) then
        principal_axes(:, axes_inds(i)) = -principal_axes(:, axes_inds(i))
      end if
    end do

    if (axes_to_find == 1) then
      principal_axes(dependent_axis_ind, dependent_axis_ind) = 1
      principal_axes(:, axes_inds(2)) = cross_product(principal_axes(:, dependent_axis_ind), principal_axes(:, axes_inds(1)))
    else
      if (max_axis_ind(1) == max_axis_ind(2)) then
        stop 'Error: the same axis is selected during labeling'
      end if
      principal_axes(:, dependent_axis_ind) = cross_product(principal_axes(:, axes_inds(1)), principal_axes(:, axes_inds(2)))
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns principal axes of inertia.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_principal_axes(cartesian, masses, zx_atom_index) result(principal_axes)
    real(real64), intent(in) :: cartesian(:, :)
    real(real64), intent(in) :: masses(:)
    integer, intent(in) :: zx_atom_index(2)
    real(real64) :: principal_axes(3, 3)
    real(real64) :: inertia_moments(3)
    real(real64) :: inertia_tensor(3, 3)

    if (zx_atom_index(1) == 0) then
      principal_axes = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
      return
    end if

    inertia_tensor = get_inertia_tensor(cartesian, masses)
    call lapack_eigensolver_real_plain(inertia_tensor, inertia_moments)
    principal_axes = rearrange_inertia_tensor(inertia_tensor, cartesian, zx_atom_index)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Shifts and rotates atoms described by a given cartesian matrix to make the center of mass coincide with the origin and the principal axes of inertia with the cartesian axes. 
! Axes labeling is according to the specified atom indices, see `rearrange_inertia_tensor`.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function normalize_cartesian_matrix(cartesian, masses, zx_atom_index) result(cartesian_normalized)
    real(real64), intent(in) :: cartesian(:, :) ! 3 x N
    real(real64), intent(in) :: masses(:)
    integer, intent(in) :: zx_atom_index(2)
    real(real64) :: cartesian_normalized(size(cartesian, 1), size(cartesian, 2))
    real(real64) :: mass_center(3), inertia_moments(3)
    real(real64) :: inertia_tensor(3, 3), principal_axes(3, 3)
    real(real64), allocatable :: cartesian_shifted(:, :)
  
    mass_center = get_mass_center(cartesian, masses)
    cartesian_shifted = cartesian - mass_center
    principal_axes = get_principal_axes(cartesian_shifted, masses, zx_atom_index)
    cartesian_normalized = matmul(transpose(principal_axes), cartesian_shifted)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns cartesian rotation matrix about a given unit axis specified by (x, y, z) by angle a (in rad).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_rotation_matrix(axis, a) result(rotation_matrix)
    real(real64), intent(in) :: axis(3)
    real(real64), intent(in) :: a
    real(real64) :: rotation_matrix(3, 3)
    real(real64) :: x, y, z

    x = axis(1)
    y = axis(2)
    z = axis(3)
    rotation_matrix(1, 1) = cos(a) + x**2 * (1 - cos(a))
    rotation_matrix(2, 1) = y * x * (1 - cos(a)) + z * sin(a)
    rotation_matrix(3, 1) = z * x * (1 - cos(a)) - y * sin(a)
    rotation_matrix(1, 2) = x * y * (1 - cos(a)) - z * sin(a)
    rotation_matrix(2, 2) = cos(a) + y**2 * (1 - cos(a))
    rotation_matrix(3, 2) = z * y * (1 - cos(a)) + x * sin(a)
    rotation_matrix(1, 3) = x * z * (1 - cos(a)) + y * sin(a)
    rotation_matrix(2, 3) = y * z * (1 - cos(a)) - x * sin(a)
    rotation_matrix(3, 3) = cos(a) + z**2 * (1 - cos(a))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Transforms intrinsic Euler angles (z-y-z) to the cartesian rotation matrix.
! The Euler angles are given in order: 1st, 2nd, 3rd.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function euler_to_rotation_matrix(euler) result(rotation_matrix)
    real(real64), intent(in) :: euler(3)
    real(real64) :: rotation_matrix(3, 3)
    integer :: i
    integer, parameter :: rot_inds(3) = [3, 2, 3]
    real(real64) :: next_rotation(3, 3)

    rotation_matrix = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
    do i = 1, size(rot_inds)
      next_rotation = get_rotation_matrix(rotation_matrix(:, rot_inds(i)), euler(i))
      rotation_matrix = matmul(next_rotation, rotation_matrix)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Floating-point aware version of acos.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function tolerant_acos(arg) result(res)
    real(real64), intent(in) :: arg
    real(real64) :: res

    if (arg > 1) then
      if (arg - zero_tolerance < 1) then
        res = acos(1d0)
      else
        stop 'Error: arg is > 1'
      end if

    elseif (arg < -1) then
      if (arg + zero_tolerance > -1) then
        res = acos(-1d0)
      else
        stop 'Error: arg is < -1'
      end if

    else
      res = acos(arg)
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Transforms cartesian rotation matrix to intrinsic Euler angles (z-y-z).
! The Euler angles are given in order: 1st, 2nd, 3rd.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function rotation_matrix_to_euler(rotation_matrix) result(euler)
    real(real64), intent(in) :: rotation_matrix(3, 3)
    real(real64) :: euler(3)
    real(real64) :: n(3)
    
    euler(2) = tolerant_acos(rotation_matrix(3, 3))
    if (euler(2) < zero_tolerance) then
      n = rotation_matrix(:, 2)
      euler(3) = 0
    else
      n = cross_product([0d0, 0d0, 1d0], rotation_matrix(:, 3)) ! Normal vector to z z' plane
      euler(3) = tolerant_acos(dot_product(n, rotation_matrix(:, 2)))
      if (rotation_matrix(3, 2) < -zero_tolerance) then
        euler(3) = 2*pi - euler(3)
      end if
    end if

    euler(1) = tolerant_acos(n(2))
    if (n(1) > zero_tolerance) then
      euler(1) = 2*pi - euler(1)
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts coordinates given w.r.t. given isotopomers to the corresponding coordinates w.r.t. the reference isotopomers.
! cartesian_iso_1 and 2 are 3 x N normalized cartesian coordinates of the two isotopomers, i.e. center of mass is at the origin and principal axes of inertia are aligned with the cartesian axes.
! Coordinates of each atom are written in columns in order: x, y, z. Use normalize_cartesian_matrix to normalize an arbitrary cartesian matrix.
! Euler coordinates are center of mass distance between the two isotopomers (i.e. origin distance between the two cartesian systems of cartesian_iso_1 and 2) 
! and 5 Euler angles corresponding to rotations around body-fixed z-y-z axes of each body (i.e. axes initially aligned with the corresponding cartesian axes).
! The euler coordinates are given in order: R, alpha1, beta1, gamma1, beta2, gamma2. 
! The first rotation on the second body (alpha2) is neglected since only relative orientation of the two bodies matters.
! masses_ref lists the new masses on both isotopomers, with respect to which the new euler_ref angles will be calculated.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_coordinates_to_reference(cartesian_iso_1, cartesian_iso_2, euler_iso, masses_ref_1, masses_ref_2, zx_atom_index_1, zx_atom_index_2) result(euler_ref)
    real(real64), intent(in) :: cartesian_iso_1(:, :), cartesian_iso_2(:, :)
    real(real64), intent(in) :: euler_iso(6)
    real(real64), intent(in) :: masses_ref_1(:), masses_ref_2(:)
    integer, intent(in) :: zx_atom_index_1(2), zx_atom_index_2(2)
    real(real64) :: euler_ref(6)
    real(real64) :: rotation_angle_combined
    real(real64) :: mass_center_ref_1(3), mass_center_ref_2(3), mass_center_shifted_ref_2(3), rotation_axis_combined(3), euler_1(3), euler_2(3)
    real(real64) :: euler_iso_rotation_1(3, 3), euler_iso_rotation_2(3, 3), combined_rotation(3, 3), principal_axes_1(3, 3), principal_axes_2(3, 3)
    real(real64), allocatable :: cartesian_iso_rotated_1(:, :), cartesian_iso_rotated_2(:, :), cartesian_iso_shifted_rotated_2(:, :), cartesian_iso_combined(:, :), cartesian_iso_shifted_combined(:, :)
    real(real64), allocatable :: cartesian_ref_combined(:, :), cartesian_ref_1(:, :), cartesian_shifted_ref_2(:, :), cartesian_ref_2(:, :)

    euler_iso_rotation_1 = euler_to_rotation_matrix([euler_iso(2), euler_iso(3), euler_iso(4)])
    euler_iso_rotation_2 = euler_to_rotation_matrix([0d0, euler_iso(5), euler_iso(6)])
    cartesian_iso_rotated_1 = matmul(euler_iso_rotation_1, cartesian_iso_1)
    cartesian_iso_rotated_2 = matmul(euler_iso_rotation_2, cartesian_iso_2)
    cartesian_iso_shifted_rotated_2 = cartesian_iso_rotated_2
    cartesian_iso_shifted_rotated_2(3, :) = cartesian_iso_shifted_rotated_2(3, :) + euler_iso(1)

    mass_center_ref_1 = get_mass_center(cartesian_iso_rotated_1, masses_ref_1)
    mass_center_ref_2 = get_mass_center(cartesian_iso_shifted_rotated_2, masses_ref_2)
    euler_ref(1) = sqrt(sum((mass_center_ref_2 - mass_center_ref_1) ** 2))

    allocate(cartesian_iso_combined(3, size(cartesian_iso_1, 2) + size(cartesian_iso_2, 2)))
    cartesian_iso_combined(:, 1:size(cartesian_iso_1, 2)) = cartesian_iso_rotated_1
    cartesian_iso_combined(:, size(cartesian_iso_1, 2) + 1 : size(cartesian_iso_combined, 2)) = cartesian_iso_shifted_rotated_2
    cartesian_iso_shifted_combined = cartesian_iso_combined
    cartesian_iso_shifted_combined = cartesian_iso_shifted_combined - mass_center_ref_1
    mass_center_shifted_ref_2 = mass_center_ref_2 - mass_center_ref_1
    rotation_angle_combined = tolerant_acos(mass_center_shifted_ref_2(3) / euler_ref(1))

    if (rotation_angle_combined > zero_tolerance) then
      rotation_axis_combined = cross_product(mass_center_shifted_ref_2, [0d0, 0d0, euler_ref(1)])
      combined_rotation = get_rotation_matrix(rotation_axis_combined, rotation_angle_combined)
      cartesian_ref_combined = matmul(combined_rotation, cartesian_iso_shifted_combined)
    else
      cartesian_ref_combined = cartesian_iso_shifted_combined
    end if

    cartesian_ref_1 = cartesian_ref_combined(:, 1:size(cartesian_iso_1, 2))
    cartesian_shifted_ref_2 = cartesian_ref_combined(:, size(cartesian_iso_1, 2) + 1 : size(cartesian_ref_combined, 2))
    cartesian_ref_2 = cartesian_shifted_ref_2
    cartesian_ref_2(3, :) = cartesian_ref_2(3, :) - euler_ref(1)

    principal_axes_1 = get_principal_axes(cartesian_ref_1, masses_ref_1, zx_atom_index_1)
    euler_1 = rotation_matrix_to_euler(principal_axes_1)
    principal_axes_2 = get_principal_axes(cartesian_ref_2, masses_ref_2, zx_atom_index_2)
    euler_2 = rotation_matrix_to_euler(principal_axes_2)

    euler_ref(2) = abs(euler_1(1) - euler_2(1))
    euler_ref(3:4) = euler_1(2:3)
    euler_ref(5:6) = euler_2(2:3)
  end function

end module
