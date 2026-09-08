module healpix_types
  implicit none

  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: i8b = selected_int_kind(18)
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0, 1.0d0))
  integer, parameter :: lgt = kind(.true.)
  real(dp), parameter :: pi = acos(-1.0_dp)
end module healpix_types
