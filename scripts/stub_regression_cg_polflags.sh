#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workdir="/tmp/opencode/dang_stub_regression_cg_polflags"
mkdir -p "$workdir"

cat > "$workdir/healpix_types.f90" <<'EOF'
module healpix_types
  implicit none
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: lgt = kind(.true.)
  real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
end module healpix_types
EOF

cat > "$workdir/dang_util_mod.f90" <<'EOF'
module dang_util_mod
  use healpix_types
  implicit none
  real(dp) :: t1, t2, t3, t4, t5, t6
  real(dp) :: missval = -1.6375d30
  integer(i4b) :: nbands = 1, npix = 2, nmaps = 3
  integer(i4b) :: ncomp = 1, ncg_groups = 1, iter = 2
  integer(i4b) :: verbosity = 0, rank = 0, master = 0
contains
  function mpi_wtime() result(t)
    real(dp) :: t
    t = 0.0_dp
  end function mpi_wtime

  function rand_normal(mean, stdev) result(c)
    real(dp), intent(in) :: mean, stdev
    real(dp) :: c
    c = mean
    if (stdev < 0.0_dp) c = mean
  end function rand_normal

  function count_delimits(string, delimiter) result(ndel)
    character(len=*), intent(in) :: string, delimiter
    integer(i4b) :: ndel, i
    ndel = 0
    do i = 1, len_trim(string)
      if (string(i:i) == delimiter(1:1)) ndel = ndel + 1
    end do
  end function count_delimits

  subroutine delimit_string(string, delimiter, list)
    character(len=*), intent(in) :: string, delimiter
    character(len=*), dimension(:), intent(out) :: list
    integer(i4b) :: i, j, k
    if (size(list) == 0) return
    list = ''
    j = 1
    k = 1
    do i = 1, len_trim(string)
      if (string(i:i) == delimiter(1:1)) then
        j = j + 1
        k = 1
      else
        list(j)(k:k) = string(i:i)
        k = k + 1
      end if
    end do
  end subroutine delimit_string

  function return_poltype_flag(string) result(flag)
    character(len=10), intent(in) :: string
    integer(i4b), allocatable :: flag(:)
    integer(i4b) :: n
    n = 1
    allocate(flag(n))
    if (trim(string) == 'T+Q+U') then
      flag(1) = 0
    else if (trim(string) == 'Q+U') then
      flag(1) = 8
    else if (trim(string) == 'Q') then
      flag(1) = 2
    else if (trim(string) == 'U') then
      flag(1) = 4
    else
      flag(1) = 1
    end if
  end function return_poltype_flag
end module dang_util_mod
EOF

cat > "$workdir/dang_param_mod.f90" <<'EOF'
module dang_param_mod
  use healpix_types
  implicit none
  type :: dang_params
    integer(i4b) :: ncomp = 1
    integer(i4b) :: ncggroup = 1
    integer(i4b), allocatable :: cg_max_iter(:)
    real(dp), allocatable :: cg_convergence(:)
    logical(lgt), allocatable :: cg_group_sample(:)
    character(len=10), allocatable :: cg_poltype(:)
    character(len=16) :: ml_mode = 'optimize'
  end type dang_params
end module dang_param_mod
EOF

cat > "$workdir/dang_bp_mod.f90" <<'EOF'
module dang_bp_mod
  implicit none
end module dang_bp_mod
EOF

cat > "$workdir/dang_component_mod.f90" <<'EOF'
module dang_component_mod
  use healpix_types
  implicit none

  type :: dang_comp
    logical(lgt) :: sample_amplitude = .true.
    integer(i4b) :: cg_group = 1
    integer(i4b) :: nfit = 1
    integer(i4b) :: nindices = 1
    character(len=16) :: type = 'powlaw'
    logical(lgt), allocatable :: corr(:)
    real(dp), allocatable :: amplitude(:,:)
    real(dp), allocatable :: indices(:,:,:)
    real(dp), allocatable :: template_amplitudes(:,:)
  contains
    procedure :: S
    procedure :: evalSignal
  end type dang_comp

  type :: comp_ptr
    class(dang_comp), pointer :: p => null()
  end type comp_ptr

  type(comp_ptr), allocatable, dimension(:) :: component_list

contains

  function S(self, band, pol, theta, pixel) result(val)
    class(dang_comp), intent(in) :: self
    integer(i4b), intent(in) :: band, pol
    real(dp), dimension(:), intent(in), optional :: theta
    integer(i4b), intent(in), optional :: pixel
    real(dp) :: val
    select case (pol)
    case (1)
      val = 1.0_dp
    case (2)
      val = 10.0_dp
    case (3)
      val = 100.0_dp
    case default
      val = -9999.0_dp
    end select
  end function S

  function evalSignal(self, band, pixel, pol) result(val)
    class(dang_comp), intent(in) :: self
    integer(i4b), intent(in) :: band, pixel, pol
    real(dp) :: val
    val = 0.0_dp
  end function evalSignal
end module dang_component_mod
EOF

cat > "$workdir/dang_data_mod.f90" <<'EOF'
module dang_data_mod
  use healpix_types
  use dang_util_mod, only: npix, nmaps, nbands
  implicit none
  type :: dang_data
    real(dp), allocatable :: sig_map(:,:,:)
    real(dp), allocatable :: rms_map(:,:,:)
    real(dp), allocatable :: masks(:,:)
    real(dp), allocatable :: gain(:)
  contains
    procedure :: update_sky_model
  end type dang_data
contains
  subroutine update_sky_model(self)
    class(dang_data), intent(inout) :: self
  end subroutine update_sky_model

  subroutine write_stats_to_term(self, iter)
    type(dang_data), intent(in) :: self
    integer(i4b), intent(in) :: iter
  end subroutine write_stats_to_term
end module dang_data_mod
EOF

cat > "$workdir/test_cg_polflags.f90" <<'EOF'
program test_cg_polflags
  use healpix_types
  use dang_util_mod
  use dang_component_mod
  use dang_data_mod
  use dang_cg_mod
  implicit none

  type(dang_cg_group) :: g
  type(dang_data) :: d
  type(dang_comp), target :: comp
  real(dp), allocatable :: b(:)

  call setup(g, d, comp)

  call run_case('T', 1_i4b, [1.0_dp, 4.0_dp], g, d, b)
  call run_case('Q', 2_i4b, [20.0_dp, 50.0_dp], g, d, b)
  call run_case('U', 4_i4b, [300.0_dp, 600.0_dp], g, d, b)
  call run_case('Q+U', 8_i4b, [20.0_dp, 50.0_dp, 300.0_dp, 600.0_dp], g, d, b)
  call run_case('T+Q+U', 0_i4b, [1.0_dp, 4.0_dp, 20.0_dp, 50.0_dp, 300.0_dp, 600.0_dp], g, d, b)

  write(*,'(a)') 'All pol_flag regression checks passed for src/dang_cg_mod.f90'

contains

  subroutine setup(group, data, c)
    type(dang_cg_group), intent(out) :: group
    type(dang_data), intent(out) :: data
    type(dang_comp), target, intent(out) :: c

    npix = 2
    nbands = 1
    nmaps = 3
    ncomp = 1
    iter = 2

    allocate(component_list(1))
    component_list(1)%p => c

    group%ncg_components = 1
    group%cg_group = 1
    group%nflag = 1
    group%ntemp = 0
    allocate(group%pol_flag(1))
    allocate(group%cg_component(1))
    group%cg_component(1)%p => c

    c%sample_amplitude = .true.
    c%cg_group = 1
    c%type = 'powlaw'
    c%nfit = 1
    c%nindices = 1
    allocate(c%corr(1))
    c%corr = .true.
    allocate(c%indices(0:npix-1,3,1))
    c%indices = 0.0_dp

    allocate(data%sig_map(0:npix-1,3,1))
    allocate(data%rms_map(0:npix-1,3,1))
    allocate(data%masks(0:npix-1,1))
    allocate(data%gain(1))

    data%sig_map(0,1,1) = 1.0_dp
    data%sig_map(0,2,1) = 2.0_dp
    data%sig_map(0,3,1) = 3.0_dp
    data%sig_map(1,1,1) = 4.0_dp
    data%sig_map(1,2,1) = 5.0_dp
    data%sig_map(1,3,1) = 6.0_dp

    data%rms_map = 1.0_dp
    data%masks = 1.0_dp
    data%gain = 1.0_dp
  end subroutine setup

  subroutine run_case(name, flag, expected, group, data, outvec)
    character(len=*), intent(in) :: name
    integer(i4b), intent(in) :: flag
    real(dp), dimension(:), intent(in) :: expected
    type(dang_cg_group), intent(inout) :: group
    type(dang_data), intent(in) :: data
    real(dp), allocatable, dimension(:), intent(inout) :: outvec

    group%pol_flag(1) = flag
    call group%compute_rhs(data, outvec, 1)
    call assert_close(name, outvec, expected)
    if (allocated(outvec)) deallocate(outvec)
  end subroutine run_case

  subroutine assert_close(name, got, expected)
    character(len=*), intent(in) :: name
    real(dp), dimension(:), intent(in) :: got, expected
    integer(i4b) :: i

    if (size(got) /= size(expected)) then
      write(*,'(a,a,a,i0,a,i0)') 'Case ', trim(name), ': size mismatch got=', size(got), ' expected=', size(expected)
      error stop 1
    end if

    do i = 1, size(got)
      if (abs(got(i) - expected(i)) > 1.0e-12_dp) then
        write(*,'(a,a,a,i0,a,f12.6,a,f12.6)') 'Case ', trim(name), ': mismatch at i=', i, ' got=', got(i), ' expected=', expected(i)
        error stop 1
      end if
    end do

    write(*,'(a,a,a)') 'Case ', trim(name), ' passed'
  end subroutine assert_close

end program test_cg_polflags
EOF

gfortran -c -J"$workdir" -I"$workdir" "$workdir/healpix_types.f90" -o "$workdir/healpix_types.o"
gfortran -c -J"$workdir" -I"$workdir" "$workdir/dang_util_mod.f90" -o "$workdir/dang_util_mod.o"
gfortran -c -J"$workdir" -I"$workdir" "$workdir/dang_param_mod.f90" -o "$workdir/dang_param_mod.o"
gfortran -c -J"$workdir" -I"$workdir" "$workdir/dang_bp_mod.f90" -o "$workdir/dang_bp_mod.o"
gfortran -c -J"$workdir" -I"$workdir" "$workdir/dang_component_mod.f90" -o "$workdir/dang_component_mod.o"
gfortran -c -J"$workdir" -I"$workdir" "$workdir/dang_data_mod.f90" -o "$workdir/dang_data_mod.o"
gfortran -c -J"$workdir" -I"$workdir" "$repo_root/src/dang_cg_mod.f90" -o "$workdir/dang_cg_mod.o"

gfortran -J"$workdir" -I"$workdir" -o "$workdir/test_cg_polflags.exe" \
  "$workdir/healpix_types.o" \
  "$workdir/dang_util_mod.o" \
  "$workdir/dang_param_mod.o" \
  "$workdir/dang_bp_mod.o" \
  "$workdir/dang_component_mod.o" \
  "$workdir/dang_data_mod.o" \
  "$workdir/dang_cg_mod.o" \
  "$workdir/test_cg_polflags.f90"

"$workdir/test_cg_polflags.exe"
