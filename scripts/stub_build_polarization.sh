#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workdir="/tmp/opencode/dang_stub_build"
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
  integer(i4b) :: nbands = 1, npix = 4, nmaps = 3
  integer(i4b) :: ncomp = 1, ncg_groups = 1, iter = 1
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
    val = 1.0_dp
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

gfortran -c -J"$workdir" -I"$workdir" "$workdir/healpix_types.f90"
gfortran -c -J"$workdir" -I"$workdir" "$workdir/dang_util_mod.f90"
gfortran -c -J"$workdir" -I"$workdir" "$workdir/dang_param_mod.f90"
gfortran -c -J"$workdir" -I"$workdir" "$workdir/dang_bp_mod.f90"
gfortran -c -J"$workdir" -I"$workdir" "$workdir/dang_component_mod.f90"
gfortran -c -J"$workdir" -I"$workdir" "$workdir/dang_data_mod.f90"

gfortran -fsyntax-only -J"$workdir" -I"$workdir" "$repo_root/src/dang_cg_mod.f90"

printf 'Stub syntax check passed for src/dang_cg_mod.f90\n'
