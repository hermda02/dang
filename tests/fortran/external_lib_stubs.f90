module pix_tools
  implicit none
end module pix_tools

module udgrade_nr
  use healpix_types
  implicit none
contains
  subroutine udgrade_ring(data_in, nside_in, data_out, nside_out)
    real(dp), dimension(:,:), intent(in) :: data_in
    integer(i4b), intent(in) :: nside_in, nside_out
    real(dp), dimension(:,:), intent(out) :: data_out
    integer(i4b) :: n1, n2

    n1 = min(size(data_in, 1), size(data_out, 1))
    n2 = min(size(data_in, 2), size(data_out, 2))
    data_out = 0.0_dp
    data_out(1:n1, 1:n2) = data_in(1:n1, 1:n2)
  end subroutine udgrade_ring
end module udgrade_nr

module head_fits
  use healpix_types
  implicit none
  interface add_card
    module procedure add_card_char
    module procedure add_card_int
  end interface add_card
contains
  subroutine add_card_char(line, key, value, comment)
    character(len=*), dimension(:), intent(inout) :: line
    character(len=*), intent(in) :: key, value, comment
    line(1) = key
  end subroutine add_card_char

  subroutine add_card_int(line, key, value, comment)
    character(len=*), dimension(:), intent(inout) :: line
    character(len=*), intent(in) :: key, comment
    integer(i4b), intent(in) :: value
    line(1) = key
  end subroutine add_card_int
end module head_fits

module fitstools
  use healpix_types
  implicit none
contains
  subroutine write_bintab(map, npix, nmaps, header, nlheader, outfile)
    real(dp), dimension(:,:), intent(in) :: map
    integer(i4b), intent(in) :: npix, nmaps, nlheader
    character(len=80), dimension(:), intent(in) :: header
    character(len=*), intent(in) :: outfile
  end subroutine write_bintab

  subroutine read_bintab(filename, map, npix, nmaps, nullval, anynull, header)
    character(len=*), intent(in) :: filename
    real(dp), dimension(:,:), intent(inout) :: map
    integer(i4b), intent(in) :: npix, nmaps
    real(dp), intent(out) :: nullval
    logical(lgt), intent(out) :: anynull
    character(len=80), dimension(:), intent(inout), optional :: header

    nullval = 0.0_dp
    anynull = .false.
  end subroutine read_bintab
end module fitstools

module mpi
  implicit none
  integer, parameter :: mpi_status_size = 5
  integer, parameter :: MPI_COMM_WORLD = 0
contains
  subroutine mpi_init(ierr)
    integer, intent(out) :: ierr
    ierr = 0
  end subroutine mpi_init

  subroutine mpi_comm_rank(comm, rank, ierr)
    integer, intent(in) :: comm
    integer, intent(out) :: rank, ierr
    rank = 0
    ierr = 0
  end subroutine mpi_comm_rank

  subroutine mpi_comm_size(comm, nprocs, ierr)
    integer, intent(in) :: comm
    integer, intent(out) :: nprocs, ierr
    nprocs = 1
    ierr = 0
  end subroutine mpi_comm_size

  real(8) function mpi_wtime()
    mpi_wtime = 0.0d0
  end function mpi_wtime
end module mpi
