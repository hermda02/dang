module init_mod
    use healpix_types
    use utility_mod
    implicit none

    ! Constants
    real(dp)           :: k_B     = 1.3806503d-23
    real(dp)           :: h       = 1.0545726691251021d-34 * 2.d0 * pi
    real(dp)           :: c       = 2.99792458d8
    real(dp)           :: T_CMB   = 2.7255d0
    real(dp)           :: t1, t2, t3
    integer(i4b)       :: nbands
    integer(i4b)       :: master      = 0 
    integer(i4b)       :: from_master = 1
    integer(i4b)       :: from_worker = 2


contains 

    subroutine init_mpi()
      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
      call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
      if (rank == 0) then
         write(*,'(a,i8)') ' The number of processors available = ', numprocs
      end if
    end subroutine init_mpi


    function getlun()
        implicit none
        integer(i4b) :: getlun
        logical(lgt) :: exists, isopen
        getlun = 9
        do
        getlun = getlun+1
        inquire(unit=getlun,exist=exists)
        if(exists) then
            inquire(unit=getlun,opened=isopen)
            if(.not. isopen) return
        end if
        end do
    end function getlun

end module init_mod
