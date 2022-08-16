module dang_util_mod
  use healpix_types
  use pix_tools
  use fitstools
  use head_fits
  use mpi
  use omp_lib
  implicit none
  
  ! Constants
  real(dp)           :: k_B     = 1.3806503d-23
  real(dp)           :: h       = 1.0545726691251021d-34 * 2.d0 * pi
  real(dp)           :: c       = 2.99792458d8
  real(dp)           :: T_CMB   = 2.7255d0
  real(dp)           :: t1, t2, t3
  real(dp)           :: nullval
  real(dp)           :: missval = -1.6375d30
  integer(i4b)       :: ierr, rank, numprocs
  integer(i4b)       :: nbands, npix, nmaps, nside, nfgs, npar
  integer(i4b)       :: npixpar, nglobalpar
  integer(i4b)       :: ncomp, ncg_groups, nsample
  integer(i4b)       :: iter, niter, ordering, nlheader
  integer(i4b)       :: proc_per_band
  integer(i4b)       :: master      = 0 
  integer(i4b)       :: from_master = 1
  integer(i4b)       :: from_worker = 2
  integer(i4b)       :: nump ! Number of unmasked pixels
  logical(lgt)       :: anynull, exist
  integer(i4b) status(mpi_status_size)
  character(len=5)                  :: iter_str
  character(len=80), dimension(180) :: header
  character(len=80), dimension(3)   :: tqu 
  character(len=128)                :: title
  character(len=10)                 :: ml_mode
  
  public    :: npix, nbands, nmaps, ordering, header, h, c, k_B, T_CMB
  public    :: ncomp, ncg_groups, nsample
  public    :: npixpar, nglobalpar, title
  public    :: iter, iter_str, exist, tqu
  public    :: ml_mode

contains 
  
  subroutine init_mpi()
    implicit none

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
    tqu(1)            = 'T'
    tqu(2)            = 'Q'
    tqu(3)            = 'U'
  end subroutine init_mpi
  
  ! Small utility for converting an integer to a string                                              
  subroutine int2string(integer, string)
    implicit none
    
    integer(i4b),     intent(in)  :: integer
    character(len=*), intent(out) :: string
    
    integer(i4b)               :: temp_int, i, k
    
    temp_int = integer
    do i = 1, len(string)
       k = temp_int / 10**(len(string)-i)
       write(string(i:i),'(I1)') k
       temp_int = temp_int - k * 10**(len(string)-i)
    end do
    
  end subroutine int2string
  
  subroutine tolower(strIn)! result(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
    ! Original author: Clive Page - switched to tolower
    
    implicit none
    
    character(len=*), intent(inout) :: strIn
    character(len=len(strIn)) :: strOut
    integer :: i,j
    
    do i = 1, len(strIn)
       j = iachar(strIn(i:i))
       if (j>= iachar("A") .and. j<=iachar("Z") ) then
          strOut(i:i) = achar(iachar(strIn(i:i))+32)
       else
          strOut(i:i) = strIn(i:i)
       end if
    end do
    
    strIn = strOut
    
  end subroutine tolower
  
  function rand_normal(mean,stdev) result(c)
    double precision :: mean,stdev,c,temp(2),theta,r
    if (stdev <= 0.0d0) then
       write(*,*) "Standard Deviation must be positive."
    else
       call RANDOM_NUMBER(temp)
       r=(-2.0d0*log(temp(1)))**0.5
       theta = 2.0d0*PI*temp(2)
       c= mean+stdev*r*sin(theta)
    end if
  end function rand_normal
  
  function eval_normal_prior(prop,mean,std) result(p)
    real(dp) :: prop, mean, std, p, var, num ,denom
    
    var   = std**2.d0
    num   = exp(-(prop-mean)**2/(2*var))
    denom = std*sqrt(2.d0*PI)
    
    p = num/denom
    
  end function eval_normal_prior
  
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
  
  subroutine write_result_map(filename, nside, ordering, header, map)
    ! Copied from map_editor - CMB/CO fortran 90 fits file tool
    implicit none
    
    character(len=*),                    intent(in)    :: filename
    integer(i4b),                        intent(in)    :: nside, ordering
    character(len=80), dimension(180),   intent(inout) :: header
    real(dp),          dimension(0:,1:), intent(in)   :: map
    
    integer(i4b) :: i, nlheader, nmaps, npix, j
    character(len=80)               :: line
    character(len=80), dimension(1) :: line2
    character(len=256)              :: outfile
    
    npix  = size(map(:,1))
    nmaps = size(map(0,:))
    
    !Modify necessary header keywords; leave all others as they were
    nlheader = size(header)
    do i = 1, nlheader
       line  = header(i)
       line2 = ''
       if (line(1:8) == 'ORDERING') then
          if (ordering == 1) then
             call add_card(line2, "ORDERING","RING", "Pixel ordering scheme, either RING or NESTED")
          else
             call add_card(line2, "ORDERING","NESTED", "Pixel ordering scheme, either RING or NESTED")
          end if
       end if
       
       if (line(1:5) == 'NSIDE') then
          call add_card(line2, "NSIDE", nside, "Resolution parameter for HEALPix")
       end if
       
       if (line(1:7) == "LASTPIX") then
          call add_card(line2, "LASTPIX", npix-1, "Last pixel # (starts at 0)")
       end if
       
       if (trim(line2(1)) /= '') header(i) = line2(1)
       
    end do
    
    outfile = trim(filename)
    
    call write_bintab(map, npix, nmaps, header, nlheader, outfile)
    
  end subroutine write_result_map
  
  function mask_avg(array,mask)
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in) :: mask
    real(dp)                           :: sum, mask_avg
    integer(i4b)                       :: mask_sum, i
    
    sum = 0.d0
    mask_sum = 0
    
    do i = 1, npix
       if (mask(i) == missval .or. mask(i) == 0.d0) then
          cycle
       else
          sum      = sum + array(i)
          mask_sum = mask_sum + 1
       end if
    end do
    
    mask_avg = sum/mask_sum
    
  end function mask_avg

  function mask_sum(array,mask)
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in) :: mask
    real(dp)                           :: sum, mask_sum
    integer(i4b)                       :: i
    
    sum = 0.d0
    
    do i = 1, npix
       if (mask(i) == missval .or. mask(i) == 0.d0) then
          cycle
       else
          sum      = sum + array(i)
       end if
    end do
    
    mask_sum = sum
    
  end function mask_sum

  function return_poltype_flag(string) result(flag)
    ! ====================================================================
    ! Here is how our bitwise poltype flagging works:
    !
    !               P  U  Q  T
    !
    ! T             0  0  0  1   = 1
    !
    !     Q + U     1  0  0  0   = 8
    !
    ! T,  Q,  U     0  1  1  1   = 7
    !
    ! T,  Q + U     1  0  0  1   = 9
    !
    ! T + Q + U     1  1  1  1   = 15
    !
    ! ====================================================================
    implicit none

    character(len=10),               intent(in) :: string
    integer(i4b)                                :: count, local_flag, nflag
    integer(i4b)                                :: i, j
    integer(i4b),     allocatable, dimension(:) :: flag
    character(len=5), allocatable, dimension(:) :: pol_list


    count = count_delimits(string,',')

    nflag = 0
    local_flag = 0
    allocate(pol_list(count+1))
    call delimit_string(string,',',pol_list)

    do i = 1, count+1
       if (pol_list(i) == 'T') then
          local_flag = local_flag + 2**0
          nflag = nflag + 1
       else if (pol_list(i) == 'Q') then
          local_flag = local_flag + 2**1
          nflag = nflag + 1
       else if (pol_list(i) == 'U') then
          local_flag = local_flag + 2**2
          nflag = nflag + 1
       else if (pol_list(i) == 'Q+U') then
          local_flag = local_flag + 2**3
          nflag = nflag + 1
       else if (pol_list(i) == 'T+Q+U') then
          local_flag = 0
          nflag = nflag + 1
       end if
    end do

    allocate(flag(nflag))
    i = 1
    do j = 0, 3
       if (iand(local_flag,2**j) .ne. 0) then
          flag(i) = 2**j
          i = i + 1
       end if
    end do
    if (iand(local_flag,0) .ne. 0) then
       flag(i) = 0
    end if

  end function return_poltype_flag

  function count_delimits(string, delimiter) result(ndel)
    implicit none
    character(len=*), intent(in) :: string, delimiter
    integer(i4b)                 :: i, j, k, ndel

    ndel = 0
    do i = 1,len_trim(string)
       if (string(i:i) == trim(delimiter)) then
          ndel = ndel + 1
       end if
    end do

  end function count_delimits

  subroutine delimit_string(string, delimiter, list)
    implicit none
    character(len=*), intent(in)                :: string, delimiter
    character(len=*), dimension(:), intent(out) :: list
    integer(i4b)                                :: i, j, k
    
    j = 1
    k = 1
    do i=1,len_trim(string)
       if (string(i:i) == trim(delimiter)) then
          list(k) = trim(string(j:i-1))
          j = i+1
          k = k + 1
       end if
    end do
    
    if (k < len(list)+1) then
       list(k) = trim(string(j:))
    end if
    
  end subroutine delimit_string

  subroutine read_map(filename,map_array)
    implicit none
    character(len=512),       intent(in)    :: filename
    real(dp), dimension(:,:), intent(inout) :: map_array

    call read_bintab(trim(filename), map_array, npix, nmaps, nullval, anynull, header=header)

  end subroutine read_map

  
end module dang_util_mod
