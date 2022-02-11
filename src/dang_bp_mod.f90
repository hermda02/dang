module dang_bp_mod
  use healpix_types
  use dang_util_mod
  use dang_param_mod
 implicit none

  type bandinfo
     character(len=20) :: label, id, unit
     integer(i4b)      :: n
     real(dp)          :: nu_c
     real(dp), allocatable, dimension(:) :: nu0, nu, tau0, tau
  end type bandinfo

  integer(i4b)                              :: numband
  type(bandinfo), allocatable, dimension(:) :: bp
  

contains

  subroutine init_bp_mod(dpar)
    implicit none
    class(dang_params), intent(inout) :: dpar
    character(len=512)                :: paramfile
    character(len=3)                  :: itext
    real(dp)                          :: threshold

    integer(i4b) :: i, j, unit
  
    call getarg(1,paramfile)
    unit = getlun()

    allocate(bp(dpar%numinc))

    do i = 1, dpar%numinc
       bp(i)%nu_c = dpar%band_nu(i)
       if (bp(i)%nu_c < 1d9) then
          bp(i)%nu_c = bp(i)%nu_c*1d9
       end if
       bp(i)%id   = dpar%bp_id(i)
       if (trim(dpar%bp_id(i)) == 'delta') then
          dpar%bp_file(j) = ''
       else if (trim(dpar%bp_id(i)) == 'LFI') then
          threshold = 0.d0
       else if (trim(dpar%bp_id(i)) == 'WMAP') then
         threshold = 0.d0
       else if (trim(dpar%bp_id(i)) == 'HFI_cmb') then
         threshold = 1.d-7
       else if (trim(dpar%bp_id(i)) == 'DIRBE') then
         threshold = 0.d0
       end if
       if (trim(dpar%bp_file(i)) /= '') then
          call read_bandpass(trim(dpar%bp_file(i)), threshold, bp(i)%n, bp(i)%nu0, bp(i)%tau0)
          bp(i)%tau0 = normalize_bandpass(bp(i)%tau0)
          allocate(bp(i)%nu(bp(i)%n), bp(i)%tau(bp(i)%n))
       end if
    end do

  end subroutine init_bp_mod

  function normalize_bandpass(tau_in) result(tau_out)
    implicit none
    
    real(dp), dimension(:), intent(in)   :: tau_in
    real(dp), dimension(:), allocatable  :: tau_out
    real(dp)     :: total
    integer(i4b) :: i, leng

    leng = size(tau_in)

    allocate(tau_out(leng))

    tau_out = tau_in

    total = sum(tau_in)
    do i = 1, leng
       tau_out(i) = tau_in(i)/total
    end do

  end function normalize_bandpass
  
  subroutine read_bandpass(filename, threshold, n, nu, tau)
    ! Copied from Commander1:
    !
    ! https://github.com/Cosmoglobe/Commander1/tree/main/src/commander/comm_bp_mod.f90
    implicit none
    
    character(len=*),                    intent(in)  :: filename
    real(dp),                            intent(in)  :: threshold
    integer(i4b),                        intent(out) :: n
    real(dp), allocatable, dimension(:), intent(out) :: nu, tau
    
    integer(i4b)        :: i, j, unit, first, last, m
    logical(lgt)        :: exist
    character(len=128)  :: string
    real(dp), allocatable, dimension(:) :: x, y
    
    unit = getlun()
    
    inquire(file=trim(filename), exist=exist)
    if (.not. exist) then
       write(*,*) 'Bandpass file does not exist = ', trim(filename)
       stop
    end if
    
    ! Find the number of entries                                                                                                                                                                               
    m = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,*,end=1) string
       if (string(1:1)=='#') cycle
       m = m + 1
    end do
1   close(unit)
    
    if (m == 0) then
       write(*,*) 'No valid data entries in spectrum file ', trim(filename)
       stop
    end if
    
    allocate(x(m), y(m))
    m = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,fmt='(a)',end=2) string
       if (string(1:1)=='#') cycle
       m = m+1
       read(string,*) x(m), y(m)
       
       ! Drop double entries                                                                                                                                                                                   
       if (m > 1) then
          if (x(m) == x(m-1)) m = m-1
       end if
    end do
2   close(unit)
    
    x = x * 1.d9 ! Convert from GHz to Hz                                                                                                                                                                      
    
    first = 1
    last  = m
    if (threshold > 0.d0) then
       do while (y(first) < threshold*maxval(y(1:m)))
          first = first+1
       end do
       do while (y(last) < threshold*maxval(y(1:m)))
          last = last-1
       end do
    end if
    
    n = last-first+1
    allocate(nu(n), tau(n))
    nu  = x(first:last)
    tau = y(first:last)
    
    deallocate(x, y)
    
    
  end subroutine read_bandpass
  
end module dang_bp_mod
