module dang_cg_mod
  use healpix_types
  use dang_util_mod
  use dang_param_mod
  use dang_bp_mod
  use dang_data_mod
  use dang_component_mod
  implicit none

  ! ====================================================================
  ! Now define cg_group class
  ! ====================================================================

  public cg_groups

  type :: dang_cg_group

     integer(i4b) :: ncg_components
     integer(i4b) :: cg_group
     integer(i4b) :: i_max
     real(dp)     :: converge
     integer(i4b) :: nflag
     logical(lgt) :: sample
     integer(i4b),            allocatable, dimension(:) :: pol_flag
     type(component_pointer), allocatable, dimension(:) :: cg_component

     real(dp), allocatable, dimension(:) :: x ! The amplitude vector for this CG group

   contains

     procedure :: compute_Ax
     procedure :: compute_rhs
     procedure :: compute_sample_vector
     procedure :: cg_search
     procedure :: unpack_amplitudes

  end type dang_cg_group

  interface dang_cg
     procedure constructor_cg
  end interface dang_cg

  type cg_pointer
     type(dang_cg_group), pointer :: p => null()
  end type cg_pointer

  type(cg_pointer), allocatable, dimension(:) :: cg_groups

  ! ====================================================================
contains

  function constructor_cg(dpar, cg_group)
    ! ============================================ |
    ! Construct the CG groups by pointing vital    |
    ! information (such as maximum # of iterations |
    ! and convergence criteria) to the group, and  |
    ! ensure that the CG group points to the       |
    ! correct components. See below for info on the|
    ! polarization flags.                          |
    !                                              |
    ! Inputs:                                      |
    ! dpar: class(dang_params)                     |
    !                                              |
    ! cg_group: integer                            |
    ! ============================================ |

    implicit none

    type(dang_params),  intent(in) :: dpar
    class(dang_cg_group),  pointer :: constructor_cg
    integer(i4b),       intent(in) :: cg_group
    integer(i4b)                   :: i, j, count, flag, nflag
    character(len=5), allocatable, dimension(:) :: pol_list

    allocate(constructor_cg)

    constructor_cg%ncg_components = 0
    constructor_cg%cg_group       = cg_group
    constructor_cg%i_max          = dpar%cg_max_iter(cg_group)
    constructor_cg%converge       = dpar%cg_convergence(cg_group)
    constructor_cg%sample         = dpar%cg_group_sample(cg_group)

    ! How many components in this CG group?
    do i = 1, dpar%ncomp
       if (component_list(i)%p%cg_group == cg_group) then
          constructor_cg%ncg_components = constructor_cg%ncg_components + 1
       end if
    end do

    ! Handle the error before we try to allocate
    if (constructor_cg%ncg_components == 0) then
       write(*,*) "Woah there, number of CG components = 0"
       write(*,*) "for CG group ", cg_group
       stop
    end if

    allocate(constructor_cg%cg_component(constructor_cg%ncg_components))

    ! Point to the correct component
    count = 1
    do i = 1, dpar%ncomp
       if (component_list(i)%p%cg_group == cg_group) then
          constructor_cg%cg_component(count)%p => component_list(i)%p
          count = count + 1
       end if
    end do

    constructor_cg%pol_flag = return_poltype_flag(dpar%cg_poltype(cg_group))
    constructor_cg%nflag = size(constructor_cg%pol_flag)

  end function constructor_cg

  subroutine initialize_cg_groups(dpar)
    ! ============================================ |
    ! Create an object for each of the CG groups   |
    ! and point to the correct components          |
    !                                              |
    ! dpar: class(dang_params)                     |
    ! ============================================ |
    implicit none
    type(dang_params)             :: dpar
    integer(i4b)                  :: i, count

    write(*,*) 'Initializing CG groups'

    allocate(cg_groups(dpar%ncggroup))
    count = 0
    do i = 1, dpar%ncggroup
       ! write(*,*) 'Initialize CG group ', i
       cg_groups(i)%p => dang_cg(dpar,i)
    end do
  end subroutine initialize_cg_groups

  subroutine sample_cg_groups(dpar, ddata)
    ! ===================================================== |
    ! Subroutine which samples over the CG groups           |
    ! by first computing the RHS of the matrix              |
    ! equation:                                             |
    !   sum_nu(T_nu^t N^-1 T_nu) = sum_nu(T_nu^t N^-1 d_nu) |
    ! and then compute a CG search given that group. Then   |
    ! unpack the amplitude vector in the correct way.       |
    !                                                       |
    ! Inputs:                                               |
    !   dpar: class(dang_params)                            |
    !   ddata: class(dang_data)                             |
    !                                                       |
    ! ===================================================== |
    implicit none
    type(dang_data),         intent(in)    :: ddata
    type(dang_params)                      :: dpar

    integer(i4b)                           :: i, f
    real(dp), allocatable, dimension(:)    :: b

    real(dp), allocatable, dimension(:,:)  :: x_map
    allocate(x_map(0:npix-1,1))

    do i = 1, ncg_groups
       if (cg_groups(i)%p%sample) then
          write(*,*) "Computing a CG search of CG group ", i
          do f = 1, cg_groups(i)%p%nflag
             call cg_groups(i)%p%compute_rhs(ddata,b,f)
             call cg_groups(i)%p%cg_search(dpar,ddata,b,f)
             ! x_map(0:,1) = cg_groups(i)%p%x(1:npix)
             ! call write_result_map('x_map_test1.fits', nside, ordering, header, x_map)
             ! x_map(0:,1) = cg_groups(i)%p%x(npix+1:)
             ! call write_result_map('x_map_test2.fits', nside, ordering, header, x_map)
             ! stop
             call cg_groups(i)%p%unpack_amplitudes(dpar,ddata,f)
             deallocate(b)
          end do
       end if
    end do

  end subroutine sample_cg_groups

  subroutine cg_search(self, dpar, ddata, b, flag_n)
    ! ========================================================= |
    ! Implementation of the canned algorithm (B2) outlined      |
    ! in Jonathan Richard Shewuck (1994) An Introduction to     |
    ! the Conjugate Gradient Method Without the Agonizing Pain" |    
    !                                                           |
    ! Inputs:                                                   |
    !   self: class(dang_cg_group) - the cg_group being sampled |
    !   dpar: class(dang_params)   - access to the param object |
    !   ddata: class(dang_data)    - access to the data object  |
    !   b: array(dp)               - RHS of matrix equation     |
    !   flag_n: integer            - poltype flag number        |
    !                                                           |
    ! ========================================================= |
    implicit none

    class(dang_cg_group)                   :: self
    type(dang_data),         intent(in)    :: ddata
    type(dang_params)                      :: dpar

    integer(i4b),            intent(in)    :: flag_n
    real(dp), dimension(:),  intent(in)    :: b

    real(dp), allocatable, dimension(:)    :: r, q, d, b2, eta
    real(dp), allocatable, dimension(:)    :: x_internal

    integer(i4b)                           :: i, j, k, l, m, n
    real(dp)                               :: alpha, beta, delta_0
    real(dp)                               :: delta_old, delta_new
    real(dp)                               :: t3, t4, t5, t6

    ! How big is the input RHS vector?
    m = size(b)
    
    ! Temporarily hard coded - size of N
    n = 2*npix

    ! To ensure that we only allocate the CG group amplitude vector once
    if (iter == 1) then
       allocate(self%x(m))
    end if

    allocate(eta(n))
    allocate(b2(m))
    allocate(x_internal(m))
    allocate(r(m))

    ! ! Testing block
    ! ! ========================================================= |
    ! x_internal(:) = 1.d0

    ! b2 = self%compute_Ax(ddata, x_internal, nbands, flag_n)
    
    ! write(*,*) 'First component'
    ! write(*,*) b(1), b(100), b(npix)
    ! write(*,*) b2(1), b2(100), b2(npix)

    ! write(*,*) 'Second component'
    ! write(*,*) b(npix+1), b(npix+101), b(npix+npix)
    ! write(*,*) b2(npix+1), b2(npix+101), b2(2*npix)
    
    ! x_map(0:npix-1,1) = b2(1:npix)
    ! call write_result_map('x_map_test1.fits', nside, ordering, header, x_map)
    ! x_map(npix:2*npix-1,1) = b2(npix+1:2*npix)
    ! call write_result_map('x_map_test1.fits', nside, ordering, header, x_map)
    ! stop

    ! ! ========================================================= |

    ! Check mode
    if (trim(dpar%ml_mode) == 'sample') then
       ! First draw a univariate for sampling if in sampling mode
       do i = 1, n
          eta(i) = rand_normal(0.d0,1.d0)
       end do
       b2 = b + self%compute_sample_vector(ddata,eta,nbands,b,flag_n)
    else if (trim(dpar%ml_mode) == 'optimize') then
       b2 = b
    end if

    ! Initial condition parameters
    !-----------------------------
    x_internal = self%x
    !-----------------------------

    ! Do CG sampling as described by Schewuck 1994
    r  = b2 - self%compute_Ax(ddata, x_internal, nbands, flag_n)
    d  = r
    delta_new = sum(r*r)
    delta_0   = delta_new
    i = 0

    write(*,fmt='(a,i4,a,e12.5)') 'CG Iter: ', i, ' | delta: ', delta_new
    if (delta_new .lt. self%converge) then
       write(*,fmt='(a,i4,a,e12.5)') 'Final CG Iter: ', i, ' | delta: ', delta_new
    end if
    do while( (i .lt. self%i_max) .and. (delta_new .gt. self%converge))

       t3         = mpi_wtime()
       q          = self%compute_Ax(ddata, d, nbands, flag_n)
       alpha      = delta_new/(sum(d*q))
       x_internal = x_internal + alpha*d

       if (mod(i,50) == 0) then
          r = b2 - self%compute_Ax(ddata, x_internal, nbands, flag_n)
       else
          r = r - alpha*q
       end if

       delta_old = delta_new
       delta_new = sum(r*r)
       beta      = delta_new/delta_old
       d         = r + beta*d
       t4        = mpi_wtime()
       i         = i + 1

       write(*,fmt='(a,i4,a,e12.5,a,e12.5,a)') 'CG Iter: ', i, ' | delta: ', delta_new, ' | time: ', t4-t3, 's.'
       if (delta_new .lt. self%converge) then
          write(*,fmt='(a,i4,a,e12.5,a,e12.5,a)') 'Final CG Iter: ', i, ' | delta: ', delta_new, ' | time: ', t4-t3, 's.'
       end if

    end do

    self%x = x_internal

    ! Deallocate
    deallocate(eta,b2,x_internal,r)

  end subroutine cg_search

  subroutine compute_rhs(self,ddata,b,flag_n)
    ! ===================================================== |
    ! Subroutine which compute RHS of the matrix            |
    ! equation:                                             |
    !   sum_nu(T_nu^t N^-1 T_nu) = sum_nu(T_nu^t N^-1 d_nu) |
    ! to initialize the CG group.                           |
    ! ----------------------------------------------------- |                   
    !                                                       |
    ! Inputs:                                               |
    !   self: class(dang_cg_group)                          |
    !   ddata: class(dang_data)                             |
    !   b: array(dp)           - RHS of matrix equation     |
    !   flag_n: integer        - poltype flag number        |
    !                                                       |
    ! ===================================================== |
    implicit none

    class(dang_cg_group)                               :: self
    type(dang_data),                     intent(in)    :: ddata
    real(dp), allocatable, dimension(:), intent(inout) :: b
    type(dang_comps),                    pointer       :: c
    integer(i4b)                                       :: component
    real(dp), allocatable, dimension(:,:,:)            :: data
    integer(i4b),                        intent(in)    :: flag_n

    integer(i4b)                            :: i, j, k, y
    integer(i4b)                            :: l, m, n
    integer(i4b)                            :: offset, z
    integer(i4b)                            :: map_n

    if (iand(self%pol_flag(flag_n),1) .ne. 0) then
       map_n = 1
    else if (iand(self%pol_flag(flag_n),2) .ne. 0) then
       map_n = 2
    else if (iand(self%pol_flag(flag_n),4) .ne. 0) then
       map_n = 3
    end if

    ! Here is an object we'll call data, which we will correct to be the                                                      
    ! portion of the sky signal we wish to fit to      
    allocate(data(0:npix-1,nmaps,nbands))
    data = ddata%sig_map

    l = ddata%npix  ! the number of pixels in our sky maps 
    m = 0           ! length of the array b
    n = nbands      ! number of bands included

    ! Iterate through CG group components to construct arrays
    ! Counting through the components to figure out array sizes 
    do i = 1, self%ncg_components
       if (self%cg_component(i)%p%type == 'hi_fit') then
          m = m + self%cg_component(i)%p%nfit
       else if (self%cg_component(i)%p%type == 'template') then
          if (iand(self%pol_flag(flag_n),8) .ne. 0) then
             m = m + self%cg_component(i)%p%nfit
          else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
             m = m + self%cg_component(i)%p%nfit
          end if
       else
          if (iand(self%pol_flag(flag_n),8) .ne. 0) then
             m = m + 2*l
          else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
             m = m + 3*l
          else
             m = m + l
          end if
       end if
    end do

    ! Allocate the RHS array, b
    allocate(b(m))
    b = 0.d0

    ! Remove other foregrounds from the data which we are fitting to
    do y = 1, ncomp
       c => component_list(y)%p
       if (c%cg_group /= self%cg_group) then
          if (c%type /= 'template') then
             do i = 0, npix-1
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                do k = 1, nmaps
                   do j = 1, nbands
                      data(i,k,j) = data(i,k,j) - c%amplitude(i,k)*c%eval_sed(j,i-1,k)
                   end do
                end do
             end do
          else if (c%type == 'hi_fit') then
             do i = 0, npix-1
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                do j = 1, nbands
                   data(i,1,j) = data(i,1,j) - c%eval_sed(j,i,1)*c%template_amplitudes(j,1)
                end do
             end do
          else
             do i = 0, npix-1
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                do k = 1, nmaps
                   do j = 1, nbands
                      data(i,k,j) = data(i,k,j) - c%template(i,k)*c%template_amplitudes(j,k)
                   end do
                end do
             end do
          end if
       end if
    end do

    ! Compute the RHS of the matrix equation Ax=b
    offset = 0
    do k = 1, self%ncg_components
       c => self%cg_component(k)%p
       if (c%type == 'hi_fit') then
          z = 1
          do j = 1, n
             if (c%corr(j)) then
                !!$OMP PARALLEL PRIVATE(i)
                !!$OMP DO SCHEDULE(static)
                do i = 1, l
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle

                   ! Bit flag selection for matrix building
                   if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,2,j)**2.d0)*&
                           & data(i-1,2,j)*c%eval_sed(j,i-1,2)
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,3,j)**2.d0)*&
                           & data(i-1,3,j)*c%eval_sed(j,i-1,3)
                   else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,1,j)**2.d0)*&
                           & data(i-1,1,j)*c%eval_sed(j,i-1,1)
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,2,j)**2.d0)*&
                           & data(i-1,2,j)*c%eval_sed(j,i-1,2)
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,3,j)**2.d0)*&
                           & data(i-1,3,j)*c%eval_sed(j,i-1,3)
                   else
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,map_n,j)**2.d0)*&
                           & data(i-1,map_n,j)*c%eval_sed(j,i-1,map_n)
                   end if
                end do
                !!$OMP END DO
                !!$OMP END PARALLEL
                z = z + 1
             end if
          end do
          offset = offset + c%nfit
       else if (c%type /= 'template') then
          !$OMP PARALLEL PRIVATE(i,j)
          !$OMP DO SCHEDULE(static)
          do i = 1, l
             do j = 1, n
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) then
                   if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                      b(i)     = 0.d0
                      b(l+i)   = 0.d0
                   else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                      b(i)     = 0.d0
                      b(l+i)   = 0.d0
                      b(2*l+i) = 0.d0
                   else 
                      b(i)     = 0.d0
                   end if
                   cycle
                else
                   ! Bit flag selection for matrix building
                   if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                      b(offset+i) = b(offset+i) + (data(i-1,2,j)*&
                           & c%eval_sed(j,i-1,2))/&
                           & (ddata%rms_map(i-1,2,j)**2.d0)
                      b(l+offset+i) = b(l+offset+i) + (data(i-1,3,j)*&
                           & c%eval_sed(j,i-1,3))/&
                           & (ddata%rms_map(i-1,3,j)**2.d0)
                   else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                      b(offset+i) = b(offset+i) + (data(i-1,1,j)*&
                           & c%eval_sed(j,i-1,1))/&
                           & (ddata%rms_map(i-1,1,j)**2.d0)
                      b(l+offset+i) = b(l+offset+i) + (data(i-1,2,j)*&
                           & c%eval_sed(j,i-1,2))/&
                           & (ddata%rms_map(i-1,2,j)**2.d0)
                      b(l+offset+i) = b(2*l+offset+i) + (data(i-1,3,j)*&
                           & c%eval_sed(j,i-1,3))/&
                           & (ddata%rms_map(i-1,3,j)**2.d0)
                   else
                      b(offset+i) = b(offset+i) + (data(i-1,map_n,j)*&
                           & c%eval_sed(j,i-1,map_n))/&
                           & (ddata%rms_map(i-1,map_n,j)**2.d0)
                   end if
                end if
             end do
          end do
          !$OMP END DO
          !$OMP END PARALLEL
          if (iand(self%pol_flag(flag_n),8) .ne. 0) then
             offset = offset + 2*l
          else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
             offset = offset + 3*l
          else
             offset = offset + l
          end if
       else
          z = 1
          do j = 1, n
             if (c%corr(j)) then
                !!$OMP PARALLEL PRIVATE(i)
                !!$OMP DO SCHEDULE(static)
                do i = 1, l
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   ! Bit flag selection for matrix building
                   if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,2,j)**2.d0)*&
                           & data(i-1,2,j)*c%template(i-1,2)
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,3,j)**2.d0)*&
                           & data(i-1,3,j)*c%template(i-1,3)
                   else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,1,j)**2.d0)*&
                           & data(i-1,1,j)*c%template(i-1,1)
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,2,j)**2.d0)*&
                           & data(i-1,2,j)*c%template(i-1,2)
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,3,j)**2.d0)*&
                           & data(i-1,3,j)*c%template(i-1,3)
                   else
                      b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,map_n,j)**2.d0)*&
                           & data(i-1,map_n,j)*c%template(i-1,map_n)
                   end if
                end do
                !!$OMP END DO
                !!$OMP END PARALLEL
                z = z + 1
             end if
          end do
          offset = offset + c%nfit
       end if
    end do

    deallocate(data)

  end subroutine compute_rhs

  function compute_Ax(self, ddata, x, nbands, flag_n) result(res)
    ! =========================================================== |
    ! How this works: result = sum_nu(T_nu^t N_nu^-1 T_nu)*x      |
    ! ----------------------------------------------------------- |                   
    ! Suppose T_nu is of size n, m, and input vector is of size m |                    
    ! Then N^-1 is of size n, n                                   |                     
    ! And T_nu^t is of size m, n                                  |                            
    !                                                             |                     
    ! The returned vector, like the input vector, is of size m    |                   
    !                                                             | 
    ! Let's compute this from right to left:                      |
    !                                                             |        
    ! temp1 = T_nu x                                              |
    ! temp2 = N_nu^-1 temp1                                       |
    ! temp3 = T_nu^t temp2                                        |
    ! res   = sum_nu(temp3)                                       |
    ! ----------------------------------------------------------- |                   
    !                                                             |        
    ! Inputs:                                                     |        
    !   self: class(dang_cg_group)  - the CG group to compute     | 
    !   ddata: class(dang_data)     - access to the data object   |        
    !   x: array(dp)  - vector we multiple the matrix to          |
    !   nbands: integer             - # of bands we evaluate over |
    !   flag_n: integer             - poltype flag number         |
    !                                                             |        
    ! =========================================================== |

    implicit none
    class(dang_cg_group)                  :: self
    type(dang_data),        intent(in)    :: ddata
    real(dp), dimension(:), intent(in)    :: x
    integer(i4b),           intent(in)    :: nbands, flag_n
    real(dp), allocatable,  dimension(:)  :: temp1, temp2, temp3, res
    type(dang_comps),       pointer       :: c

    integer(i4b)                          :: i, j, k, map_n
    integer(i4b)                          :: l, m, n, offset

    ! Initialize CG group stuff based off of pol_flags
    if (iand(self%pol_flag(flag_n),1) .ne. 0) then
       map_n = 1
    else if (iand(self%pol_flag(flag_n),2) .ne. 0) then
       map_n = 2
    else if (iand(self%pol_flag(flag_n),4) .ne. 0) then
       map_n = 3
    end if

    n      = size(x) ! Size of the input array
    m      = 0 ! Size of the number of pixels in the noise covariance matrix
    l      = 1 ! Used for template fit counting
    offset = 0 ! Usage below
    
    ! Count up the offset for the template handling at the end of the 
    do i = 1, self%ncg_components                                                     
       if (self%cg_component(i)%p%type == 'hi_fit') then
          offset = 0
          m      = m + npix
       else if (self%cg_component(i)%p%type /= 'template') then                           
          if (iand(self%pol_flag(flag_n),8) .ne. 0) then
             offset = offset + 2*npix
             m      = m + 2*npix
          else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
             offset = offset + 3*npix
             m      = m + 3*npix
          else
             offset = offset + npix
             m      = m + npix
          end if
       end if
    end do

    allocate(temp1(m))
    allocate(temp2(m))
    allocate(temp3(n))
    allocate(res(n))

    res = 0.d0

    do j = 1, nbands
       temp1 = 0.d0
       temp2 = 0.d0
       temp3 = 0.d0
       ! Reset the offset for each band
       offset = 0
       ! Solving temp1 = (T_nu)x
       do k = 1, self%ncg_components
          c => self%cg_component(k)%p
          if (c%type /= 'template' .and. c%type /= 'hi_fit') then
             !$OMP PARALLEL PRIVATE(i)
             !$OMP DO SCHEDULE(static)
             ! do i = 1, npix
             !    if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
             !    if (iand(self%pol_flag(flag_n),8) .ne. 0) then
             !       temp1(i)        = x(i)       *c%eval_sed(j,i-1,2)
             !       temp1(npix+i)   = x(npix+i)  *c%eval_sed(j,i-1,3)
             !    else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
             !       temp1(i)        = x(i)       *c%eval_sed(j,i-1,1)
             !       temp1(npix+i)   = x(npix+i)  *c%eval_sed(j,i-1,2)
             !       temp1(2*npix+i) = x(2*npix+i)*c%eval_sed(j,i-1,3)
             !    else
             !       temp1(i)        = x(i)       *c%eval_sed(j,i-1,map_n)
             !    end if
             ! end do
             do i = 1, npix
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                   temp1(i)        = temp1(i) + x(offset+i)       *c%eval_sed(j,i-1,2)
                   temp1(npix+i)   = temp1(i) + x(offset+npix+i)  *c%eval_sed(j,i-1,3)
                else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                   temp1(i)        = temp1(i) + x(offset+i)       *c%eval_sed(j,i-1,1)
                   temp1(npix+i)   = temp1(i) + x(offset+npix+i)  *c%eval_sed(j,i-1,2)
                   temp1(2*npix+i) = temp1(i) + x(offset+2*npix+i)*c%eval_sed(j,i-1,3)
                else
                   temp1(i)        = temp1(i) + x(offset+i)       *c%eval_sed(j,i-1,map_n)
                end if
             end do
             !$OMP END DO
             !$OMP END PARALLEL
             if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                offset = offset + 2*npix
             else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                offset = offset + 3*npix
             else
                offset = offset + npix
             end if
          else if (c%type == 'hi_fit') then
             if (c%corr(j)) then
                do i = 1, npix
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   temp1(i) = temp1(i) + x(j)*c%eval_sed(j,i-1,map_n)
                end do
             end if
          else if (c%type == 'template') then
             if (c%corr(j)) then
                !$OMP PARALLEL PRIVATE(i)
                !$OMP DO SCHEDULE(static)
                do i = 1, npix
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                      temp1(i)      = temp1(i)      + x(offset+l)*c%template(i-1,2)
                      temp1(npix+i) = temp1(npix+i) + x(offset+l)*c%template(i-1,3)
                   else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                      temp1(i)      = temp1(i)      + x(offset+l)*c%template(i-1,1)
                      temp1(npix+i) = temp1(npix+i) + x(offset+l)*c%template(i-1,2)
                      temp1(npix+i) = temp1(npix+i) + x(offset+l)*c%template(i-1,3)
                   else
                      temp1(i)      = temp1(i) + x(offset+l)*c%template(i-1,map_n)
                   end if
                end do
                !$OMP END DO
                !$OMP END PARALLEL
             end if
          end if
       end do

       ! Solving temp2 = (N^-1)temp1
       !$OMP PARALLEL PRIVATE(i)
       !$OMP DO SCHEDULE(static)
       do i = 1, npix
          if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
          if (iand(self%pol_flag(flag_n),8) .ne. 0) then
             temp2(i)        = temp1(i)/(ddata%rms_map(i-1,2,j)**2.d0)
             temp2(npix+i)   = temp1(npix+i)/(ddata%rms_map(i-1,3,j)**2.d0)
          else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
             temp2(i)        = temp1(i)/(ddata%rms_map(i-1,1,j)**2.d0)
             temp2(npix+i)   = temp1(npix+i)/(ddata%rms_map(i-1,2,j)**2.d0)
             temp2(2*npix+i) = temp1(npix+i)/(ddata%rms_map(i-1,3,j)**2.d0)
          else
             temp2(i)        = temp1(i)/(ddata%rms_map(i-1,map_n,j)**2.d0)
          end if
       end do
       !$OMP END DO
       !$OMP END PARALLEL

       ! Reset the offset for each multiplication
       offset = 0

       ! Solving (T_nu^t)temp2
       do k = 1, self%ncg_components
          c => self%cg_component(k)%p
          if (c%type /= 'template' .and. c%type /= 'hi_fit') then
             !$OMP PARALLEL PRIVATE(i)
             !$OMP DO SCHEDULE(static)
             ! do i = 1, npix
             !    if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
             !    if (iand(self%pol_flag(flag_n),8) .ne. 0) then
             !       temp3(i)        = temp2(i)*c%eval_sed(j,i-1,2)
             !       temp3(npix+i)   = temp2(npix+i)*c%eval_sed(j,i-1,3)
             !    else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
             !       temp3(i)        = temp2(i)*c%eval_sed(j,i-1,1)
             !       temp3(npix+i)   = temp2(npix+i)*c%eval_sed(j,i-1,2)
             !       temp3(2*npix+i) = temp2(2*npix+i)*c%eval_sed(j,i-1,3)
             !    else
             !       temp3(i)        = temp2(i)*c%eval_sed(j,i-1,map_n)
             !    end if
             ! end do
             do i = 1, npix
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                   temp3(offset+i)        = temp2(i)*c%eval_sed(j,i-1,2)
                   temp3(offset+npix+i)   = temp2(npix+i)*c%eval_sed(j,i-1,3)
                else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                   temp3(offset+i)        = temp2(i)*c%eval_sed(j,i-1,1)
                   temp3(offset+npix+i)   = temp2(npix+i)*c%eval_sed(j,i-1,2)
                   temp3(offset+2*npix+i) = temp2(2*npix+i)*c%eval_sed(j,i-1,3)
                else
                   temp3(offset+i)        = temp2(i)*c%eval_sed(j,i-1,map_n)
                end if
             end do
             !$OMP END DO
             !$OMP END PARALLEL
             if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                offset = offset + 2*npix
             else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                offset = offset + 3*npix
             else
                offset = offset + npix
             end if
          else if (c%type == 'hi_fit') then
             if (c%corr(j)) then
                do i = 1, npix
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   temp3(j) = temp3(j) + temp2(i)*c%eval_sed(j,i-1,1)
                end do
             end if
          else if (c%type == 'template') then
             if (c%corr(j)) then
                !!$OMP PARALLEL PRIVATE(i)
                !!$OMP DO SCHEDULE(static)
                do i = 1, npix
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                      temp3(offset+l) = temp3(offset+l) + temp2(i)*c%template(i-1,2)
                      temp3(offset+l) = temp3(offset+l) + temp2(npix+i)*&
                           & c%template(i-1,3)
                   else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                      temp3(offset+l) = temp3(offset+l) + temp2(i)*c%template(i-1,1)
                      temp3(offset+l) = temp3(offset+l) + temp2(npix+i)*&
                           & c%template(i-1,2)
                      temp3(offset+l) = temp3(offset+l) + temp2(npix+i)*&
                           & c%template(i-1,3)
                   else
                      temp3(offset+l) = temp3(offset+l) + temp2(i)*c%template(i-1,map_n)
                   end if
                end do
                !!$OMP END DO
                !!$OMP END PARALLEL
                l = l + 1
             end if
          end if
       end do
       res = res + temp3
    end do
  end function compute_Ax

  function compute_sample_vector(self, ddata, eta, nbands, b, flag_n) result(res)
    ! =========================================================== |
    ! How this works: result = sum_nu(T_nu^t N_nu^{-1/2})eta      | 
    ! ----------------------------------------------------------- |        
    ! Suppose T_nu is of size n, m, and input vector is of size m |                                
    ! Then N^{-1/2} is of size n, n                               |                                
    ! And T_nu^t is of size m, n                                  |                                
    !                                                             |                                
    ! The returned vector, like the input vector, is of size m    |   
    !                                                             |                                
    ! Let's compute this from right to left:                      |
    !                                                             |                                
    ! temp1 = N_nu^{-1/2}vec                                      |
    ! temp2 = T_nu^t temp1                                        |
    ! res   = sum_nu(temp2)                                       |
    ! ----------------------------------------------------------- |                                
    !                                                             |                                
    ! Inputs:                                                     |        
    !   self: class(dang_cg_group)  - the CG group to compute     | 
    !   ddata: class(dang_data)     - access to the data object   |        
    !   eta: array(dp)  - vector of univariates for sampling term |
    !   nbands: integer             - # of bands we evaluate over |
    !   b: array(dp)                - RHS of matrix equation      |
    !   flag_n: integer             - poltype flag number         |
    !                                                             |                                
    ! Outputs:                                                    |                                
    !                                                             |                                
    !   res: array(dp)              - resulting vector            |                                
    !                                                             |                                
    ! =========================================================== |
    implicit none
    class(dang_cg_group)                  :: self
    type(dang_data),        intent(in)    :: ddata
    
    real(dp), dimension(:), intent(in)    :: eta, b
    integer(i4b),           intent(in)    :: nbands, flag_n
    
    real(dp), allocatable, dimension(:)   :: temp1, temp2, res
    type(dang_comps),       pointer       :: c

    integer(i4b)                          :: i, j, k, map_n
    integer(i4b)                          :: offset, l, m, n


    ! Initialize CG group stuff based off of pol_flags
    if (iand(self%pol_flag(flag_n),1) .ne. 0) then
       map_n = 1
    else if (iand(self%pol_flag(flag_n),2) .ne. 0) then
       map_n = 2
    else if (iand(self%pol_flag(flag_n),4) .ne. 0) then
       map_n = 3
    end if

    n      = size(eta)
    m      = size(b)
    offset = 0
    l      = 1
    
    ! Count up the offset for the template handling at the end of the 
    do i = 1, self%ncg_components                                                     
       if (self%cg_component(i)%p%type /= 'template') then                           
          if (iand(self%pol_flag(flag_n),8) .ne. 0) then
             offset = offset + 2*npix
          else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
             offset = offset + 3*npix
          else
             offset = offset + npix
          end if
       end if
    end do

    ! offset = offset + 2*npix

    allocate(temp1(n))
    allocate(temp2(m))
    allocate(res(m))

    res = 0.d0

    do j = 1, nbands
       temp1 = 0.d0
       temp2 = 0.d0

       !!$OMP PARALLEL PRIVATE(i)
       !!$OMP DO SCHEDULE(static)
       ! Solving temp1 = N^{-1/2}eta
       do i = 1, npix
          if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
          if (iand(self%pol_flag(flag_n),8) .ne. 0) then
             temp1(i)        = eta(i)/(ddata%rms_map(i-1,2,j))
             temp1(npix+i)   = eta(npix+i)/(ddata%rms_map(i-1,3,j))
          else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
             temp1(i)        = eta(i)/(ddata%rms_map(i-1,1,j))
             temp1(npix+i)   = eta(npix+i)/(ddata%rms_map(i-1,2,j))
             temp1(2*npix+i) = eta(2*npix+i)/(ddata%rms_map(i-1,3,j))
          else
             temp1(i)        = eta(i)/(ddata%rms_map(i-1,map_n,j))
          end if
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL

       ! Solving temp2 = T^t temp1
       do k = 1, self%ncg_components
          c => self%cg_component(k)%p
          if (c%type /= 'template' .and. c%type /= 'hi_fit') then
             !!$OMP PARALLEL PRIVATE(i)
             !!$OMP DO SCHEDULE(static)
             do i = 1, npix
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                   temp2(i)        = temp1(i)*c%eval_sed(j,i-1,2)
                   temp2(npix+i)   = temp1(npix+i)*c%eval_sed(j,i-1,3)
                else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                   temp2(i)        = temp1(i)*c%eval_sed(j,i-1,1)
                   temp2(npix+i)   = temp1(npix+i)*c%eval_sed(j,i-1,2)
                   temp2(2*npix+i) = temp1(2*npix+i)*c%eval_sed(j,i-1,3)
                else
                   temp2(i)        = temp1(i)*c%eval_sed(j,i-1,map_n)
                end if
             end do
             !!$OMP END DO
             !!$OMP END PARALLEL
          else if (c%type == 'hi_fit') then
             if (c%corr(j)) then
                do i = 1, npix
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   temp2(j) = temp2(j) + temp1(i)*c%eval_sed(j,i-1,map_n)
                end do
             end if
          else if (c%type == 'template') then
             if (c%corr(j)) then
                !!$OMP PARALLEL PRIVATE(i)
                !!$OMP DO SCHEDULE(static)
                do i = 1, npix
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                      temp2(offset+l) = temp2(offset+l) + temp1(i)*c%template(i-1,2)
                      temp2(offset+l) = temp2(offset+l) + temp1(npix+i)*c%template(i-1,3)
                   else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                      temp2(offset+l) = temp2(offset+l) + temp1(i)*c%template(i-1,1)
                      temp2(offset+l) = temp2(offset+l) + temp1(npix+i)*c%template(i-1,2)
                      temp2(offset+l) = temp2(offset+l) + temp1(2*npix+i)*c%template(i-1,3)
                   else
                      temp2(offset+l) = temp2(offset+l) + temp1(i)*c%template(i-1,map_n)
                   end if
                end do
                !!$OMP END DO
                !!$OMP END PARALLEL
                l = l + 1
             end if
          end if
       end do
       res = res + temp2
    end do
 end function compute_sample_vector

 subroutine unpack_amplitudes(self, dpar, ddata, flag_n)
    ! =========================================================== |
    !  This routine walks through the resulting amplitude vector, |                                
    !  storing the amplitude results in the correct foregrounds.  | 
    !                                                             |                                
    ! ----------------------------------------------------------- |                                
    !                                                             |                                
    ! Inputs:                                                     |                                
    !   self: class(dang_cg_group) - the cg_group being sampled   |
    !   dpar: class(dang_params)   - access to the param object   |
    !   ddata: class(dang_data)    - access to the data object    |
    !   flag_n: integer            - poltype flag number          |
    !                                                             |                                
    ! =========================================================== |
    implicit none
    class(dang_cg_group)               :: self
    type(dang_data),     intent(in)    :: ddata
    integer(i4b),        intent(in)    :: flag_n
    type(dang_params)                  :: dpar
    type(dang_comps),    pointer       :: c
    
    integer(i4b)                       :: offset, i, j, k, l
    integer(i4b)                       :: map_n 


    ! Initialize CG group stuff based off of pol_flags
    if (iand(self%pol_flag(flag_n),1) .ne. 0) then
       map_n = 1
    else if (iand(self%pol_flag(flag_n),2) .ne. 0) then
       map_n = 2
    else if (iand(self%pol_flag(flag_n),4) .ne. 0) then
       map_n = 3
    end if

    offset = 0

    ! Let's unravel self%x such that the foreground amplitudes
    ! are properly stored
    do k = 1, self%ncg_components
       c => self%cg_component(k)%p
       if (c%type /= 'template' .and. c%type /= 'hi_fit') then
          ! Unravel according to polarization flags
          if (iand(self%pol_flag(flag_n),8) .ne. 0) then
             do i = 1, npix
                c%amplitude(i-1,2) = self%x(offset+i)
             end do
             offset = offset + npix
             do i = 1, npix
                c%amplitude(i-1,3) = self%x(offset+i)
             end do
             offset = offset + npix
          else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
             do i = 1, npix
                c%amplitude(i-1,1) = self%x(offset+i)
             end do
             offset = offset + npix
             do i = 1, npix
                c%amplitude(i-1,2) = self%x(offset+i)
             end do
             offset = offset + npix
             do i = 1, npix
                c%amplitude(i-1,3) = self%x(offset+i)
             end do
             offset = offset + npix
          else
             do i = 1, npix
                c%amplitude(i-1,map_n) = self%x(offset+i)
             end do
             offset = offset + npix
          end if
       else if (c%type == 'hi_fit') then
          l = 1
          do while (l .lt. c%nfit)
             do j = 1, nbands
                if (c%corr(j)) then
                   c%template_amplitudes(j,map_n) = self%x(offset+l)
                   l = l + 1
                end if
             end do
          end do
       else if (c%type == 'template') then
          l = 1
          do while (l .lt. c%nfit)
             do j = 1, nbands
                if (c%corr(j)) then
                   if (iand(self%pol_flag(flag_n),8) .ne. 0) then
                      c%template_amplitudes(j,2)     = self%x(offset+l)
                      c%template_amplitudes(j,3)     = self%x(offset+l)
                   else if (iand(self%pol_flag(flag_n),0) .ne. 0) then
                      c%template_amplitudes(j,1)     = self%x(offset+l)
                      c%template_amplitudes(j,2)     = self%x(offset+l)
                      c%template_amplitudes(j,3)     = self%x(offset+l)
                   else
                      c%template_amplitudes(j,map_n) = self%x(offset+l)
                   end if
                   l = l + 1
                end if
             end do
          end do
          offset = offset + l - 1
       end if
    end do
  end subroutine unpack_amplitudes

end module dang_cg_mod
