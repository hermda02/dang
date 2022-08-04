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
    implicit none

    type(dang_params), intent(in) :: dpar
    class(dang_cg_group), pointer :: constructor_cg
    integer(i4b),      intent(in) :: cg_group
    integer(i4b)                  :: i, count

    allocate(constructor_cg)

    constructor_cg%ncg_components = 0
    constructor_cg%cg_group       = cg_group
    constructor_cg%i_max          = dpar%cg_max_iter(cg_group)
    constructor_cg%converge       = dpar%cg_convergence(cg_group)

    do i = 1, dpar%ncomp
       if (component_list(i)%p%cg_group == cg_group) then
          constructor_cg%ncg_components = constructor_cg%ncg_components + 1
       end if
    end do

    if (constructor_cg%ncg_components == 0) then
       write(*,*) "Woah there, number of CG components = 0"
       write(*,*) "for CG group ", cg_group
       stop
    end if

    allocate(constructor_cg%cg_component(constructor_cg%ncg_components))

    count = 1
    do i = 1, dpar%ncomp
       if (component_list(i)%p%cg_group == cg_group) then
          constructor_cg%cg_component(count)%p => component_list(i)%p
          count = count + 1
       end if
    end do

  end function constructor_cg

  subroutine initialize_cg_groups(dpar)
    implicit none
    type(dang_params)             :: dpar
    integer(i4b)                  :: i, count

    allocate(cg_groups(dpar%ncggroup))

    count = 0
    do i = 1, dpar%ncggroup
       write(*,*) 'Initialize CG group ', i
       cg_groups(i)%p => dang_cg(dpar,i)
    end do
  end subroutine initialize_cg_groups

  subroutine sample_cg_groups(dpar, ddata, map_n)
    implicit none
    type(dang_data),         intent(in)    :: ddata
    type(dang_params)                      :: dpar
    integer(i4b),            intent(in)    :: map_n 

    integer(i4b)                           :: i, j
    real(dp), allocatable, dimension(:)    :: b

    do i = 1, ncg_groups
       write(*,*) "Computing a CG search of CG group ", i
       call cg_groups(i)%p%compute_rhs(ddata,b)
       call cg_groups(i)%p%cg_search(dpar,ddata,map_n,b)
       call cg_groups(i)%p%unpack_amplitudes(dpar,ddata,map_n)

       if (i == 1) then
          open(55,file='sampling_group_rhs_testing.txt')
          do j = 1, size(b)
             write(55,fmt='(2(E16.8))') b(j), cg_groups(1)%p%x(j)
          end do
          close(55)
       end if
       deallocate(b)
    end do

  end subroutine sample_cg_groups

  subroutine cg_search(self, dpar, ddata, map_n, b)

    ! Implementation of the canned algorithm (B2) outlined in Jonathan Richard Shewuck (1994)                   
    ! "An introduction to the Conjugate Gradient Method Without the Agonizing Pain"      

    implicit none

    class(dang_cg_group)                   :: self
    type(dang_data),         intent(in)    :: ddata
    type(dang_params)                      :: dpar

    integer(i4b),            intent(in)    :: map_n 
    real(dp), dimension(:),  intent(in)    :: b

    real(dp), allocatable, dimension(:)    :: r, q, d, b2, eta
    real(dp), allocatable, dimension(:)    :: x_internal

    integer(i4b)                           :: i, j, k, l, m, n
    real(dp)                               :: alpha, beta, delta_0
    real(dp)                               :: delta_old, delta_new
    real(dp)                               :: t3, t4, t5, t6

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

    ! Check mode
    if (trim(dpar%ml_mode) == 'sample') then
       ! First draw a univariate for sampling if in sampling mo
       do i = 1, n
          eta(i) = rand_normal(0.d0,1.d0)
       end do
       b2 = b + self%compute_sample_vector(ddata,eta,nbands,map_n,b)
    else if (trim(dpar%ml_mode) == 'optimize') then
       b2 = b
    end if

    ! Initial condition parameters
    !-----------------------------
    x_internal = self%x
    !-----------------------------

    r  = b2 - self%compute_Ax(ddata, x_internal, nbands, map_n)
    d  = r
    delta_new = sum(r*r)
    delta_0   = delta_new
    i = 0
    do while( (i .lt. self%i_max) .and. (delta_new .gt. self%converge))
       t3         = mpi_wtime()
       q          = self%compute_Ax(ddata, d, nbands, map_n)
       alpha      = delta_new/(sum(d*q))
       x_internal = x_internal + alpha*d

       if (mod(i,50) == 0) then
          r = b2 - self%compute_Ax(ddata, x_internal, nbands, map_n)
       else
          r = r - alpha*q
       end if

       delta_old = delta_new
       delta_new = sum(r*r)
       beta      = delta_new/delta_old
       d         = r + beta*d
       t4        = mpi_wtime()
       i         = i + 1

       ! write(*,fmt='(a,i4,a,e12.5,a,e12.5,a)') 'CG Iter: ', i, ' | delta: ', delta_new, ' | time: ', t4-t3, 's.'
       ! if (delta_new .lt. self%converge) then
       !    write(*,fmt='(a,i4,a,e12.5,a,e12.5,a)') 'Final CG Iter: ', i, ' | delta: ', delta_new, ' | time: ', t4-t3, 's.'
       ! end if

    end do

    self%x = x_internal

    deallocate(eta,b2,x_internal,r)

  end subroutine cg_search

  subroutine compute_rhs(self,ddata,b)
    implicit none

    class(dang_cg_group)                               :: self
    type(dang_data),                     intent(in)    :: ddata
    real(dp), allocatable, dimension(:), intent(inout) :: b
    integer(i4b)                                       :: component
    real(dp), allocatable, dimension(:,:,:)            :: data

    integer(i4b)                            :: i, j, k
    integer(i4b)                            :: l, m, n
    integer(i4b)                            :: offset, z

    integer(i4b)                            :: map_n

    map_n = 2 ! Dummy variable for now

    allocate(data(0:npix-1,nmaps,nbands))

    data = ddata%sig_map

    l = ddata%npix  ! the number of pixels in our sky maps 
    m = 0           ! length of the array b
    n = nbands      ! number of bands included

    ! Iterate through CG group components to construct arrays
    do i = 1, self%ncg_components
       if (self%cg_component(i)%p%type == 'template') then
          if (self%cg_component(i)%p%polfit) then
             m = m + self%cg_component(i)%p%nfit
          else
             m = m + 2*self%cg_component(i)%p%nfit
          end if
       else
          ! This needs to be adjusted to properly account for solving for
          ! different combinations of poltypes
          m = m + 2*l
       end if
    end do

    ! Allocate the RHS array, b
    allocate(b(m))
    b = 0.d0

    ! Remove other foregrounds from the data which we are fitting to
    do i = 1, ncomp
       if (component_list(i)%p%cg_group /= self%cg_group) then
          data(:,:,:) = data(:,:,:) - ddata%fg_map(:,:,1:,i)
       end if
    end do

    offset = 0
    do k = 1, self%ncg_components
       if (self%cg_component(k)%p%type /= 'template') then
          do j = 1, n
             do i = 1, l
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) then
                   b(i)   = 0.d0
                   b(l+i) = 0.d0
                   cycle
                else
                   b(i) = b(i) + (data(i-1,map_n,j)*&
                        & self%cg_component(k)%p%eval_sed(j,i-1,map_n))/&
                        & (ddata%rms_map(i-1,map_n,j)**2.d0)
                   b(l+i) = b(l+i) + (data(i-1,map_n+1,j)*&
                        & self%cg_component(k)%p%eval_sed(j,i-1,map_n+1))/&
                        & (ddata%rms_map(i-1,map_n+1,j)**2.d0)
                end if
             end do
          end do
          offset = offset + 2*l
       else
          z = 1
          do j = 1, n
             if (self%cg_component(k)%p%corr(j)) then
                do i = 1, l
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,map_n,j)**2.d0)*&
                        & data(i-1,map_n,j)*self%cg_component(k)%p%template(i-1,map_n)
                   b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,map_n+1,j)**2.d0)*&
                        & data(i-1,map_n+1,j)*self%cg_component(k)%p%template(i-1,map_n+1)
                end do
                z = z + 1
             end if
          end do
          offset = offset + self%cg_component(k)%p%nfit
       end if
    end do
  end subroutine compute_rhs

  function compute_Ax(self, ddata, x, nbands, map_n) result(res)
    implicit none
    class(dang_cg_group)                  :: self
    type(dang_data),        intent(in)    :: ddata
    real(dp), dimension(:), intent(in)    :: x
    integer(i4b),           intent(in)    :: nbands, map_n

    real(dp), allocatable, dimension(:)   :: temp1, temp2, temp3, res

    integer(i4b)                          :: i, j, k
    integer(i4b)                          :: l, m, n, offset

    ! How this works: result = sum_nu(T_nu^t N_nu^-1 T_nu)*x       
    !--------------------------------------------------------------                   
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
    !--------------------------------------------------------------                   

    n      = size(x)
    offset = 0
    
    ! do i = 1, self%ncg_components                                                     
    !    if (self%cg_component(k)%p%type /= 'template') then                           
    !       if (self%cg_component(k)%p%joint) then                                     
    !          offset = offset + 2*npix                                                 
    !       else                                                                        
    !          offset = offset + npix                                                  
    !       end if                                                                     
    !    end if                                                                        
    ! end do 

    ! THERE IS CURRENTLY NO SUPPORT FOR MULTI-TEMPLATE HANDLING

    offset = offset + 2*npix
    l = 1

    allocate(temp1(offset))
    allocate(temp2(offset))
    allocate(temp3(n))
    allocate(res(n))

    res = 0.d0

    do j = 1, nbands
       temp1 = 0.d0
       temp2 = 0.d0
       temp3 = 0.d0
       ! Solving temp1 = (T_nu)x
       do k = 1, self%ncg_components
          if (self%cg_component(k)%p%type /= 'template') then
             do i = 1, npix
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                temp1(i)      = x(i)*self%cg_component(k)%p%eval_sed(j,i-1,map_n)
                ! This is just for the joint stuff - add modularity later
                temp1(npix+i) = x(npix+i)*self%cg_component(k)%p%eval_sed(j,i-1,map_n+1)
             end do
          else
             if (self%cg_component(k)%p%corr(j)) then
                do i = 1, npix
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   temp1(i)      = temp1(i) + x(offset+l)*self%cg_component(k)%p%template(i-1,map_n)
                   temp1(npix+i) = temp1(npix+i) + x(offset+l)*&
                        & self%cg_component(k)%p%template(i-1,map_n+1)
                end do
             end if
          end if
       end do
       ! res = res + temp1
    ! end do
       ! Solving temp2 = (N^-1)temp1
       do i = 1, npix
          if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
          temp2(i)      = temp1(i)/(ddata%rms_map(i-1,map_n,j)**2.d0)
          temp2(npix+i) = temp1(npix+i)/(ddata%rms_map(i-1,map_n+1,j)**2.d0)
       end do

       ! Solving (T_nu^t)temp2
       do k = 1, self%ncg_components
          if (self%cg_component(k)%p%type /= 'template') then
             do i = 1, npix
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                temp3(i)      = temp2(i)*self%cg_component(k)%p%eval_sed(j,i-1,map_n)
                temp3(npix+i) = temp2(npix+i)*self%cg_component(k)%p%eval_sed(j,i-1,map_n)
             end do
          else
             if (self%cg_component(k)%p%corr(j)) then
                do i = 1, npix
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   temp3(offset+l) = temp3(offset+l) + temp2(i)*self%cg_component(k)%p%template(i-1,map_n)
                   temp3(offset+l) = temp3(offset+l) + temp2(npix+i)*&
                        & self%cg_component(k)%p%template(i-1,map_n+1)
                end do
                l = l + 1
             end if
          end if
       end do
       res = res + temp3
    end do
  end function compute_Ax

  function compute_sample_vector(self, ddata, eta, nbands, map_n, b) result(res)
    implicit none
    class(dang_cg_group)                  :: self
    type(dang_data),        intent(in)    :: ddata
    
    real(dp), dimension(:), intent(in)    :: eta, b
    integer(i4b),           intent(in)    :: nbands, map_n
    
    real(dp), allocatable, dimension(:)   :: temp1, temp2, res

    integer(i4b)                          :: i, j, k
    integer(i4b)                          :: offset, l, m, n

    ! How this works: result = sum_nu(T_nu^t N_nu^{-1/2})vec                                       
    !--------------------------------------------------------------                                
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
    ! res   = sum_nu(temp2)                                      |
    !--------------------------------------------------------------                                

    n      = size(eta)
    m      = size(b)
    offset = 0
    l      = 1
    
    ! do i = 1, self%ncg_components
    !    if (self%cg_component(k)%p%type /= 'template') then
    !       if (self%cg_component(k)%p%joint) then
    !          offset = offset + 2*npix
    !       else
    !          offset = offset + npix
    !       end if
    !    end if
    ! end do

    offset = offset + 2*npix

    allocate(temp1(n))
    allocate(temp2(m))
    allocate(res(m))

    res = 0.d0

    do j = 1, nbands
       temp1 = 0.d0
       temp2 = 0.d0

       do i = 1, npix
          if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
          temp1(i)      = eta(i)/(ddata%rms_map(i-1,map_n,j))
          ! This is just for the joint stuff - add modularity later
          temp1(npix+i) = eta(npix+i)/(ddata%rms_map(i-1,map_n+1,j))
       end do

       do k = 1, self%ncg_components
          if (self%cg_component(k)%p%type /= 'template') then
             do i = 1, npix
                if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                temp2(i) = temp1(i)*self%cg_component(k)%p%eval_sed(j,i-1,map_n)
                ! This is just for the joint stuff - add modularity later
                temp2(npix+i) = temp1(npix+i)*self%cg_component(k)%p%eval_sed(j,i-1,map_n+1)
             end do
          else
             if (self%cg_component(k)%p%corr(j)) then
                do i = 1, npix
                   if (ddata%masks(i-1,1) == 0.d0 .or. ddata%masks(i-1,1) == missval) cycle
                   temp2(offset+l) = temp2(offset+l) + temp1(i)*self%cg_component(k)%p%template(i-1,map_n)
                   ! This is just for the joint stuff - add modularity later
                   temp2(offset+l) = temp2(offset+l) + temp1(npix+i)*&
                        & self%cg_component(k)%p%template(i-1,map_n+1)
                end do
                l = l + 1
             end if
          end if
       end do
       res = res + temp2
    end do
 end function compute_sample_vector

  subroutine unpack_amplitudes(self, dpar, ddata, map_n)
    implicit none
    class(dang_cg_group)               :: self
    type(dang_data),     intent(in)    :: ddata
    type(dang_params)                  :: dpar
    integer(i4b),        intent(in)    :: map_n 
    
    integer(i4b)                       :: offset, i, j, k, l

    ! Let's unravel self%x such that the foreground amplitudes
    ! are properly stored
    do k = 1, self%ncg_components
       if (self%cg_component(k)%p%type /= 'template') then
          do i = 1, npix
             self%cg_component(k)%p%amplitude(i-1,map_n) = self%x(offset+i)
          end do
          offset = offset + npix
          do i = 1, npix
             self%cg_component(k)%p%amplitude(i-1,map_n+1) = self%x(offset+i)
          end do
          offset = offset + npix
       else
          l = 1
          do while (l .lt. self%cg_component(k)%p%nfit)
             do j = 1, nbands
                if (self%cg_component(k)%p%corr(j)) then
                   self%cg_component(k)%p%template_amplitudes(j,map_n)   = self%x(offset+l)
                   self%cg_component(k)%p%template_amplitudes(j,map_n+1) = self%x(offset+l)
                   l = l + 1
                end if
             end do
          end do
          offset = offset + l - 1
       end if
    end do
  end subroutine unpack_amplitudes

end module dang_cg_mod
