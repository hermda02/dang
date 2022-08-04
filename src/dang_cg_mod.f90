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
     type(component_pointer), allocatable, dimension(:) :: cg_component

   contains

     procedure :: compute_rhs
     procedure :: compute_sample_vector

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

  subroutine compute_rhs(self,ddata,b)
    implicit none

    class(dang_cg_group)                               :: self
    type(dang_data),                     intent(in)    :: ddata
    type(dang_params)                                  :: dpar
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
       write(*,*) self%cg_component(i)%p%type
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
    do i = 1, dpar%ncomp
       if (component_list(i)%p%cg_group /= self%cg_group) then
          data(:,:,:) = data(:,:,:) - ddata%fg_map(:,:,1:,i)
       end if
    end do

    write(*,*) data(0,2,1)
    write(*,*) self%cg_component(1)%p%eval_sed(1,0,2)
    write(*,*) ddata%rms_map(0,2,1)
    write(*,*) ddata%temps(0,2,1)

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
                        & data(i-1,map_n,j)*ddata%temps(i-1,map_n,1)
                   b(offset+z) = b(offset+z) + 1.d0/(ddata%rms_map(i-1,map_n+1,j)**2.d0)*&
                        & data(i-1,map_n+1,j)*ddata%temps(i-1,map_n+1,1)
                end do
                z = z + 1
             end if
          end do
          offset = offset + self%cg_component(k)%p%nfit
       end if
    end do
  end subroutine compute_rhs

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
                   temp2(offset+l) = temp2(offset+l) + temp1(i)*ddata%temps(i-1,map_n,1)
                   ! This is just for the joint stuff - add modularity later
                   temp2(offset+l) = temp2(offset+l) + temp1(npix+i)*ddata%temps(i-1,map_n+1,1)
                end do
                l = l + 1
             end if
          end if
       end do
       res = res + temp2
    end do


  end function compute_sample_vector

  function constructor_cg(dpar, cg_group)
    implicit none

    type(dang_params), intent(in) :: dpar
    class(dang_cg_group), pointer :: constructor_cg
    integer(i4b),      intent(in) :: cg_group

    integer(i4b)                  :: i, count

    allocate(constructor_cg)

    constructor_cg%ncg_components = 0
    constructor_cg%cg_group       = cg_group

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


end module dang_cg_mod
