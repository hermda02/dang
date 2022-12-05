module dang_data_mod
  use healpix_types
  use dang_util_mod
  use dang_param_mod
  use dang_component_mod
  use dang_bp_mod
  implicit none

  type, public                                :: dang_data

     ! Storage
    character(len=255)                        :: label
    
    integer(i4b)                              :: npix
    integer(i4b)                              :: nfgs
    integer(i4b)                              :: nmaps
    real(dp)                                  :: chisq

    real(dp), allocatable, dimension(:,:,:)   :: sig_map    ! data maps 
    real(dp), allocatable, dimension(:,:,:)   :: rms_map    ! noise maps
    real(dp), allocatable, dimension(:,:,:)   :: res_map    ! Residual maps
    real(dp), allocatable, dimension(:,:)     :: chi_map    ! Chisq map (for outputs)
    real(dp), allocatable, dimension(:,:,:)   :: sky_model  ! Total sky model

    real(dp), allocatable, dimension(:)       :: conversion ! Store the conversion factor from native to uK_RJ

    real(dp), allocatable, dimension(:)       :: gain       ! Where band gains are stored
    real(dp), allocatable, dimension(:)       :: offset     ! Where band offsets are stored
    
    real(dp), allocatable, dimension(:,:)     :: masks      ! Where masks are stored

    real(dp), allocatable, dimension(:,:,:)   :: temps      ! Where template maps are stored   
    real(dp), allocatable, dimension(:,:,:)   :: temp_amps  ! Where template amplitudes are stored
    real(dp), allocatable, dimension(:,:)     :: temp_norm  ! Where we store template normalizations

    real(dp), allocatable, dimension(:)       :: band_chisq ! A variable to store the chisq of each band
    integer(i4b), allocatable, dimension(:)   :: band_calib

    logical(lgt), allocatable, dimension(:)   :: fit_gain
    logical(lgt), allocatable, dimension(:)   :: fit_offset

  contains
    procedure :: init_data_maps
    procedure :: read_data_maps
    procedure :: update_sky_model
    procedure :: convert_maps
    procedure :: mask_hi_threshold
    procedure :: read_band_offsets
  end type dang_data

  private :: i, j, k, l
  integer(i4b) :: i, j, k, l

contains

  subroutine init_data_maps(self,dpar)
    implicit none
    class(dang_data),       intent(inout) :: self
    type(dang_params),      intent(in)    :: dpar

    integer(i4b)                          :: i

    write(*,*) "Initializing data maps"
    
    allocate(self%masks(0:npix-1,nmaps))

    allocate(self%sig_map(0:npix-1,nmaps,nbands))
    allocate(self%rms_map(0:npix-1,nmaps,nbands))
    allocate(self%res_map(0:npix-1,nmaps,nbands))
    allocate(self%chi_map(0:npix-1,nmaps))

    allocate(self%sky_model(0:npix-1,nmaps,nbands))
    self%sky_model(:,:,:) = 0.d0

    allocate(self%temps(0:npix-1,nmaps,dpar%ntemp))
    allocate(self%temp_amps(nbands,nmaps,dpar%ntemp))
    allocate(self%temp_norm(nmaps,dpar%ntemp))
    self%temp_amps(:,:,:) = 0.d0

    allocate(self%conversion(nbands))
   
    allocate(self%gain(nbands))
    allocate(self%offset(nbands))
    allocate(self%fit_gain(nbands))
    allocate(self%fit_offset(nbands))
    allocate(self%band_chisq(nbands))
    
    self%gain   = 1.d0
    self%offset = 0.d0

    self%fit_gain = dpar%fit_gain
    self%fit_offset = dpar%fit_offs

  end subroutine init_data_maps

  subroutine read_data_maps(self, dpar)
    implicit none
    class(dang_data),       intent(inout) :: self
    type(dang_params),      intent(in)    :: dpar
    integer(i4b)                          :: i, j
    real(dp), allocatable, dimension(:,:) :: map, rms
    
    allocate(map(0:npix-1,nmaps))
    allocate(rms(0:npix-1,nmaps))
    
    ! Read maskfile
    call read_bintab(dpar%mask_file,self%masks,npix,1,nullval,anynull,header=header)
    do i = 0, npix-1
       do j = 1, nmaps
          if (self%masks(i,j) == 0.d0 .or. self%masks(i,j) == missval) then
             self%masks(i,j) = missval
          else 
             nump = nump + 1
          end if
       end do
    end do
    
    ! Read maps in
    do j = 1, nbands
       ! If not supposed to be swapped BP maps, load that map
       if (.not. dpar%bp_map(j)) then
          call read_bintab(trim(dpar%datadir) // trim(dpar%band_noisefile(j)), &
               rms,npix,nmaps,nullval,anynull,header=header)
          self%rms_map(:,:,j) = rms
          call read_bintab(trim(dpar%datadir) // trim(dpar%band_mapfile(j)), &
               map,npix,nmaps,nullval,anynull,header=header)
          self%sig_map(:,:,j) = map
       end if
    end do
    
    deallocate(map,rms)
  
    ! Read templates
    do i = 1, dpar%ntemp
       call read_bintab(dpar%temp_file(i), self%temps(:,:,i), npix,nmaps,nullval,anynull,header=header)
       ! Normalize templates
       do k = 1, nmaps
          self%temp_norm(k,i) = 1.d0!maxval(self%temps(:,k,i))
          self%temps(:,k,i)   = self%temps(:,k,i)/self%temp_norm(k,i)
       end do
    end do
  end subroutine read_data_maps

  subroutine read_band_offsets(self,dpar)
    implicit none
    class(dang_data),       intent(inout) :: self
    type(dang_params),      intent(in)    :: dpar
    character(len=512)                    :: file
    character(len=128)                    :: fmt, nband_str, band

    real(dp)                              :: offset

    logical(lgt)                          :: exist
    logical(lgt), allocatable, dimension(:) :: loaded

    integer(i4b)                          :: i, j, k
    integer(i4b)                          :: unit, ios, ierror

    allocate(loaded(nbands))

    write(nband_str, '(i4)') nbands

    file = trim(dpar%offset_file)

    unit = getlun()
    ierror  = 0
    fmt  = '('//trim(nband_str)//'(E17.8))'

    if (trim(file) == '') then
       write(*,*) 'No BAND_OFFSET_FILE -- setting all offsets to 0'
       self%offset(:) = 0.d0
    else
       open(unit,file=trim(dpar%datadir)//file)
       do while (ierror .eq. 0) 
          read(unit=unit,fmt=*,iostat=ierror) band, offset
          do j = 1, nbands
             if (trim(band) == trim(dpar%band_label(j))) then
                self%offset(j) = offset
                loaded(j) = .true.
             end if
          end do
       end do
    end if
    do j = 1, nbands
       if (.not. loaded(j)) then
          write(*,*) trim(dpar%band_label(j))//' offset not loaded -- set to 0'
       end if
    end do
    
  end subroutine read_band_offsets

  subroutine update_sky_model(self)
    implicit none
    class(dang_data),  intent(inout) :: self
    type(dang_comps),  pointer       :: c
    integer(i4b)                     :: i, j, k, l

    self%sky_model(:,:,:) = 0.d0
    do l = 1, ncomp
       c => component_list(l)%p
       do i = 0, npix-1
          do k = 1, nmaps
             do j = 1, nbands
                self%sky_model(i,k,j) = self%sky_model(i,k,j) + c%eval_signal(j,i,k)
             end do
          end do
       end do
    end do

    ! Calculate the residual
    do i = 0, npix-1
       do k = 1, nmaps
          do j = 1, nbands
             if (k == 1) then
                ! For intensity
                self%res_map(i,1,j) = (self%sig_map(i,1,j)-self%offset(j))/self%gain(j)-self%sky_model(i,1,j)
             else
                ! and polarization
                self%res_map(i,2:3,j) = self%sig_map(i,2:3,j)-self%sky_model(i,2:3,j)
             end if
          end do
       end do
    end do

  end subroutine update_sky_model

  subroutine mask_hi_threshold(self, dpar)
    implicit none
    class(dang_data),            intent(inout) :: self
    type(dang_params)                          :: dpar
    type(dang_comps),   pointer                :: c
    integer(i4b)                               :: i, j

    do i = 1, ncomp
       if (component_list(i)%p%type == 'hi_fit') then
          c => component_list(i)%p
          ! call ddata%mask_hi_threshold(dpar,component_list(i)%p)
       end if
    end do
    do i = 0, npix-1
       if (c%template(i,1) > dpar%thresh) then
          self%masks(i,1) = 0.d0
       else if (self%masks(i,1) == missval) then
          self%masks(i,1) = 0.d0
       else if (self%rms_map(i,1,1) == 0.d0) then
          self%masks(i,1) = 0.d0
       else
          self%masks(i,1) = 1.d0
       end if
    end do
    ! nump = 0
    ! do i = 0, npix-1
    !    do j = 1, nmaps
    !       if (self%masks(i,j) == 0.d0 .or. self%masks(i,j) == missval) then
    !          self%masks(i,j) = 0.d0
    !       else 
    !          nump = nump + 1
    !       end if
    !    end do
    ! end do
  end subroutine mask_hi_threshold
  
  ! subroutine dust_correct_band(self,dpar,comp,band,iter)
  !   implicit none
  !   type(dang_data),             intent(inout) :: self
  !   type(dang_params)                          :: dpar
  !   type(dang_comps)                           :: comp
  !   integer(i4b),                intent(in)    :: band
  !   integer(i4b), optional,      intent(in)    :: iter
  !   real(dp), allocatable, dimension(:,:,:)    :: thermal_map
  !   integer(i4b)                               :: i, j, k
  !   character(len=256)                         :: title

  !   real(dp)                                   :: T_d, beta

  !   allocate(thermal_map(0:npix-1,nmaps,nbands))

  !   write(*,'(a,a)') 'Dust correcting band ', trim(dpar%band_label(band))
  !   if (trim(dpar%dust_corr_type) == 'uniform') then
  !      T_d    = dpar%mbb_gauss(1,1)
  !      beta   = dpar%mbb_gauss(2,1)
  !      do i = 0, npix-1
  !         do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
  !            thermal_map(i,k,band) = self%temps(i,k,1)*evaluate_mbb(bp(band),353.d9,T_d,beta)
  !            self%sig_map(i,k,band) = self%sig_map(i,k,band) - thermal_map(i,k,band)
  !         end do
  !      end do
  !   else if (trim(dpar%dust_corr_type) == 'sample') then
  !      do i = 0, npix-1
  !         do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
  !            if (dpar%mbb_gauss(1,2) .gt. 0.d0) then
  !               T_d    = rand_normal(dpar%mbb_gauss(1,1),dpar%mbb_gauss(1,2))
  !            else 
  !               T_d    = dpar%mbb_gauss(1,1)
  !            end if
  !            if (dpar%mbb_gauss(2,2) .gt. 0.d0) then
  !               beta = rand_normal(dpar%mbb_gauss(2,1),dpar%mbb_gauss(2,2))
  !            else
  !               beta = dpar%mbb_gauss(2,1)
  !            end if
  !            thermal_map(i,k,band) = self%temps(i,k,1)*evaluate_mbb(bp(band),353.d9,T_d,beta)
  !            self%sig_map(i,k,band) = self%sig_map(i,k,band) - thermal_map(i,k,band)
  !         end do
  !      end do
  !   else if (trim(dpar%dust_corr_type) == 'planck') then
  !      stop
  !   end if


  !   if (present(iter)) then
  !      write(iter_str, '(i0.5)') iter
  !      title = trim(dpar%outdir)//trim(dpar%band_label(band))//'_thermal_map_k' // trim(iter_str) // '.fits'
  !   else
  !      title = trim(dpar%outdir)//trim(dpar%band_label(band))//'_thermal_map.fits'
  !   end if
  !   call write_result_map(trim(title), nside, ordering, header, thermal_map(:,:,band))
  ! end subroutine dust_correct_band

  subroutine convert_maps(self,dpar)
    ! We want to run everything in uK_RJ (at least for compsep), yeah?
    implicit none
    class(dang_data)                 :: self
    type(dang_params), intent(inout) :: dpar
    integer(i4b)                :: j
    
    do j = 1, nbands
       if (.not. dpar%bp_map(j)) then
          if (trim(dpar%band_unit(j)) == 'uK_RJ') then
             self%conversion(j)  = 1.0
          else if (trim(dpar%band_unit(j)) == 'uK_cmb') then
             ! uK_cmb -> uK_RJ
             write(*,*) 'Putting band ', trim(dpar%band_label(j)), ' from uK_cmb to uK_RJ.'
             self%conversion(j)  = 1.0/a2t(bp(j))
          else if (trim(dpar%band_unit(j)) == 'MJy/sr') then
             ! MJy/sr -> uK_RJ
             write(*,*) 'Putting band ', trim(dpar%band_label(j)), ' from MJy/sr to uK_RJ'
             self%conversion(j)  = 1.0/a2f(bp(j))
          else
             write(*,*) 'Not a unit, dumbass! '//dpar%band_unit(j)
             stop
          end if
          self%sig_map(:,:,j) = self%sig_map(:,:,j)*self%conversion(j)
          self%rms_map(:,:,j) = self%rms_map(:,:,j)*self%conversion(j)
       end if
    end do

  end subroutine convert_maps
  
  subroutine convert_bp_maps(self,dpar)
    implicit none
    type(dang_data),   intent(inout) :: self
    type(dang_params), intent(inout) :: dpar
    integer(i4b)                     :: j
    
    do j = 1, nbands
       if (dpar%bp_map(j)) then
          if (trim(dpar%band_unit(j)) == 'uK_RJ') then
             self%conversion(j)  = 1.0
          else if (trim(dpar%band_unit(j)) == 'uK_cmb') then
             ! uK_cmb -> uK_RJ
             write(*,*) 'Putting band ', trim(dpar%band_label(j)), ' from uK_cmb to uK_RJ.'
             self%conversion(j)  = 1.0/a2t(bp(j))
          else if (trim(dpar%band_unit(j)) == 'MJy/sr') then
             ! MJy/sr -> uK_RJ
             write(*,*) 'Putting band ', trim(dpar%band_label(j)), ' from MJy/sr to uK_RJ'
             self%conversion(j)  = 1.0/a2f(bp(j))
          else
             write(*,*) 'Not a unit, dumbass! '//dpar%band_unit(j)
             stop
          end if
          self%sig_map(:,:,j) = self%sig_map(:,:,j)*self%conversion(j)
          self%rms_map(:,:,j) = self%rms_map(:,:,j)*self%conversion(j)
       end if
    end do
    
  end subroutine convert_bp_maps

  subroutine compute_chisq(self,dpar)
    use healpix_types
    implicit none
    type(dang_data),                              intent(inout) :: self
    type(dang_params)                                           :: dpar
    real(dp)                                                    :: s, signal
    integer(i4b)                                                :: i, j, k

    self%chisq = 0.d0
    do i = 0, npix-1
       if (self%masks(i,1) == missval .or. self%masks(i,1) == 0.d0) cycle
       do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
          if (k == 1) then
             do j = 1, nbands
                self%chisq = self%chisq + ((self%sig_map(i,k,j)/self%gain(j)-self%offset(j)) - & 
                     self%sky_model(i,k,j))**2.d0/(self%rms_map(i,k,j)**2.d0)
             end do
          else
             do j = 1, nbands
                self%chisq = self%chisq + (self%sig_map(i,k,j) - self%sky_model(i,k,j))**2.d0 / &
                     & (self%rms_map(i,k,j)**2.d0)
             end do
          end if
       end do
    end do
    self%chisq = self%chisq
    
  end subroutine compute_chisq

  subroutine write_stats_to_term(self,dpar,iter)
    implicit none

    type(dang_data)                     :: self
    type(dang_params)                   :: dpar
    type(dang_comps),        pointer    :: c
    integer(i4b),            intent(in) :: iter
    integer(i4b)                        :: i, j, k

    write(*,fmt='(a)') '---------------------------------------------'
    ! do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
    !    if (rank == master) then
    !       if (mod(iter, 1) == 0 .or. iter == 1) then
    !          write(*,fmt='(i6, a, a, f7.3, a, f8.4, a, 10e10.3)')&
    !               iter, " - Poltype: "//trim(tqu(k)), " - A_s: ",&
    !               component_list(1)%p%amplitude(23000,k),  " - beta_s: ",&
    !               mask_avg(component_list(1)%p%indices(:,k,1),self%masks(:,1)), ' - A_d: ', &
    !               component_list(2)%p%template_amplitudes(:,k)
    !          write(*,fmt='(a)') '---------------------------------------------'
    !       end if
    !    end if
    ! end do
    call compute_chisq(self,dpar)
    if (rank == master) then
       if (mod(iter, 1) == 0 .or. iter == 1) then
          write(*,fmt='(i6,a,E16.5)') iter, " - Chisq: ", self%chisq
          do i = 1, ncomp
             c => component_list(i)%p
             do j = 1, c%nindices
                if (c%sample_index(j)) then
                   write(*,fmt='(a,a,a,a,a,e12.5)')  '     ',trim(c%label), ' ', trim(c%ind_label(j)), ' mean:   ',&
                        mask_avg(c%indices(:,1,j),self%masks(:,1))
                end if
             end do
          end do
          write(*,fmt='(a)') '---------------------------------------------'
       end if
    end if

  end subroutine write_stats_to_term

  ! Data output routines
  subroutine write_maps(dpar,ddata)
    implicit none
    
    type(dang_params)                   :: dpar
    type(dang_data)                     :: ddata
    type(dang_comps),         pointer   :: c
    real(dp), dimension(0:npix-1,nmaps) :: map
    real(dp)                            :: s, signal
    integer(i4b)                        :: n, mn, i, j, k, l
    character(len=128)                  :: title

    logical(lgt)                        :: output

    write(*,*) 'Output data maps'
    
    write(iter_str, '(i0.5)') iter
    ! If we ask to output all components for each band (in band units):
    if (dpar%output_fg .eqv. .true.) then
       do j = 1, nbands
          do n = 1, ncomp
             c => component_list(n)%p
             title = trim(dpar%outdir) // trim(dpar%band_label(j)) //'_'// trim(c%label) //&
                  '_k' // trim(iter_str) // '.fits'
             do i = 0, npix-1
                do k = 1, nmaps
                   map(i,k) = c%eval_signal(j,i,k)/ddata%conversion(j)
                end do
             end do
             ! Mask it!
             do i = 0, npix-1
                if (ddata%masks(i,1) == 0.d0 .or. ddata%masks(i,1) == missval) then
                   map(i,:) = missval
                end if
             end do
             call write_result_map(trim(title), nside, ordering, header, map)
          end do
          title = trim(dpar%outdir) // trim(dpar%band_label(j)) // '_sky_model_k' // trim(iter_str) // '.fits'
          map(:,:)   = ddata%sky_model(:,:,j)/ddata%conversion(j)
          do i = 0, npix-1
             if (ddata%masks(i,1) == 0.d0 .or. ddata%masks(i,1) == missval) then
                map(i,:) = missval
             end if
          end do
          call write_result_map(trim(title), nside, ordering, header, map)
       end do
    end if

    ! Write residual and sky model for each band (output in native band units)
    do j = 1, nbands
       title = trim(dpar%outdir) // trim(dpar%band_label(j)) // '_residual_k' // trim(iter_str) // '.fits'
       map(:,:)   = ddata%res_map(:,:,j)/ddata%conversion(j)
       do i = 0, npix-1
          if (ddata%masks(i,1) == 0.d0 .or. ddata%masks(i,1) == missval) then
             map(i,:) = missval
          end if
       end do
       call write_result_map(trim(title), nside, ordering, header, map)
    end do

    ! Write component maps (all output in uK_RJ)
    do n = 1, ncomp
       output = .true.
       c => component_list(n)%p
       title = trim(dpar%outdir) // trim(c%label) // '_c001_k' // trim(iter_str) // '.fits'
       map(:,:)   = c%amplitude
       do i = 0, npix-1
          if (ddata%masks(i,1) == 0.d0 .or. ddata%masks(i,1) == missval) then
             map(i,:) = missval
          end if
       end do
       call write_result_map(trim(title), nside, ordering, header, map)
       do l = 1, c%nindices
          title = trim(dpar%outdir) // trim(c%label) //&
               '_' // trim(c%ind_label(l))//'_k' // trim(iter_str) // '.fits'
          do i = 0, npix-1
             if (ddata%masks(i,1) == 0.d0 .or. ddata%masks(i,1) == missval) then
                map(i,:) = missval
                cycle
             end if
             do k = 1, nmaps
                map(i,k) = c%indices(i,k,l)
             end do
          end do
          call write_result_map(trim(title),nside,ordering,header,map)
       end do
    end do
    ! Write the chisquare map
    do mn = 1, nmaps
       ddata%chi_map(:,mn) = 0.d0
       ! For intensity
       if (mn == 1) then
          do i = 0, npix-1
             do j = 1, nbands
                ddata%chi_map(i,mn) = ddata%chi_map(i,mn) + ddata%masks(i,1)*(ddata%res_map(i,mn,j)**2)/ddata%rms_map(i,mn,j)**2.d0
             end do
          end do
          
       else
          do i = 0, npix-1
             do j = 1, nbands
                ddata%chi_map(i,mn) = ddata%chi_map(i,mn) + ddata%masks(i,1)*(ddata%sig_map(i,mn,j) - ddata%sky_model(i,mn,j))**2.d0&
                     & /ddata%rms_map(i,mn,j)**2.d0
             end do
          end do
       end if
    end do
    ddata%chi_map(:,:) = ddata%chi_map(:,:)/(nbands)
    title = trim(dpar%outdir) // 'chisq_k'// trim(iter_str) // '.fits'
    map(:,:)   = ddata%chi_map(:,:)
    do i = 0, npix-1
       if (ddata%masks(i,1) == 0.d0 .or. ddata%masks(i,1) == missval) then
          map(i,:) = missval
          ddata%chi_map(i,:) = missval
       end if
    end do
    call write_result_map(trim(title), nside, ordering, header, map)
    
  end subroutine write_maps
  
  subroutine write_data(dpar,ddata,map_n)
    implicit none
    type(dang_params)                   :: dpar
    type(dang_data)                     :: ddata
    type(dang_comps),         pointer   :: c
    integer(i4b),            intent(in) :: map_n
    integer(i4b)                        :: i, n
    character(len=2)                    :: temp_n
    character(len=128)                  :: title, fmt, str_fmt
    character(len=4)                    :: nband_str
    character(len=5)                    :: iter_str

    write(*,*) 'Output data files'

    ! Select output formatting
    write(nband_str, '(i4)') nbands
    write(iter_str, '(i0.5)') iter
    
    ! Output total chisquare
    title = trim(dpar%outdir) // 'total_chisq_' // trim(tqu(map_n)) // '.dat'
    inquire(file=title,exist=exist)
    if (exist) then
       open(33,file=title, status="old",position="append", action="write")
    else
       open(33,file=title, status="new", action="write")
    endif
    call compute_chisq(ddata,dpar)
    write(33,*) ddata%chisq
    close(33)
    
    ! Output template amplitudes - if applicable
    fmt = '('//trim(nband_str)//'(E17.8))'
    do n = 1, ncomp
       c => component_list(n)%p
       if (trim(c%type) == 'template' .or. trim(c%type) == 'hi_fit') then
          title = trim(dpar%outdir) //  trim(c%label) // '_' //trim(tqu(map_n)) // '_amplitudes.dat'
          inquire(file=title,exist=exist)
          if (exist) then
             open(34,file=title, status="old", &
                  position="append", action="write")
          else
             open(34,file=title, status="new", action="write")
             write(34,fmt='('//trim(nband_str)//'(A17))') dpar%band_label
          endif
          write(34,fmt=fmt) c%template_amplitudes(:,map_n)
          close(34)
       end if
    end do
    

    ! And finally band calibration values
    fmt = '(a12,f12.8)'
    
    title = trim(dpar%outdir)//'band_gains_k'//iter_str//'.dat'
    inquire(file=title,exist=exist)
    if (exist) then
       open(37,file=title,status="old",position="append",action="write") 
    else
       open(37,file=title,status="new",action="write")
    end if
    do j = 1, nbands
       write(37,fmt=fmt) trim(dpar%band_label(j)), ddata%gain(j)
    end do
    close(37)
    
    title = trim(dpar%outdir)//'band_offsets_k'//iter_str//'.dat'
    inquire(file=title,exist=exist)
    if (exist) then
       open(38,file=title,status="old",position="append",action="write") 
    else
       open(38,file=title,status="new",action="write")
    end if
    do j = 1, nbands
       write(38,fmt=fmt) trim(dpar%band_label(j)), ddata%offset(j)
    end do
    close(38)
    
  end subroutine write_data

end module dang_data_mod
