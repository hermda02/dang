module dang_data_mod
  use healpix_types
  use dang_util_mod
  use dang_param_mod
  use dang_component_mod
  use dang_bp_mod
  implicit none

  type, public                                :: dang_data

     ! Storage
    
    integer(i4b)                              :: npix
    integer(i4b)                              :: nfgs
    integer(i4b)                              :: nmaps
    real(dp)                                  :: chisq

    character(len=32), allocatable, dimension(:) :: label   ! store the band label 
    character(len=32), allocatable, dimension(:) :: unit    ! store the band unit 

    integer(i4b), allocatable, dimension(:)   :: pol_type         ! Points above to map number

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
    logical(lgt), allocatable, dimension(:)   :: cg_map

  contains
    procedure :: initialize_data_module
    procedure :: init_data_maps
    procedure :: read_data_maps
    procedure :: update_sky_model
    procedure :: convert_maps
    procedure :: mask_hi_threshold
    procedure :: read_band_offsets
    procedure :: read_band_gains
    procedure :: swap_bp_maps
    procedure :: convert_bp_maps
    procedure :: write_data
    procedure :: write_maps
  end type dang_data

  private :: i, j, k, l
  integer(i4b) :: i, j, k, l

contains

  subroutine initialize_data_module(self,dpar)
    implicit none
    class(dang_data),       intent(inout) :: self
    type(dang_params)                     :: dpar

    call self%init_data_maps(dpar)
    call self%read_data_maps(dpar)
    call self%read_band_offsets(dpar)
    call self%read_band_gains(dpar)
    if (dpar%bp_swap) then
       call self%swap_bp_maps(dpar)
       write(*,*) ''
       call self%convert_bp_maps
       write(*,*) ''
    end if
    call self%convert_maps
    
  end subroutine initialize_data_module

  subroutine init_data_maps(self,dpar)
   ! Allocate and intialize arrays - pass dpar info into ddata objects
    implicit none
    class(dang_data),       intent(inout) :: self
    type(dang_params),      intent(in)    :: dpar

    integer(i4b)                          :: i

    write(*,*) "Initializing data maps"
    
    allocate(self%masks(0:npix-1,nmaps))

    allocate(self%label(nbands))
    allocate(self%unit(nbands))
    allocate(self%cg_map(nbands))

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

    allocate(self%pol_type(size(dpar%pol_type)))

    self%gain   = 1.d0
    self%offset = 0.d0

    self%fit_gain = dpar%fit_gain
    self%fit_offset = dpar%fit_offs

    self%pol_type = dpar%pol_type

    self%label  = dpar%band_label
    self%unit   = dpar%band_unit
    self%cg_map = dpar%bp_map

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
       do k = 1, nmaps
          if (self%masks(i,k) == 0.d0 .or. self%masks(i,k) == missval) then
             self%masks(i,k) = 0.d0
          else 
             nump = nump + 1
          end if
       end do
    end do
    
    ! Read maps in
    do j = 1, nbands
       ! If not supposed to be swapped BP maps, load that map
       if (.not. self%cg_map(j)) then
          call read_bintab(trim(dpar%datadir) // trim(dpar%band_noisefile(j)), &
               rms,npix,nmaps,nullval,anynull,header=header)
          self%rms_map(:,:,j) = rms
          call read_bintab(trim(dpar%datadir) // trim(dpar%band_mapfile(j)), &
               map,npix,nmaps,nullval,anynull,header=header)
          self%sig_map(:,:,j) = map
       end if
    end do
    
    deallocate(map,rms)
  end subroutine read_data_maps

  subroutine swap_bp_maps(self,dpar)
    type(dang_params)                               :: dpar
    class(dang_data),                 intent(inout) :: self
    character(len=512)                              :: chain_c
    character(len=300), allocatable, dimension(:,:) :: bp_maps
    integer(i4b)                                    :: i, j, iter_i, chain_i
    character(len=6)                                :: iter_str
    real(dp), allocatable, dimension(:,:)           :: map, rms
    real(dp)                                        :: norm
    double precision                                :: temp(1)

    allocate(map(0:npix-1,3))
    allocate(rms(0:npix-1,3))
    allocate(bp_maps(dpar%numband,2))

    do j = 1, dpar%numband
       if (self%cg_map(j)) then

          call RANDOM_SEED()
          call RANDOM_NUMBER(temp)
          
          norm    = temp(1)*dpar%num_chains
          chain_i = int(norm)+1
          chain_c = dpar%bp_chain_list(chain_i)
          
          call RANDOM_SEED()
          call RANDOM_NUMBER(temp)
          
          norm    = temp(1)*(dpar%bp_max-dpar%bp_burnin)
          iter_i  = int(norm)+1+dpar%bp_burnin
          
          write(iter_str,'(i0.6)') iter_i

          ! Normal BP switching here
          bp_maps(j,1) = trim(dpar%bp_dir) // trim(self%label(j))//'_map_'//trim(chain_c)//&
               '_n0064_60arcmin_k'//trim(iter_str) // '.fits'

          ! This is for switching to the 'plus_cmb' maps:
          bp_maps(j,2) = trim(dpar%bp_dir) // trim(self%label(j))//'_rms_'//trim(chain_c)//&
               '_n0064_60arcmin_k'//trim(iter_str) // '.fits'
          write(*,'(a,a,a)') 'Swapping band ', trim(self%label(j)), '.'
          call read_bintab(trim(bp_maps(j,1)),map,self%npix,3,nullval,anynull,header=header)
          self%sig_map(:,:,j) = map
          call read_bintab(trim(bp_maps(j,2)),rms,self%npix,3,nullval,anynull,header=header)
          self%rms_map(:,:,j) = rms
       end if
    end do

  end subroutine swap_bp_maps

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

    loaded(:) = .false.

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
             if (trim(band) == trim(self%label(j))) then
                self%offset(j) = offset
                loaded(j) = .true.
             end if
          end do
       end do
    end if
    do j = 1, nbands
       if (.not. loaded(j)) then
          write(*,*) trim(self%label(j))//' offset not loaded -- set to 0'
       end if
    end do
    write(*,*) ''

  end subroutine read_band_offsets

  subroutine read_band_gains(self,dpar)
    implicit none
    class(dang_data),       intent(inout) :: self
    type(dang_params),      intent(in)    :: dpar
    character(len=512)                    :: file
    character(len=128)                    :: fmt, nband_str, band

    real(dp)                              :: gain

    logical(lgt)                          :: exist
    logical(lgt), allocatable, dimension(:) :: loaded

    integer(i4b)                          :: i, j, k
    integer(i4b)                          :: unit, ios, ierror

    allocate(loaded(nbands))

    loaded(:) = .false.

    write(nband_str, '(i4)') nbands

    file = trim(dpar%gain_file)

    unit = getlun()
    ierror  = 0
    fmt  = '('//trim(nband_str)//'(E17.8))'

    if (trim(file) == '') then
       write(*,*) 'No BAND_GAIN_FILE -- setting all gains to 1.0'
       self%gain(:) = 1.d0
    else
       open(unit,file=trim(dpar%datadir)//file)
       do while (ierror .eq. 0) 
          read(unit=unit,fmt=*,iostat=ierror) band, gain
          do j = 1, nbands
             if (trim(band) == trim(self%label(j))) then
                self%gain(j) = gain
                loaded(j) = .true.
             end if
          end do
       end do
    end if
    do j = 1, nbands
       if (.not. loaded(j)) then
          write(*,*) trim(self%label(j))//' gain not loaded -- set to 1.0'
       end if
    end do
    write(*,*) ''
    
  end subroutine read_band_gains

  subroutine update_sky_model(self)
    implicit none
    class(dang_data),  intent(inout) :: self
    type(dang_comps),  pointer       :: c
    integer(i4b)                     :: i, j, k, l

    self%sky_model(:,:,:) = 0.d0
    do l = 1, ncomp
       c => component_list(l)%p
       if (c%type == 'monopole') then
          self%offset = c%template_amplitudes(:,1)
          cycle
       end if
       !$OMP PARALLEL PRIVATE(i,j,k)
       !$OMP DO SCHEDULE(static)
       do i = 0, npix-1
          do k = 1, nmaps
             do j = 1, nbands
                self%sky_model(i,k,j) = self%sky_model(i,k,j) + c%eval_signal(j,i,k)
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    end do
    !$OMP BARRIER

    ! Calculate the residual
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(static)
    do i = 0, npix-1
       do k = 1, nmaps
          do j = 1, nbands
             if (k == 1) then
                ! For intensity
                self%res_map(i,1,j) = (self%sig_map(i,1,j)-self%offset(j))/self%gain(j)-self%sky_model(i,1,j)
             else
                ! and polarization
                self%res_map(i,k,j) = self%sig_map(i,k,j)-self%sky_model(i,k,j)
             end if
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP BARRIER

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
       end if
    end do
    do i = 0, npix-1
       if (c%template(i,1) > dpar%thresh) then
          self%masks(i,1) = 0.d0
       else if (self%masks(i,1) == missval .or. self%masks(i,1) == 0.d0) then
          self%masks(i,1) = 0.d0
       else if (self%rms_map(i,1,1) == 0.d0) then
          self%masks(i,1) = 0.d0
       else
          self%masks(i,1) = 1.d0
       end if
    end do
  end subroutine mask_hi_threshold

  subroutine convert_maps(self)!,dpar)
    ! We want to run everything in uK_RJ (at least for compsep), yeah?
    implicit none
    class(dang_data)                 :: self
   !  type(dang_params), intent(inout) :: dpar
    integer(i4b)                :: j
    
    do j = 1, nbands
       if (.not. self%cg_map(j)) then
          if (trim(self%unit(j)) == 'uK_RJ') then
             self%conversion(j)  = 1.0
          else if (trim(self%unit(j)) == 'uK_cmb') then
             ! uK_cmb -> uK_RJ
             write(*,*) 'Putting band ', trim(self%label(j)), ' from uK_cmb to uK_RJ.'
             self%conversion(j)  = 1.0/a2t(bp(j))
          else if (trim(self%unit(j)) == 'MJy/sr') then
             ! MJy/sr -> uK_RJ
             write(*,*) 'Putting band ', trim(self%label(j)), ' from MJy/sr to uK_RJ.'
             self%conversion(j)  = 1.0/a2f(bp(j))
          else
             write(*,*) 'Not a unit, dumbass! '//self%unit(j)
             stop
          end if
          self%sig_map(:,:,j) = self%sig_map(:,:,j)*self%conversion(j)
          self%rms_map(:,:,j) = self%rms_map(:,:,j)*self%conversion(j)
          self%offset(j)      = self%offset(j)*self%conversion(j)
          ! Set the loaded monopole values into the monopole component
         do i = 1, size(component_list)
            if (trim(component_list(i)%p%type) == 'monopole') then
               component_list(i)%p%template_amplitudes(:,1) =self%offset
            end if
         end do
       end if
    end do

  end subroutine convert_maps
  
  subroutine convert_bp_maps(self)
    implicit none
    class(dang_data),   intent(inout) :: self
    integer(i4b)                     :: j
    
    do j = 1, nbands
       if (self%cg_map(j)) then
          if (trim(self%unit(j)) == 'uK_RJ') then
             self%conversion(j)  = 1.0
          else if (trim(self%unit(j)) == 'uK_cmb') then
             ! uK_cmb -> uK_RJ
             write(*,*) 'Putting band ', trim(self%label(j)), ' from uK_cmb to uK_RJ.'
             self%conversion(j)  = 1.0/a2t(bp(j))
          else if (trim(self%unit(j)) == 'MJy/sr') then
             ! MJy/sr -> uK_RJ
             write(*,*) 'Putting band ', trim(self%label(j)), ' from MJy/sr to uK_RJ'
             self%conversion(j)  = 1.0/a2f(bp(j))
          else
             write(*,*) 'Not a unit, dumbass! '//self%unit(j)
             stop
          end if
          self%sig_map(:,:,j) = self%sig_map(:,:,j)*self%conversion(j)
          self%rms_map(:,:,j) = self%rms_map(:,:,j)*self%conversion(j)
          self%offset(j)      = self%offset(j)*self%conversion(j)
       end if
    end do
    
  end subroutine convert_bp_maps

  subroutine compute_chisq(self)
    use healpix_types
    implicit none
    type(dang_data),                              intent(inout) :: self
   !  type(dang_params)                                           :: dpar
    real(dp)                                                    :: s, signal, chisq
    integer(i4b)                                                :: i, j, k

    self%chisq = 0.d0
    self%chi_map(:,:) = 0.d0
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(static)
    do i = 0, npix-1
       if (self%masks(i,1) == missval .or. self%masks(i,1) == 0.d0) cycle
       do k = self%pol_type(1), self%pol_type(size(self%pol_type))
          if (k == 1) then
             do j = 1, nbands
                self%chi_map(i,k) = self%chi_map(i,k) + ((self%sig_map(i,k,j)-self%offset(j))/self%gain(j) - & 
                     self%sky_model(i,k,j))**2.d0/(self%rms_map(i,k,j)**2.d0)
             end do
          else
             do j = 1, nbands
                self%chi_map(i,k) = self%chi_map(i,k) + (self%sig_map(i,k,j) - self%sky_model(i,k,j))**2.d0 / &
                     & (self%rms_map(i,k,j)**2.d0)
             end do
          end if
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    self%chisq = self%chisq + sum(self%chi_map)/nump
    self%chi_map(:,:) = self%chi_map(:,:)/nbands
    
  end subroutine compute_chisq

  subroutine write_stats_to_term(self,iter)
    implicit none

    type(dang_data)                     :: self
    type(dang_comps),        pointer    :: c
    integer(i4b),            intent(in) :: iter
    integer(i4b)                        :: i, j, k

    write(*,fmt='(a)') '---------------------------------------------'
    call compute_chisq(self)
    if (rank == master) then
       if (mod(iter, 1) == 0 .or. iter == 1) then
          write(*,fmt='(i6,a,E16.5)') iter, " - Chisq: ", self%chisq
          do i = 1, ncomp
             c => component_list(i)%p
             do j = 1, c%nindices
                if (c%sample_index(j)) then
                   do k = 1, c%nflag(j)
                      if (iand(c%pol_flag(j,k),1) .ne. 0) then
                         write(*,fmt='(a,a,a,a,a,f12.5)')  '     ',trim(c%label), ' ', trim(c%ind_label(j)), ' I mean:   ',&
                              mask_avg(c%indices(:,1,j),self%masks(:,1))
                      else if (iand(c%pol_flag(j,k),2) .ne. 0) then
                         write(*,fmt='(a,a,a,a,a,f12.5)')  '     ',trim(c%label), ' ', trim(c%ind_label(j)), ' Q mean:   ',&
                              mask_avg(c%indices(:,2,j),self%masks(:,1))
                      else if (iand(c%pol_flag(j,k),4) .ne. 0) then
                         write(*,fmt='(a,a,a,a,a,f12.5)')  '     ',trim(c%label), ' ', trim(c%ind_label(j)), ' U mean:   ',&
                              mask_avg(c%indices(:,3,j),self%masks(:,1))
                      else if (iand(c%pol_flag(j,k),8) .ne. 0) then
                         write(*,fmt='(a,a,a,a,a,f12.5)')  '     ',trim(c%label), ' ', trim(c%ind_label(j)), ' Q+U mean:   ',&
                              mask_avg(c%indices(:,2,j),self%masks(:,1))
                      else if (iand(c%pol_flag(j,k),0) .ne. 0) then
                         write(*,fmt='(a,a,a,a,a,f12.5)')  '     ',trim(c%label), ' ', trim(c%ind_label(j)), ' I+Q+U mean:   ',&
                              mask_avg(c%indices(:,1,j),self%masks(:,1))
                      end if
                   end do
                end if
             end do
          end do
          write(*,fmt='(a)') '---------------------------------------------'
       end if
    end if

  end subroutine write_stats_to_term

  ! Data output routines
  subroutine write_maps(self,dpar,suff)
    implicit none
    
    type(dang_params)                   :: dpar
    class(dang_data)                    :: self
    type(dang_comps),         pointer   :: c
    character(len=*), intent(in), optional :: suff
    real(dp), dimension(0:npix-1,nmaps) :: map
    real(dp)                            :: s, signal
    integer(i4b)                        :: n, mn, i, j, k, l
    character(len=128)                  :: title, suffix

    logical(lgt)                        :: output

    write(*,*) 'Output data maps'

    suffix = '.fits'
    if (present(suff)) then
       suffix = '_' // trim(suff) // trim(suffix)
    end if
    
    write(iter_str, '(i0.5)') iter
    ! If we ask to output all components for each band (in band units):
    if (dpar%output_fg .eqv. .true.) then
       do j = 1, nbands
          do n = 1, ncomp
             c => component_list(n)%p
             title = trim(dpar%outdir) // trim(self%label(j)) //'_'// trim(c%label) //&
                  '_k' // trim(iter_str) // suffix
             !$OMP PARALLEL PRIVATE(i,k)
             !$OMP DO
             do i = 0, npix-1
                do k = 1, nmaps
                   map(i,k) = c%eval_signal(j,i,k)/self%conversion(j)
                end do
             end do
             !$OMP END DO
             !$OMP END PARALLEL
             !$OMP BARRIER
             ! Mask it!
             call apply_dang_mask(map,self%masks(:,1),missing=.true.)
             call write_result_map(trim(title), nside, ordering, header, map)
          end do
          title = trim(dpar%outdir) // trim(self%label(j)) // '_sky_model_k' // trim(iter_str) // suffix
          map(:,:)   = self%sky_model(:,:,j)/self%conversion(j)
          !$OMP PARALLEL PRIVATE(i)
          !$OMP DO
          do i = 0, npix-1
             if (self%masks(i,1) == 0.d0 .or. self%masks(i,1) == missval) then
                map(i,:) = missval
             end if
          end do
          !$OMP END DO
          !$OMP END PARALLEL
          call write_result_map(trim(title), nside, ordering, header, map)
       end do
    end if

    ! Write residual and sky model for each band (output in native band units)
    do j = 1, nbands
       title = trim(dpar%outdir) // trim(self%label(j)) // '_residual_k' // trim(iter_str) // suffix
       map(:,:)   = self%res_map(:,:,j)/self%conversion(j)
       ! Mask it!
       call apply_dang_mask(map,self%masks(:,1),missing=.true.)
       call write_result_map(trim(title), nside, ordering, header, map)
    end do

    ! Write component maps (all output in uK_RJ)
    do n = 1, ncomp
       output = .true.
       c => component_list(n)%p
       title = trim(dpar%outdir) // trim(c%label) // '_c001_k' // trim(iter_str) // suffix
       map(:,:)   = c%amplitude
       ! Mask it!
       call apply_dang_mask(map,self%masks(:,1),missing=.true.)
       call write_result_map(trim(title), nside, ordering, header, map)

       do l = 1, c%nindices
          title = trim(dpar%outdir) // trim(c%label) //&
               '_' // trim(c%ind_label(l))//'_k' // trim(iter_str) // suffix
          map(:,:) = c%indices(:,:,l)
          ! Mask it!
          call apply_dang_mask(map,self%masks(:,1),missing=.true.)
          call write_result_map(trim(title),nside,ordering,header,map)
       end do
    end do
    ! Write the chisquare map
    call compute_chisq(self)
    title = trim(dpar%outdir) // 'chisq_k'// trim(iter_str) // suffix
    map(:,:)   = self%chi_map(:,:)
    ! Mask it!
    call apply_dang_mask(map,self%masks(:,1),missing=.true.)
    call write_result_map(trim(title), nside, ordering, header, map)
    
  end subroutine write_maps
  
  subroutine write_data(self,dpar,map_n)
    implicit none
    type(dang_params)                   :: dpar
    class(dang_data)                    :: self
    type(dang_comps),         pointer   :: c
    integer(i4b),            intent(in) :: map_n
    integer(i4b)                        :: i, j, n, unit
    character(len=2)                    :: temp_n
    character(len=128)                  :: title, fmt, str_fmt
    character(len=1)                    :: nmaps_str
    character(len=4)                    :: nband_str
    character(len=5)                    :: iter_str

    write(*,*) 'Output data files'

    ! Select output formatting
    write(nband_str, '(i4)') nbands
    write(iter_str, '(i0.5)') iter
    write(nmaps_str, '(i1)') nmaps
    
    ! Output total chisquare
    title = trim(dpar%outdir) // 'total_chisq_' // trim(tqu(map_n)) // '.dat'
    inquire(file=title,exist=exist)
    if (exist) then
       open(33,file=title, status="old",position="append", action="write")
    else
       open(33,file=title, status="new", action="write")
    endif
    call compute_chisq(self)
    write(33,*) self%chisq
    close(33)

    ! Output template amplitudes - if applicable
    do n = 1, ncomp
       c => component_list(n)%p
       if (trim(c%type) == 'template' .or. trim(c%type) == 'hi_fit') then
          unit = getlun()
          title = trim(dpar%outdir) //  trim(c%label) // '_' //trim(tqu(map_n)) // '_amplitudes.dat'
          inquire(file=title,exist=exist)
          if (exist) then
             open(unit,file=title, status="old", &
                  position="append", action="write")
          else
             open(unit,file=title, status="new", action="write")
             write(unit,fmt='('//trim(nband_str)//'(A17))') dpar%band_label
          endif
          write(unit,fmt='('//trim(nband_str)//'(E17.8))') c%template_amplitudes(:,map_n)
          close(unit)
       end if
       
       ! Output mean spectral index values
       do j = 1, c%nindices
          if (c%sample_index(j)) then
             fmt = '('//nmaps_str//'(f12.8))'
             unit = getlun()
             title = trim(dpar%outdir) // trim(c%label) //'_' // trim(c%ind_label(j))//'_mean_'// trim(tqu(map_n)) // '.dat'
             inquire(file=title,exist=exist)
             if (exist) then
                open(unit,file=title,status="old",position="append",action="write")
             else
                open(unit,file=title,status="new",action="write")
             end if
             write(unit,fmt=fmt) mask_avg(c%indices(:,map_n,j),self%masks(:,1))
          end if
       end do
    end do

    ! And finally band calibration values
    fmt = '(a12,E16.8)'
    
    title = trim(dpar%outdir)//'band_gains_k'//iter_str//'.dat'
    inquire(file=title,exist=exist)
    if (exist) then
       open(37,file=title,status="old",position="append",action="write") 
    else
       open(37,file=title,status="new",action="write")
    end if
    do j = 1, nbands
       write(37,fmt=fmt) trim(self%label(j)), self%gain(j)
    end do
    close(37)
    
    title = trim(dpar%outdir)//'band_offsets_k'//iter_str//'.dat'
    inquire(file=title,exist=exist)
    if (exist) then
       open(38,file=title,status="old",position="append",action="write") 
    else
       open(38,file=title,status="new",action="write")
    end if
    do n = 1, ncomp
       c => component_list(n)%p
       if (trim(c%type) == 'monopole') then
          do j = 1, nbands
             write(38,fmt=fmt) trim(self%label(j)), c%template_amplitudes(j,1)/self%conversion(j)
          end do
       end if
    end do
    close(38)

  end subroutine write_data

end module dang_data_mod
