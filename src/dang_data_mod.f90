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
    procedure :: mask_hi
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

  ! subroutine init_band_calibration(self)
  !   !==================================================|
  !   ! Write a routine that tells ddata which component | 
  !   ! we want to calibrate that band to.               |
  !   !==================================================|
  !   ! Right now we leave this out and calibrate via    |
  !   ! broadband (i.e. full sky model).
  !   implicit none
  !   class(dang_data),       intent(inout) :: self

  !   type(dang_comps),   pointer             :: c2
        


  ! end subroutine init_band_calibration

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
       if (trim(dpar%mode) == 'comp_sep') then
          if (.not. dpar%bp_map(j)) then
             call read_bintab(trim(dpar%datadir) // trim(dpar%band_noisefile(j)), &
                  rms,npix,nmaps,nullval,anynull,header=header)
             self%rms_map(:,:,j) = rms
             call read_bintab(trim(dpar%datadir) // trim(dpar%band_mapfile(j)), &
                  map,npix,nmaps,nullval,anynull,header=header)
             self%sig_map(:,:,j) = map
          end if
       else 
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
    if (trim(dpar%mode) == 'comp_sep') then
       do i = 1, dpar%ntemp
          call read_bintab(dpar%temp_file(i), self%temps(:,:,i), npix,nmaps,nullval,anynull,header=header)
          ! Normalize templates
          do k = 1, nmaps
             self%temp_norm(k,i) = 1.d0!maxval(self%temps(:,k,i))
             self%temps(:,k,i)   = self%temps(:,k,i)/self%temp_norm(k,i)
          end do
       end do
    end if
  end subroutine read_data_maps

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

  subroutine mask_hi(self, dpar, c)
    implicit none
    class(dang_data),            intent(inout) :: self
    type(dang_params)                          :: dpar
    type(dang_comps),   pointer, intent(in)    :: c
    integer(i4b)                               :: i, j

    do i = 0, npix-1
       if (c%template(i,1) > dpar%thresh) then
          self%masks(i,1) = missval
       else if (self%masks(i,1) == missval) then
          self%masks(i,1) = missval
       else if (self%rms_map(i,1,1) == 0.d0) then
          self%masks(i,1) = missval
       else
          self%masks(i,1) = 1.d0
       end if
    end do
    nump = 0
    do i = 0, npix-1
       do j = 1, nmaps
          if (self%masks(i,j) == 0.d0 .or. self%masks(i,j) == missval) then
             self%masks(i,j) = missval
          else 
             nump = nump + 1
          end if
       end do
    end do
  end subroutine mask_hi
  
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
    
    if (dpar%mode == 'comp_sep') then

       do j = 1, nbands
          if (.not. dpar%bp_map(j)) then
             if (trim(dpar%band_unit(j)) == 'uK_RJ') then
                cycle
             else if (trim(dpar%band_unit(j)) == 'uK_cmb') then
                ! uK_cmb -> uK_RJ
                write(*,*) 'Putting band ', trim(dpar%band_label(j)), ' from uK_cmb to uK_RJ.'
                self%sig_map(:,:,j) = self%sig_map(:,:,j)/a2t(bp(j))
                self%rms_map(:,:,j) = self%rms_map(:,:,j)/a2t(bp(j))
             else if (trim(dpar%band_unit(j)) == 'MJy/sr') then
                ! MJy/sr -> uK_RJ
                write(*,*) 'Putting band ', trim(dpar%band_label(j)), ' from MJy/sr to uK_RJ'
                self%sig_map(:,:,j) = self%sig_map(:,:,j)/a2f(bp(j))
                self%rms_map(:,:,j) = self%rms_map(:,:,j)/a2f(bp(j))
             else
                write(*,*) 'Not a unit, dumbass! '//dpar%band_unit(j)
                stop
             end if
          end if
       end do
       
    else if (dpar%mode == 'hi_fit') then

       do j = 1, nbands
          if (.not. dpar%bp_map(j)) then
             if (trim(dpar%band_unit(j)) == 'MJy/sr') then
                cycle
             else if (trim(dpar%band_unit(j)) == 'uK_cmb') then
                ! uK_cmb -> MJy/sr
                ! Check bandpass type
                write(*,*) 'Putting band ', trim(dpar%band_label(j)), ' from uK_cmb to MJy/sr.'
                self%sig_map(:,:,j) = self%sig_map(:,:,j)/f2t(bp(j))
                self%rms_map(:,:,j) = self%rms_map(:,:,j)/f2t(bp(j))
             else if (trim(dpar%band_unit(j)) == 'MJy/sr') then
                ! uK_RJ -> MJy/sr
                write(*,*) 'Putting band ', trim(dpar%band_label(j)), ' from uK_RJ to MJy/sr'
                self%sig_map(:,:,j) = self%sig_map(:,:,j)*a2f(bp(j))
                self%rms_map(:,:,j) = self%rms_map(:,:,j)*a2f(bp(j))
             else
                write(*,*) 'Not a unit, dumbass! '//dpar%band_unit(j)
                stop
             end if
          end if
       end do
    end if

  end subroutine convert_maps
  
  subroutine convert_bp_maps(self,dpar)
    implicit none
    type(dang_data),   intent(inout) :: self
    type(dang_params), intent(inout) :: dpar
    integer(i4b)                     :: j
    
    do j = 1, nbands
       if (dpar%bp_map(j)) then
          if (trim(dpar%band_unit(j)) == 'uK_RJ') then
             cycle
          else if (trim(dpar%band_unit(j)) == 'uK_cmb') then
             ! uK_cmb -> uK_RJ
             write(*,*) 'Putting band ', trim(dpar%band_label(j)), ' from uK_cmb to uK_RJ.'
             self%sig_map(:,:,j) = self%sig_map(:,:,j)/a2t(bp(j))
             self%rms_map(:,:,j) = self%rms_map(:,:,j)/a2t(bp(j))
          else if (trim(dpar%band_unit(j)) == 'MJy/sr') then
             ! MJy/sr -> uK_RJ
             write(*,*) 'Putting band ', trim(dpar%band_label(j)), ' from MJy/sr to uK_RJ'
             self%sig_map(:,:,j) = self%sig_map(:,:,j)/a2f(bp(j))
             self%rms_map(:,:,j) = self%rms_map(:,:,j)/a2f(bp(j))
          else
             write(*,*) 'Not a unit, dumbass!'
             stop
          end if
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

    if (trim(dpar%mode) == 'comp_sep') then
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
       
    else if (trim(dpar%mode) == 'hi_fit') then
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
       ! if (rank == master) then
       !    if (mod(iter, 1) == 0 .or. iter == 1) then
       !       if (nbands .lt. 10) then
       !          write(*,fmt='(i6, a, f16.5, a, f10.3, a, 10e10.3)')&
       !               iter, " - chisq: " , self%chisq, " - T_d: ",&
       !               mask_avg(component_list(1)%p%indices(:,1,1),self%masks(:,1)),&
       !               ' - A_HI: ', component_list(1)%p%template_amplitudes(:,1)
       !          write(*,fmt='(a)') '---------------------------------------------'
       !       else
       !          write(*,fmt='(i6, a, E10.3, a, e10.3)')&
       !               iter, " - chisq: " , self%chisq, " - T_d: ",&
       !               mask_avg(dcomps%T_d(:,1),self%masks(:,1))
       !          write(*,fmt='(a)') '---------------------------------------------'
       !       end if
       !    end if
       ! end if
    end if

  end subroutine write_stats_to_term

  ! Data output routines
  subroutine write_maps(dpar,dat)
    implicit none
    
    type(dang_params)                   :: dpar
    type(dang_data)                     :: dat
    type(dang_comps),         pointer   :: c
    real(dp), dimension(0:npix-1,nmaps) :: map
    real(dp)                            :: s, signal
    integer(i4b)                        :: n, mn, i, j, k, l
    character(len=128)                  :: title

    logical(lgt)                        :: output

    write(*,*) 'Output data maps'
    
    if (trim(dpar%mode) == 'comp_sep') then
       
       write(iter_str, '(i0.5)') iter
       ! If we ask to output all components for each band:
       if (dpar%output_fg .eqv. .true.) then
          do j = 1, nbands
             do n = 1, ncomp
                c => component_list(n)%p
                title = trim(dpar%outdir) // trim(dpar%band_label(j)) //'_'// trim(c%label) //&
                     '_k' // trim(iter_str) // '.fits'
                if (c%type /= 'template') then
                   do i = 0, npix-1
                      do k = 1, nmaps
                         map(i,k) = c%amplitude(i,k)*c%eval_sed(j,i,k)
                      end do
                   end do
                else
                   do i = 0, npix-1
                      do k = 1, nmaps
                         map(i,k) = c%template(i,k)*c%template_amplitudes(j,k)
                      end do
                   end do
                end if
                do i = 0, npix-1
                   if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                      map(i,:) = missval
                   end if
                end do
                call write_result_map(trim(title), nside, ordering, header, map)
             end do
          end do
       end if
       ! Write residual and sky model for each band
       do j = 1, nbands
          title = trim(dpar%outdir) // trim(dpar%band_label(j)) // '_residual_k' // trim(iter_str) // '.fits'
          map(:,:)   = dat%res_map(:,:,j)
          do i = 0, npix-1
             if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                map(i,:) = missval
             end if
          end do
          call write_result_map(trim(title), nside, ordering, header, map)


          title = trim(dpar%outdir) // trim(dpar%band_label(j)) // '_sky_model_k' // trim(iter_str) // '.fits'
          map(:,:)   = dat%sky_model(:,:,j)
          do i = 0, npix-1
             if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                map(i,:) = missval
             end if
          end do
          call write_result_map(trim(title), nside, ordering, header, map)
       end do
       ! Write component maps
       do n = 1, ncomp
          output = .true.
          c => component_list(n)%p
          title = trim(dpar%outdir) // trim(c%label) // '_c001_k' // trim(iter_str) // '.fits'
          map(:,:)   = c%amplitude
          do i = 0, npix-1
             if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                map(i,:) = missval
             end if
          end do
          call write_result_map(trim(title), nside, ordering, header, map)
          do l = 1, c%nindices
             title = trim(dpar%outdir) // trim(c%label) //&
                  '_' // trim(c%ind_label(l))//'_k' // trim(iter_str) // '.fits'
             do i = 0, npix-1
                do k = 1, nmaps
                   map(i,k) = c%indices(i,k,l)
                end do
             end do
             call write_result_map(trim(title),nside,ordering,header,map)
          end do
       end do
       ! Write the chisquare map
       do mn = 1, nmaps
          dat%chi_map(:,mn) = 0.d0
          ! For intensity
          if (mn == 1) then
             do i = 0, npix-1
                do j = 1, nbands
                   dat%chi_map(i,mn) = dat%chi_map(i,mn) + dat%masks(i,1)*(dat%res_map(i,mn,j)**2)/dat%rms_map(i,mn,j)**2.d0
                end do
             end do

          else
             do i = 0, npix-1
                do j = 1, nbands
                   dat%chi_map(i,mn) = dat%chi_map(i,mn) + dat%masks(i,1)*(dat%sig_map(i,mn,j) - dat%sky_model(i,mn,j))**2.d0&
                        & /dat%rms_map(i,mn,j)**2.d0
                end do
             end do
          end if
       end do
       dat%chi_map(:,:) = dat%chi_map(:,:)/(nbands)
       title = trim(dpar%outdir) // 'chisq_k'// trim(iter_str) // '.fits'
       map(:,:)   = dat%chi_map(:,:)
       do i = 0, npix-1
          if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
             map(i,:) = missval
             dat%chi_map(i,:) = missval
          end if
       end do
       call write_result_map(trim(title), nside, ordering, header, map)

    ! For hi_fit maps
    ! Commented out this whole section for now
    else if (trim(dpar%mode) == 'hi_fit') then
       
       write(iter_str, '(i0.5)') iter
       do j = 1, nbands
          title = trim(dpar%outdir) // trim(dpar%band_label(j)) //'_hi_amplitude_k'// trim(iter_str) // '.fits'
          map(:,1)   = dat%sky_model(:,1,j)!component_list(1)%p%template_amplitudes(j,1)*component_list(1)%p%template(:,1)
          do i = 0, npix-1
             if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                map(i,1) = missval
             end if
          end do
          call write_bintab(map,npix,1, header, nlheader, trim(title))
       end do
       ! Write out residual maps
       do j = 1, nbands
          title = trim(dpar%outdir) // trim(dpar%band_label(j)) // '_residual_k' // trim(iter_str) // '.fits'
          map(:,:)   = dat%res_map(:,:,j)
          do i = 0, npix-1
             if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                map(i,1) = missval
             end if
          end do
          call write_bintab(map,npix,1, header, nlheader, trim(title))
       end do

       ! Write out T_d map
       title = trim(dpar%outdir) // 'T_d_k'// trim(iter_str) // '.fits'
       map(:,1)   = component_list(1)%p%indices(:,1,1)
       do i = 0, npix-1
          if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
             map(i,1) = missval
          end if
       end do
       call write_bintab(map,npix,1, header, nlheader, trim(title))

       ! Compute and write out \chi^2 map
       dat%chi_map = 0.d0
       do i = 0, npix-1
          do j = 1, nbands
             s = dat%gain(j)*dat%sky_model(i,1,j)+dat%offset(j)
             dat%chi_map(i,1) = dat%chi_map(i,1) + dat%masks(i,1)*(dat%sig_map(i,1,j) - s)**2.d0/dat%rms_map(i,1,j)**2.d0
          end do
       end do
       dat%chi_map(:,1) = dat%chi_map(:,1)!/(nbands)
       title = trim(dpar%outdir) // 'chisq_k' // trim(iter_str) // '.fits'
       map(:,1)   = dat%chi_map(:,1)
       do i = 0, npix-1
          if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
             map(i,1) = missval
          end if
       end do
       call write_bintab(map,npix,1, header, nlheader, trim(title))
    end if
    
  end subroutine write_maps
  
  subroutine write_data(dpar,dat,map_n)
    implicit none
    type(dang_params)                   :: dpar
    type(dang_data)                     :: dat
    type(dang_comps),         pointer   :: c
    integer(i4b),            intent(in) :: map_n
    integer(i4b)                        :: i, n
    character(len=2)                    :: temp_n
    character(len=128)                  :: title, fmt
    character(len=4)                    :: nband_str

    write(*,*) 'Output data files'

    if (trim(dpar%mode) == 'comp_sep') then
       
       title = trim(dpar%outdir) // 'total_chisq_' // trim(tqu(map_n)) // '.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(33,file=title, status="old",position="append", action="write")
       else
          open(33,file=title, status="new", action="write")
       endif
       call compute_chisq(dat,dpar)!,dcomps,map_n)
       write(33,*) dat%chisq
       close(33)
       
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
             endif
             write(34,'(10(E17.8))') c%template_amplitudes(:,map_n)
             close(34)
          ! else
          !    do i = 1, c%nindices
          !       title = trim(dpar%outdir) // trim(c%label) // '_' // trim(c%ind_label(i)) // '.dat'
                
          !    end do

          end if
       end do

       write(nband_str, '(i4)') nbands

       fmt = '('//trim(nband_str)//'(E17.8))'

       title = trim(dpar%outdir)//'band_gains.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(37,file=title,status="old",position="append",action="write") 
       else
          open(37,file=title,status="new",action="write")
       end if
       write(37,fmt=fmt) dat%gain
       close(37)

       title = trim(dpar%outdir)//'band_offsets.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(38,file=title,status="old",position="append",action="write") 
       else
          open(38,file=title,status="new",action="write")
       end if
       write(38,fmt=fmt) dat%offset
       close(38)


       
    else if (trim(dpar%mode) == 'hi_fit') then

       write(nband_str, '(i4)') nbands

       fmt = '('//trim(nband_str)//'(E17.8))'
       title = trim(dpar%outdir) // 'HI_amplitudes.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(35,file=title,status="old",position="append",action="write") 
       else
          open(35,file=title,status="new",action="write")
       end if
       write(35,fmt=fmt) component_list(1)%p%template_amplitudes(:,1)
       close(35)
       
       title = trim(dpar%outdir)//'HI_chisq.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(36,file=title,status="old",position="append",action="write") 
       else
          open(36,file=title,status="new",action="write")
       end if
       call compute_chisq(dat,dpar)
       write(36,'(E17.8)') dat%chisq
       close(36)

       ! title = trim(dpar%outdir)//'band_gains.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(37,file=title,status="old",position="append",action="write") 
       ! else
       !    open(37,file=title,status="new",action="write")
       ! end if
       ! write(37,fmt=fmt) dat%gain
       ! close(37)

       ! title = trim(dpar%outdir)//'band_offsets.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(38,file=title,status="old",position="append",action="write") 
       ! else
       !    open(38,file=title,status="new",action="write")
       ! end if
       ! write(38,fmt=fmt) dat%offset
       ! close(38)

       title = trim(dpar%outdir)//'HI_Td_mean.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(39,file=title,status="old",position="append",action="write") 
       else
          open(39,file=title,status="new",action="write")
       end if
       write(39,'(E17.8)') mask_avg(component_list(1)%p%indices(:,1,1),dat%masks(:,1))
       close(39)

       ! title = trim(dpar%outdir)//'band_chisq.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(40,file=title,status="old",position="append",action="write") 
       ! else
       !    open(40,file=title,status="new",action="write")
       ! end if
       ! write(40,fmt=fmt) dat%band_chisq
       ! close(40)

    end if
    
  end subroutine write_data

end module dang_data_mod
