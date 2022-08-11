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
    real(dp), allocatable, dimension(:,:,:,:) :: fg_map     ! Component maps (npix, nmaps, nbands, nfgs)
    real(dp), allocatable, dimension(:,:,:)   :: sky_model  ! Total sky model

    real(dp), allocatable, dimension(:)       :: gain       ! Where band gains are stored
    real(dp), allocatable, dimension(:)       :: offset     ! Where band offsets are stored
    
    real(dp), allocatable, dimension(:,:)     :: masks      ! Where masks are stored

    real(dp), allocatable, dimension(:,:,:)   :: temps      ! Where template maps are stored   
    real(dp), allocatable, dimension(:,:,:)   :: temp_amps  ! Where template amplitudes are stored
    real(dp), allocatable, dimension(:,:)     :: temp_norm  ! Where we store template normalizations

    real(dp), allocatable, dimension(:)       :: amp_vec    ! Amplitude vector returned from the CG solver
    real(dp), allocatable, dimension(:)       :: band_chisq ! A variable to store the chisq of each band

  contains
    procedure :: init_data_maps
    procedure :: read_data_maps

  end type dang_data

  private :: i, j, k, l
  integer(i4b) :: i, j, k, l

contains

  ! I think I want a constructor that creates a data object for each band


  subroutine init_data_maps(self,dpar)
    implicit none
    class(dang_data),       intent(inout) :: self
    type(dang_params),      intent(in)    :: dpar

    write(*,*) "Init data maps"
    
    allocate(self%fg_map(0:npix-1,nmaps,0:nbands,nfgs))    
    self%fg_map(:,:,:,:) = 0.d0

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
    allocate(self%band_chisq(nbands))
    
    self%gain   = 1.d0
    self%offset = 0.d0

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

  subroutine update_sky_model(ddata)
    implicit none
    type(dang_data),   intent(inout) :: ddata
    type(dang_comps),  pointer       :: c
    integer(i4b)                     :: i, j, k, l

    ddata%sky_model(:,:,:) = 0.d0
    do l = 1, ncomp
       c => component_list(l)%p
       if (c%type /= 'template') then
          do i = 0, npix-1
             do k = 1, nmaps
                do j = 1, nbands
                   ddata%sky_model(i,k,j) = ddata%sky_model(i,k,j) + &
                        & c%amplitude(i,k)*c%eval_sed(j,i,k)
                end do
             end do
          end do
       else
          do i = 0, npix-1
             do k = 1, nmaps
                do j = 1, nbands
                   ddata%sky_model(i,k,j) = ddata%sky_model(i,k,j) + &
                        & c%template(i,k)*c%template_amplitudes(j,k)
                end do
             end do
          end do
       end if
    end do

    ddata%res_map = ddata%sig_map - ddata%sky_model

  end subroutine update_sky_model

  subroutine mask_hi(self, dpar, dcomps)
    implicit none
    type(dang_data),             intent(inout) :: self
    type(dang_params)                          :: dpar
    type(dang_comps)                           :: dcomps
    integer(i4b)                               :: i, j

    do i = 0, npix-1
       if (dcomps%HI(i,1) > dpar%thresh) then
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
  
  subroutine dust_correct_band(self,dpar,comp,band,iter)
    implicit none
    type(dang_data),             intent(inout) :: self
    type(dang_params)                          :: dpar
    type(dang_comps)                           :: comp
    integer(i4b),                intent(in)    :: band
    integer(i4b), optional,      intent(in)    :: iter
    real(dp), allocatable, dimension(:,:,:)    :: thermal_map
    integer(i4b)                               :: i, j, k
    character(len=256)                         :: title

    allocate(thermal_map(0:npix-1,nmaps,nbands))

    if (trim(dpar%dust_corr_type) == 'uniform') then
       comp%T_d    = dpar%mbb_gauss(1,1)
       comp%beta_d = dpar%mbb_gauss(2,1)
    else if (trim(dpar%dust_corr_type) == 'sample') then
       if (dpar%mbb_gauss(1,2) .gt. 0.d0) then
          comp%T_d    = rand_normal(dpar%mbb_gauss(1,1),dpar%mbb_gauss(1,2))
       else 
          comp%T_d    = dpar%mbb_gauss(1,1)
       end if
       if (dpar%mbb_gauss(2,2) .gt. 0.d0) then
          comp%beta_d = rand_normal(dpar%mbb_gauss(2,1),dpar%mbb_gauss(2,2))
       else
          comp%beta_d = dpar%mbb_gauss(2,1)
       end if
    else if (trim(dpar%dust_corr_type) == 'planck') then
       stop
    end if
    write(*,'(a,a)') 'Dust correcting band ', trim(dpar%band_label(band))
    do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
       do i = 0, npix-1
          thermal_map(i,k,band) = self%temps(i,k,1)*compute_spectrum(dpar,comp,bp(band),2,i,k)
          self%sig_map(i,k,band) = self%sig_map(i,k,band) - thermal_map(i,k,band)
       end do
    end do

    if (present(iter)) then
       write(iter_str, '(i0.5)') iter
       title = trim(dpar%outdir)//trim(dpar%band_label(band))//'_thermal_map_k' // trim(iter_str) // '.fits'
    else
       title = trim(dpar%outdir)//trim(dpar%band_label(band))//'_thermal_map.fits'
    end if
    call write_result_map(trim(title), nside, ordering, header, thermal_map(:,:,band))
  end subroutine dust_correct_band

  subroutine extrapolate_foreground(dpar, dat, comp, ind, map_n)
    implicit none
    type(dang_data), intent(inout) :: dat
    type(dang_params)              :: dpar
    type(dang_comps)               :: comp
    integer(i4b),    intent(in)    :: ind, map_n
    integer(i4b)                   :: i, j, k

    do i = 0, npix-1
       do j = 1, nbands
          dat%fg_map(i,map_n,j,ind) = dat%fg_map(i,map_n,0,ind)*compute_spectrum(dpar,comp,bp(j),ind,i,map_n)!ind,dpar%band_nu(j),i,map_n)
       end do
    end do

  end subroutine extrapolate_foreground

  subroutine extrapolate_template(dpar, dat, comp, ind, map_n)
    implicit none
    type(dang_data), intent(inout) :: dat
    type(dang_params)              :: dpar
    type(dang_comps)               :: comp
    integer(i4b),    intent(in)    :: ind, map_n
    integer(i4b)                   :: i, j, k

    do i = 0, npix-1
       do j = 1, nbands
          dat%fg_map(i,map_n,j,dpar%ncomp+ind) = dat%temp_amps(j,map_n,ind)*dat%temps(i,map_n,ind)
       end do
    end do

  end subroutine extrapolate_template

  function compute_bnu_prime_RJ(nu)
    ! Assume that nu is in GHz
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: compute_bnu_prime_RJ

    compute_bnu_prime_RJ = 2.d0*k_B*nu**2.d0/c**2.d0

  end function compute_bnu_prime_RJ

  function compute_bnu_prime(nu)
    ! Assume that nu is in GHz
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: compute_bnu_prime, y

    y = h*nu/(k_B*T_CMB)
    compute_bnu_prime = (2.d0*h*nu**3)/(c**2.d0*(exp(y)-1))*(exp(y)/(exp(y)-1))*h*nu/(k_B*T_CMB**2)

  end function compute_bnu_prime

  function a2f(bp)
    ! [MJy/sr / uK_RJ]
    ! Assume that nu is in GHz
    implicit none
    type(bandinfo), intent(in) :: bp
    real(dp)                   :: a2f, y, sum
    integer(i4b)               :: i

    sum = 0.d0

    if (bp%id == 'delta') then
       sum = compute_bnu_prime_RJ(bp%nu_c*1d9)
    else
       do i = 1, bp%n
          sum = sum + bp%tau0(i)*compute_bnu_prime_RJ(bp%nu0(i)*1d9)
       end do
    end if
    a2f = sum

  end function a2f

  function a2t(bp)
    ! [uK_cmb/uK_RJ]
    ! Assume that nu is in GHz
    implicit none
    type(bandinfo), intent(in) :: bp
    real(dp)                   :: a2t, y, sum 
    integer(i4b)               :: i

    sum = 0.d0

    if (bp%id == 'delta') then
       if (bp%nu_c > 1e7) then
          y = (h*bp%nu_c)/(k_B*T_CMB)
       else
          y = (h*bp%nu_c*1d9)/(k_B*T_CMB)
       end if
       sum = (exp(y)-1.d0)**2/(y**2*exp(y))
    else
       do i = 1, bp%n
          if (bp%nu0(i) > 1e7) then
             y = (h*bp%nu0(i))/(k_B*T_CMB)
             sum = sum + bp%tau0(i)*(exp(y)-1.d0)**2/(y**2*exp(y))
          else
             y = (h*bp%nu0(i)*1d9)/(k_B*T_CMB)
             sum = sum + bp%tau0(i)*(exp(y)-1.d0)**2/(y**2*exp(y))
          end if
       end do
    end if

    a2t = sum

  end function a2t

  function f2t(bp)
    ! [uK_cmb/MJysr-1]
    ! Assume that nu is in GHz
    implicit none
    type(bandinfo), intent(in) :: bp
    real(dp)                   :: f2t, sum
    integer(i4b)               :: i

    sum = 0.d0

    if (bp%id == 'delta') then
       if (bp%nu_c > 1e7) then
          sum = 1.d0/(compute_bnu_prime(bp%nu_c))*1.0d-14
       else
          sum = 1.d0/(compute_bnu_prime(bp%nu_c*1d9))*1.0d-14
       end if
    else
       do i = 1, bp%n 
          if (bp%nu0(i) > 1e7) then
             sum = sum + bp%tau0(i)/(compute_bnu_prime(bp%nu0(i)))*1.0d-14
          else
             sum = sum + bp%tau0(i)/(compute_bnu_prime(bp%nu0(i)*1d9))*1.0d-14
          end if
       end do
    end if

    f2t = sum

  end function f2t

  subroutine convert_maps(self,dpar)
    ! We want to run everything in uK_RJ (at least for compsep), yeah?
    implicit none
    type(dang_data),   intent(inout) :: self
    type(dang_params), intent(inout) :: dpar
    integer(i4b)                :: j
    
    if (dpar%mode == 'comp_sep') then

       do j = 1, nbands
          if (.not. dpar%bp_map(j)) then
             if (trim(dpar%band_unit(j)) == 'uK_RJ') then
                cycle
             else if (trim(dpar%band_unit(j)) == 'uK_cmb') then
                ! uK_cmb -> uK_RJ
                ! Check bandpass type
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
  
  subroutine convert_bp_maps(self,dpar)!,bp)
    implicit none
    type(dang_data),   intent(inout) :: self
    type(dang_params), intent(inout) :: dpar
    ! type(bandinfo),    intent(in)    :: bp
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

  subroutine compute_chisq(self,dpar)!,comp,map_n)
    use healpix_types
    implicit none
    type(dang_data),                              intent(inout) :: self
    type(dang_params)                                           :: dpar
    ! type(dang_comps)                                            :: comp
    ! integer(i4b),                                 intent(in)    :: map_n
    real(dp)                                                    :: s, signal
    integer(i4b)                                                :: i, j, k

    self%chisq = 0.d0
    do i = 0, npix-1
       if (self%masks(i,1) == missval .or. self%masks(i,1) == 0.d0) cycle
       do k = 2,3!dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
          do j = 1, nbands
             self%chisq = self%chisq + (self%sig_map(i,k,j) - self%sky_model(i,k,j))**2.d0 / &
                  & (self%rms_map(i,k,j)**2.d0)
          end do
       end do
    end do
    self%chisq = self%chisq/(size(dpar%pol_type)*(nump*nbands)-npixpar-nglobalpar)
    ! if (trim(dpar%mode) == 'comp_sep') then
    !    self%chisq = 0.d0
    !    do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
    !       do i = 0, npix-1
    !          if (self%masks(i,1) == missval .or. self%masks(i,1) == 0.d0) cycle
    !          do j = 1, nbands
    !             s = 0.d0
    !             do w = 1, nfgs
    !                signal = self%fg_map(i,k,j,w)
    !                s = s + signal
    !             end do
    !             self%chisq = self%chisq + (((self%sig_map(i,k,j) - s)**2))/(self%rms_map(i,k,j)**2)
    !          end do
    !       end do
    !    end do
    !    self%chisq = self%chisq/(size(dpar%pol_type)*(nump*nbands)-npixpar-nglobalpar)
    ! else if (trim(dpar%mode) == 'hi_fit') then
    !    self%chisq = 0.d0
    !    do i = 0, npix-1
    !       if (self%masks(i,1) == missval .or. self%masks(i,1) == 0.d0) cycle
    !       do j = 1, nbands    
    !          s = 0.0
    !          s = self%gain(j)*comp%HI_amps(j)*comp%HI(i,1)*planck(bp(j),comp%T_d(i,1))+self%offset(j)
    !          self%chisq = self%chisq + (self%sig_map(i,map_n,j)-s)**2.d0/(self%rms_map(i,map_n,j)**2.d0)
    !       end do
    !    end do
    !    self%chisq = self%chisq!/((nump*nbands)-npixpar-nglobalpar)
    ! end if
    
  end subroutine compute_chisq

  subroutine write_stats_to_term(self,dpar,dcomps,iter)
    implicit none

    type(dang_data)                     :: self
    type(dang_params)                   :: dpar
    type(dang_comps)                    :: dcomps
    integer(i4b),            intent(in) :: iter
    integer(i4b)                        :: k

    if (trim(dpar%mode) == 'comp_sep') then
       write(*,fmt='(a)') '---------------------------------------------'
       do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
          if (rank == master) then
             if (mod(iter, 1) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, a, f7.3, a, f8.4, a, 10e10.3)')&
                     iter, " - Poltype: "//trim(tqu(k)), " - A_s: ",&
                     component_list(1)%p%amplitude(23000,k),  " - beta_s: ",&
                     mask_avg(component_list(1)%p%indices(:,k,1),self%masks(:,1)), ' - A_d: ', &
                     component_list(2)%p%template_amplitudes(:,k)
                write(*,fmt='(a)') '---------------------------------------------'
             end if
          end if
       end do
       call compute_chisq(self,dpar)
       if (rank == master) then
          if (mod(iter, 1) == 0 .or. iter == 1) then
             write(*,fmt='(i6,a,f16.5)') iter, " - Chisq: ", self%chisq
             write(*,fmt='(a)') '---------------------------------------------'
          end if
       end if
       
    else if (trim(dpar%mode) == 'hi_fit') then
       call compute_chisq(self,dpar)!,dcomps,1)
       if (rank == master) then
          if (mod(iter, 1) == 0 .or. iter == 1) then
             if (nbands .lt. 10) then
                write(*,fmt='(i6, a, f10.5, a, f10.3, a, 10e10.3)')&
                     iter, " - chisq: " , self%chisq, " - T_d: ",&
                     mask_avg(dcomps%T_d(:,1),self%masks(:,1)), ' - A_HI: ', dcomps%HI_amps
                write(*,fmt='(a)') '---------------------------------------------'
             else
                write(*,fmt='(i6, a, E10.3, a, e10.3)')&
                     iter, " - chisq: " , self%chisq, " - T_d: ",&
                     mask_avg(dcomps%T_d(:,1),self%masks(:,1))
                write(*,fmt='(a)') '---------------------------------------------'
             end if
          end if
       end if
    end if

  end subroutine write_stats_to_term

  ! Data output routines
  subroutine write_maps(dpar,dat)!,dcomps)
    implicit none
    
    type(dang_params)                   :: dpar
    type(dang_data)                     :: dat
    type(dang_comps),         pointer   :: c
    real(dp), dimension(0:npix-1,nmaps) :: map
    real(dp)                            :: s, signal
    integer(i4b)                        :: n, mn, i, j, k, l
    character(len=128)                  :: title

    write(*,*) 'Output data maps'
    
    if (trim(dpar%mode) == 'comp_sep') then
       
       write(iter_str, '(i0.5)') iter
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
       else 
          do n = 1, dpar%ncomp
             title = trim(dpar%outdir) // 'reference_'// trim(dpar%fg_label(n)) //&
                  '_amplitude_k' // trim(iter_str) // '.fits'
             map(:,:)   = component_list(n)%p%amplitude!dat%fg_map(:,:,0,n)
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                   map(i,:) = missval
                end if
             end do
             call write_result_map(trim(title), nside, ordering, header, map)
          end do
       end if
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
       do n = 1, ncomp
          c => component_list(n)%p
          if (trim(c%type) == 'template') then
             cycle
          else if (trim(c%type) == 'power-law') then
             title = trim(dpar%outdir) // trim(c%label) // '_beta_k' // trim(iter_str) // '.fits'
             map(:,:)   = c%indices(:,:,1)
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                   map(i,:) = missval
                end if
             end do
          else if (trim(c%type) == 'mbb') then
             title = trim(dpar%outdir) // trim(c%label) // '_beta_k' // trim(iter_str) // '.fits'
             map(:,:)   = c%indices(:,:,1)
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                   map(i,:) = missval
                end if
             end do
             title = trim(dpar%outdir) // trim(c%label) // '_T_k' // trim(iter_str) // '.fits'
             map(:,:)   = c%indices(:,:,2)
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                   map(i,:) = missval
                end if
             end do
          end if
       end do
       call write_result_map(trim(title), nside, ordering, header, map)
       do mn = 1, nmaps
          dat%chi_map(:,mn) = 0.d0
          do i = 0, npix-1
             do j = 1, nbands
                s      = 0.d0
                dat%chi_map(i,mn) = dat%chi_map(i,mn) + dat%masks(i,1)*(dat%sig_map(i,mn,j) - dat%sky_model(i,mn,j))**2.d0&
                     & /dat%rms_map(i,mn,j)**2.d0
             end do
          end do
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
       

       ! write(iter_str, '(i0.5)') iter
       ! do j = 1, nbands
       !    title = trim(dpar%outdir) // trim(dpar%band_label(j)) //'_hi_amplitude_k'// trim(iter_str) // '.fits'
       !    map(:,1)   = dcomps%HI_amps(j)*dcomps%HI(:,1)
       !    do i = 0, npix-1
       !       if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
       !          map(i,1) = missval
       !       end if
       !    end do
       !    call write_bintab(map,npix,1, header, nlheader, trim(title))
       ! end do

       ! ! Write out residual maps
       ! do j = 1, nbands
       !    title = trim(dpar%outdir) // trim(dpar%band_label(j)) // '_residual_k' // trim(iter_str) // '.fits'
       !    map(:,:)   = dat%res_map(:,:,j)
       !    do i = 0, npix-1
       !       if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
       !          map(i,1) = missval
       !       end if
       !    end do
       !    call write_bintab(map,npix,1, header, nlheader, trim(title))
       ! end do

       ! ! Write out T_d map
       ! title = trim(dpar%outdir) // 'T_d_k'// trim(iter_str) // '.fits'
       ! map(:,1)   = dcomps%T_d(:,1)
       ! do i = 0, npix-1
       !    if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
       !       map(i,1) = missval
       !    end if
       ! end do
       ! call write_bintab(map,npix,1, header, nlheader, trim(title))

       ! ! Compute and write out \chi^2 map
       ! dat%chi_map = 0.d0
       ! do i = 0, npix-1
       !    do j = 1, nbands
       !       s = dat%gain(j)*dcomps%HI_amps(j)*dcomps%HI(i,1)*planck(bp(j),dcomps%T_d(i,1))+dat%offset(j)
       !       dat%chi_map(i,1) = dat%chi_map(i,1) + dat%masks(i,1)*(dat%sig_map(i,1,j) - s)**2.d0/dat%rms_map(i,1,j)**2.d0
       !    end do
       ! end do
       ! dat%chi_map(:,1) = dat%chi_map(:,1)!/(nbands)
       ! title = trim(dpar%outdir) // 'chisq_k' // trim(iter_str) // '.fits'
       ! map(:,1)   = dat%chi_map(:,1)
       ! do i = 0, npix-1
       !    if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
       !       map(i,1) = missval
       !    end if
       ! end do
       ! call write_bintab(map,npix,1, header, nlheader, trim(title))
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

    if (trim(dpar%mode) == 'comp_sep') then
       
       ! title = trim(dpar%outdir) // 'pixel_23000_A_d_' // trim(tqu(map_n)) // '.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(30,file=title, status="old",position="append", action="write")
       ! else
       !    open(30,file=title, status="new", action="write")
       ! endif
       ! write(30,*) dat%fg_map(23000,map_n,2,2)
       ! close(30)
       
       ! title = trim(dpar%outdir) // 'pixel_23000_A_s_' // trim(tqu(map_n)) // '.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(31,file=title, status="old",position="append", action="write")
       ! else
       !    open(31,file=title, status="new", action="write")
       ! endif
       ! write(31,*) dat%fg_map(23000,map_n,2,1)
       ! close(31)
       
       ! title = trim(dpar%outdir) // 'pixel_23000_beta_s_' // trim(tqu(map_n)) // '.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(32,file=title, status="old",position="append", action="write")
       ! else
       !    open(32,file=title, status="new", action="write")
       ! endif
       ! write(32,*) dcomps%beta_s(23000,map_n)
       ! close(32)
       
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
          if (c%type == 'template') then
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
          end if
       end do
       
    else if (trim(dpar%mode) == 'hi_fit') then

       write(nband_str, '(i4)') nbands

       fmt = '('//trim(nband_str)//'(E17.8))'

       ! title = trim(dpar%outdir) // 'HI_amplitudes.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(35,file=title,status="old",position="append",action="write") 
       ! else
       !    open(35,file=title,status="new",action="write")
       ! end if
       ! write(35,fmt=fmt) dcomps%HI_amps
       ! close(35)
       
       ! title = trim(dpar%outdir)//'HI_chisq.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(36,file=title,status="old",position="append",action="write") 
       ! else
       !    open(36,file=title,status="new",action="write")
       ! end if
       ! call compute_chisq(dat,dpar)!,dcomps,1)
       ! write(36,'(E17.8)') dat%chisq
       ! close(36)

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

       ! title = trim(dpar%outdir)//'HI_Td_mean.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(39,file=title,status="old",position="append",action="write") 
       ! else
       !    open(39,file=title,status="new",action="write")
       ! end if
       ! write(39,'(E17.8)') mask_avg(dcomps%T_d(:,1),dat%masks(:,1))
       ! close(39)

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

  subroutine update_ddata(self,dpar,dcomps)
    implicit none
    type(dang_data), intent(inout) :: self
    type(dang_params)              :: dpar
    type(dang_comps)               :: dcomps
    
    integer(i4b)                   :: i,j
    real(dp)                       :: s

    call compute_chisq(self,dpar)!,dcomps,1)
    self%chi_map = 0.d0
    do j = 1, nbands
       self%band_chisq(j) = 0.d0
       ! Compute the residual for each map
       do i = 0, npix-1
          if (self%masks(i,1) == missval .or. self%masks(i,1) == 0.d0) then
             self%res_map(i,1,j) = missval
             self%chi_map(i,1)   = missval
             cycle
          end if
          ! Declare the model value
          s = self%gain(j)*dcomps%HI_amps(j)*dcomps%HI(i,1)*planck(bp(j),dcomps%T_d(i,1))+self%offset(j)
          ! Compute residual per band
          self%res_map(i,1,j) = (self%sig_map(i,1,j)-s)
          ! Compute the summed total chisq
          self%chi_map(i,1) = self%chi_map(i,1) + self%masks(i,1)*(self%sig_map(i,1,j) - s)**2.d0/self%rms_map(i,1,j)**2.d0
          ! Compute the chisq for each band
          self%band_chisq(j) = self%band_chisq(j) + (self%res_map(i,1,j)/self%rms_map(i,1,j))**2
       end do
    end do

  end subroutine update_ddata

end module dang_data_mod
