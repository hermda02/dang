module dang_data_mod
  use healpix_types
  use dang_util_mod
  use dang_param_mod
  use dang_component_mod
  implicit none

  type, public                                :: data

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
    real(dp), allocatable, dimension(:,:,:,:) :: fg_map     ! Component maps

    real(dp), allocatable, dimension(:)       :: gain       ! Where band gains are stored
    real(dp), allocatable, dimension(:)       :: offset     ! Where band offsets are stored
    
    real(dp), allocatable, dimension(:,:)     :: masks      ! Where masks are stored

    real(dp), allocatable, dimension(:,:,:)   :: temps      ! Where template maps are stored   
    real(dp), allocatable, dimension(:,:,:)   :: temp_amps  ! Where template amplitudes are stored
    real(dp), allocatable, dimension(:,:)     :: temp_norm  ! Where we store template normalizations

  end type data

  private :: i, j, k, l
  integer(i4b) :: i, j, k, l

contains

  subroutine init_fg_map(self,npix,nmaps,nbands,nfgs)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: npix, nmaps, nbands, nfgs
    
    allocate(self%fg_map(0:npix-1,nmaps,nbands,nfgs))
    
  end subroutine init_fg_map

  subroutine init_mask_maps(self,npix,nmaps)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: npix, nmaps
    
    allocate(self%masks(0:npix-1,nmaps))
    
  end subroutine init_mask_maps

  subroutine init_data_maps(self,npix,nmaps,nbands)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: npix, nmaps, nbands

    allocate(self%sig_map(0:npix-1,nmaps,nbands))
    allocate(self%rms_map(0:npix-1,nmaps,nbands))
    allocate(self%res_map(0:npix-1,nmaps,nbands))
    allocate(self%chi_map(0:npix-1,nmaps))

    allocate(self%gain(nbands))
    allocate(self%offset(nbands))

    self%gain   = 1.d0
    self%offset = 0.d0

  end subroutine init_data_maps

  subroutine init_template(self,npix,nmaps,ntemp,nbands)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: npix, nmaps, ntemp, nbands

    allocate(self%temps(0:npix-1,nmaps,ntemp))
    allocate(self%temp_amps(nbands,nmaps,ntemp))
    allocate(self%temp_norm(nmaps,ntemp))

  end subroutine init_template

  subroutine dust_correct_band(dat,param,comp,band)
    implicit none
    type(data),   intent(inout) :: dat
    type(params)                :: param
    type(component)             :: comp
    integer(i4b), intent(in)    :: band
    real(dp), allocatable, dimension(:,:,:) :: thermal_map
    integer(i4b)                :: i, j, k
    character(len=256)          :: title

    allocate(thermal_map(0:npix-1,nmaps,nbands))

    if (trim(param%dust_corr_type) == 'uniform') then
       comp%T_d    = param%mbb_gauss(1,1)
       comp%beta_d = param%mbb_gauss(2,1)
    else if (trim(param%dust_corr_type) == 'sample') then
       if (param%mbb_gauss(1,2) .gt. 0.d0) then
          comp%T_d    = rand_normal(param%mbb_gauss(1,1),param%mbb_gauss(1,2))
       else 
          comp%T_d    = param%mbb_gauss(1,1)
       end if
       if (param%mbb_gauss(2,2) .gt. 0.d0) then
          comp%beta_d = rand_normal(param%mbb_gauss(2,1),param%mbb_gauss(2,2))
       else
          comp%beta_d = param%mbb_gauss(2,1)
       end if
    else if (trim(param%dust_corr_type) == 'planck') then
       stop
    end if
    write(*,'(a,a)') 'Dust correcting band ', trim(param%dat_label(band))
    do k = param%pol_type(1), param%pol_type(size(param%pol_type))
       do i = 0, npix-1
          thermal_map(i,k,band) = dat%temps(i,k,1)*compute_spectrum(param,comp,2,param%dat_nu(band),i,k)
          dat%sig_map(i,k,band) = dat%sig_map(i,k,band) - thermal_map(i,k,band)
       end do
    end do
    title = trim(param%outdir)//trim(param%dat_label(band))//'_thermal_map.fits'
    call write_result_map(trim(title), nside, ordering, header, thermal_map(:,:,band))
  end subroutine dust_correct_band

  subroutine extrapolate_foreground(param, dat, comp, ind, map_n)
    implicit none
    type(data),   intent(inout) :: dat
    type(params)                :: param
    type(component)             :: comp
    integer(i4b), intent(in)    :: ind, map_n
    integer(i4b)                :: i, j, k

    do i = 0, npix-1
       do j = 1, nbands
          dat%fg_map(i,map_n,j,ind) = dat%fg_map(i,map_n,param%fg_ref_loc(ind),ind)*compute_spectrum(param,comp,ind,param%dat_nu(j),i,map_n)
       end do
    end do

  end subroutine extrapolate_foreground

  subroutine extrapolate_template(param, dat, comp, ind, map_n)
    implicit none
    type(data),   intent(inout) :: dat
    type(params)                :: param
    type(component)             :: comp
    integer(i4b), intent(in)    :: ind, map_n
    integer(i4b)                :: i, j, k

    do i = 0, npix-1
       do j = 1, nbands
          dat%fg_map(i,map_n,j,param%ncomp+ind) = dat%temp_amps(j,map_n,ind)*dat%temps(i,map_n,ind)
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

  function a2f(nu)
    ! [MJy/sr / uK_RJ]
    ! Assume that nu is in GHz
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: a2f, y

    a2f = compute_bnu_prime_RJ(nu)

  end function a2f

  function a2t(nu)
    ! [uK_cmb/uK_RJ]
    ! Assume that nu is in GHz
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: a2t, y

    y = (h*nu)/(k_B*T_CMB)

    a2t = (exp(y)-1.d0)**2/(y**2*exp(y))

  end function a2t

  function f2t(nu)
    ! [uK_cmb/MJysr-2]
    ! Assume that nu is in GHz
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: f2t, y

    f2t = 1.d0/(compute_bnu_prime(nu))*1.0d-14

  end function f2t

  subroutine convert_maps(self,dat)
    ! We want to run everything in uK_RJ, yeah?
    implicit none
    type(params), intent(inout) :: self
    type(data),   intent(inout) :: dat
    integer(i4b)                :: j
    
    do j = 1, nbands
       if (.not. self%bp_map(j)) then
          if (trim(self%band_unit(j)) == 'uK_RJ') then
             cycle
          else if (trim(self%band_unit(j)) == 'uK_cmb') then
             ! uK_cmb -> uK_RJ
             write(*,*) 'Putting band ', trim(self%band_label(j)), ' from uK_cmb to uK_RJ.'
             dat%sig_map(:,:,j) = dat%sig_map(:,:,j)/a2t(self%band_nu(j)*1.0d9)
             dat%rms_map(:,:,j) = dat%rms_map(:,:,j)/a2t(self%band_nu(j)*1.0d9)
          else if (trim(self%band_unit(j)) == 'MJy/sr') then
             ! MJy/sr -> uK_RJ
             write(*,*) 'Putting band ', trim(self%band_label(j)), ' from MJy/sr to uK_RJ'
             dat%sig_map(:,:,j) = dat%sig_map(:,:,j)/a2f(self%band_nu(j)*1.0d9)
             dat%rms_map(:,:,j) = dat%rms_map(:,:,j)/a2f(self%band_nu(j)*1.0d9)
          else
             write(*,*) 'Not a unit, dumbass!'
             stop
          end if
       end if
    end do
    
  end subroutine convert_maps
  
  subroutine convert_maps_bp(self,dat)
    implicit none
    type(params), intent(inout) :: self
    type(data),   intent(inout) :: dat
    integer(i4b)                :: j
    
    do j = 1, nbands
       if (self%bp_map(j)) then
          if (trim(self%band_unit(j)) == 'uK_RJ') then
             cycle
          else if (trim(self%band_unit(j)) == 'uK_cmb') then
             ! uK_cmb -> uK_RJ
             write(*,*) 'Putting band ', trim(self%band_label(j)), ' from uK_cmb to uK_RJ.'
             dat%sig_map(:,:,j) = dat%sig_map(:,:,j)/a2t(self%band_nu(j)*1.0d9)
             dat%rms_map(:,:,j) = dat%rms_map(:,:,j)/a2t(self%band_nu(j)*1.0d9)
          else if (trim(self%band_unit(j)) == 'MJy/sr') then
             ! MJy/sr -> uK_RJ
             write(*,*) 'Putting band ', trim(self%band_label(j)), ' from MJy/sr to uK_RJ'
             dat%sig_map(:,:,j) = dat%sig_map(:,:,j)/a2f(self%band_nu(j)*1.0d9)
             dat%rms_map(:,:,j) = dat%rms_map(:,:,j)/a2f(self%band_nu(j)*1.0d9)
          else
             write(*,*) 'Not a unit, dumbass!'
             stop
          end if
       end if
    end do
    
  end subroutine convert_maps_bp

  subroutine compute_chisq(param,dat,comp,map_n)
    use healpix_types
    implicit none
    type(params)                                                :: param
    type(data),                                   intent(inout) :: dat
    type(component)                                             :: comp
    integer(i4b),                                 intent(in)    :: map_n
    real(dp)                                                    :: s, signal
    integer(i4b)                                                :: i, j, w
    
    if (trim(param%mode) == 'comp_sep') then
       dat%chisq = 0.d0
       do k = param%pol_type(1), param%pol_type(size(param%pol_type))
          do i = 0, npix-1
             if (dat%masks(i,1) == missval .or. dat%masks(i,1) == 0.d0) cycle
             do j = 1, nbands
                s = 0.d0
                do w = 1, nfgs
                   signal = dat%fg_map(i,k,j,w)
                   s = s + signal
                end do
                dat%chisq = dat%chisq + (((dat%sig_map(i,k,j) - s)**2))/(dat%rms_map(i,k,j)**2)
             end do
          end do
       end do
       dat%chisq = dat%chisq/(size(param%pol_type)*(nump*nbands)-npixpar-nglobalpar)
    else if (trim(param%mode) == 'hi_fit') then
       dat%chisq = 0.d0
       do i = 0, npix-1
          if (dat%masks(i,1) == missval .or. dat%masks(i,1) == 0.d0) cycle
          do j = 1, nbands    
             s = 0.0
             s = dat%gain(j)*comp%HI_amps(j)*comp%HI(i,1)*planck(param%dat_nu(j)*1d9,comp%T_d(i,1))+dat%offset(j)
             dat%chisq = dat%chisq + (dat%sig_map(i,map_n,j)-s)**2.d0/(dat%rms_map(i,map_n,j)**2.d0)
          end do
       end do
       dat%chisq = dat%chisq/(nump*nbands)
    end if
    
  end subroutine compute_chisq

  subroutine write_stats_to_term(param,dat,comp,iter)
    implicit none

    type(params)                        :: param
    type(data)                          :: dat
    type(component)                     :: comp
    integer(i4b),            intent(in) :: iter
    integer(i4b)                        :: k

    write(*,fmt='(a)') '---------------------------------------------'
    do k = param%pol_type(1), param%pol_type(size(param%pol_type))
       if (rank == master) then
          if (mod(iter, 1) == 0 .or. iter == 1) then
             write(*,fmt='(i6, a, a, f7.3, a, f8.4, a, 10e10.3)')&
                  iter, " - Poltype: "//trim(tqu(k)), " - A_s: ",&
                  dat%fg_map(23000,k,param%fg_ref_loc(1),1),  " - beta_s: ",&
                  mask_avg(comp%beta_s(:,k),dat%masks(:,1)), ' - A_d: ', &
                  dat%temp_amps(:,k,1)/dat%temp_norm(k,1)
             write(*,fmt='(a)') '---------------------------------------------'
          end if
       end if
    end do
    call compute_chisq(param,dat,comp,k)
    if (rank == master) then
       if (mod(iter, 1) == 0 .or. iter == 1) then
          write(*,fmt='(i6,a,e10.3)') iter, " - Chisq: ", dat%chisq
          write(*,fmt='(a)') '---------------------------------------------'
       end if
    end if

  end subroutine write_stats_to_term

  ! Data output routines
  subroutine write_maps(param,dat,comp)
    implicit none
    
    type(params)                        :: param
    type(data)                          :: dat
    type(component)                     :: comp
    real(dp), dimension(0:npix-1,nmaps) :: map
    real(dp)                            :: s, signal
    integer(i4b)                        :: n, mn, i, j, k, l
    character(len=128)                  :: title

    write(*,*) 'Output data maps'
    
    if (trim(param%mode) == 'comp_sep') then
       
       write(iter_str, '(i0.5)') iter
       if (param%output_fg .eqv. .true.) then
          do j = 1, nbands
             do n = 1, param%ntemp
                title = trim(param%outdir) // trim(param%dat_label(j)) //'_'// trim(param%temp_label(n)) //&
                     '_k' // trim(iter_str) // '.fits'
                map(:,:)   = dat%fg_map(:,:,j,n+param%ncomp)
                do i = 0, npix-1
                   if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                      map(i,:) = missval
                   end if
                end do
                call write_result_map(trim(title), nside, ordering, header, map)
             end do
             do n = 1, param%ncomp
                title = trim(param%outdir) // trim(param%dat_label(j)) //'_'// trim(param%fg_label(n)) //&
                     '_amplitude_k' // trim(iter_str) // '.fits'
                map(:,:)   = dat%fg_map(:,:,j,n)
                do i = 0, npix-1
                   if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                      map(i,:) = missval
                   end if
                end do
                call write_result_map(trim(title), nside, ordering, header, map)
             end do
          end do
       else 
          do n = 1, param%ncomp
             title = trim(param%outdir) // trim(param%dat_label(param%fg_ref_loc(n))) //'_'// trim(param%fg_label(n)) //&
                  '_amplitude_k' // trim(iter_str) // '.fits'
             map(:,:)   = dat%fg_map(:,:,param%fg_ref_loc(n),n)
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                   map(i,:) = missval
                end if
             end do
             call write_result_map(trim(title), nside, ordering, header, map)
          end do
       end if
       do j = 1, nbands
          title = trim(param%outdir) // trim(param%dat_label(j)) // '_residual_k' // trim(iter_str) // '.fits'
          map(:,:)   = dat%res_map(:,:,j)
          do i = 0, npix-1
             if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                map(i,:) = missval
             end if
          end do
          call write_result_map(trim(title), nside, ordering, header, map)
       end do
       title = trim(param%outdir) // 'synch_beta_k' // trim(iter_str) // '.fits'
       map(:,:)   = comp%beta_s(:,:)
       do i = 0, npix-1
          if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
             map(i,:) = missval
          end if
       end do
       call write_result_map(trim(title), nside, ordering, header, map)
       do mn = 1, nmaps
          dat%chi_map(:,mn) = 0.d0
          do i = 0, npix-1
             do j = 1, nbands
                s      = 0.d0
                do l = 1, nfgs
                   signal = dat%fg_map(i,mn,j,l)
                   s      = s + signal
                end do
                dat%chi_map(i,mn) = dat%chi_map(i,mn) + dat%masks(i,1)*(dat%sig_map(i,mn,j) - s)**2.d0/dat%rms_map(i,mn,j)**2.d0
             end do
          end do
       end do
       dat%chi_map(:,:) = dat%chi_map(:,:)/(nbands)
       title = trim(param%outdir) // 'chisq_k'// trim(iter_str) // '.fits'
       map(:,:)   = dat%chi_map(:,:)
       do i = 0, npix-1
          if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
             map(i,:) = missval
             dat%chi_map(i,:) = missval
          end if
       end do
       call write_result_map(trim(title), nside, ordering, header, map)

    ! For hi_fit maps
    else if (trim(param%mode) == 'hi_fit') then

       write(iter_str, '(i0.5)') iter
       n = 1
       do j = 1, param%numband
          if (param%band_inc(j)) then
             title = trim(param%outdir) // trim(param%dat_label(j)) //'_hi_amplitude_k'// trim(iter_str) // '.fits'
             map(:,1)   = comp%HI_amps(n)*comp%HI(:,1)
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                   map(i,1) = missval
                end if
             end do
             call write_bintab(map,npix,1, header, nlheader, trim(title))
             n = n + 1
          end if
       end do

       n = 1
       do j = 1, param%numband
          if (param%band_inc(j)) then
             title = trim(param%outdir) // trim(param%dat_label(j)) // '_residual_k' // trim(iter_str) // '.fits'
             map(:,:)   = dat%res_map(:,:,n)
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                   map(i,1) = missval
                end if
             end do
             call write_bintab(map,npix,1, header, nlheader, trim(title))
             n = n + 1
          end if
       end do
       title = trim(param%outdir) // 'T_d_k'// trim(iter_str) // '.fits'
       map(:,1)   = comp%T_d(:,1)
       do i = 0, npix-1
          if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
             map(i,1) = missval
          end if
       end do
       call write_bintab(map,npix,1, header, nlheader, trim(title))
       dat%chi_map = 0.d0
       do i = 0, npix-1
          do j = 1, nbands
             s =  comp%HI_amps(j)*comp%HI(i,1)*planck(param%dat_nu(j)*1d9,comp%T_d(i,1))
             dat%chi_map(i,1) = dat%chi_map(i,1) + dat%masks(i,1)*(dat%sig_map(i,1,j) - s)**2.d0/dat%rms_map(i,1,j)**2.d0
          end do
       end do
       dat%chi_map(:,1) = dat%chi_map(:,1)/(nbands+nfgs)
       title = trim(param%outdir) // 'chisq_k' // trim(iter_str) // '.fits'
       map(:,1)   = dat%chi_map(:,1)
       do i = 0, npix-1
          if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
             map(i,1) = missval
          end if
       end do
       call write_bintab(map,npix,1, header, nlheader, trim(title))
          
    end if
    
  end subroutine write_maps
  
  subroutine write_data(param,dat,comp,map_n)
    implicit none
    type(params)                        :: param
    type(data)                          :: dat
    type(component)                     :: comp
    integer(i4b),            intent(in) :: map_n
    integer(i4b)                        :: i
    character(len=2)                    :: temp_n
    character(len=128)                  :: title

    if (trim(param%mode) == 'comp_sep') then
       
       title = trim(param%outdir) // 'pixel_23000_A_d_' // trim(tqu(map_n)) // '.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(30,file=title, status="old",position="append", action="write")
       else
          open(30,file=title, status="new", action="write")
       endif
       write(30,*) dat%fg_map(23000,map_n,param%fg_ref_loc(1),2)
       close(30)
       
       title = trim(param%outdir) // 'pixel_23000_A_s_' // trim(tqu(map_n)) // '.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(31,file=title, status="old",position="append", action="write")
       else
          open(31,file=title, status="new", action="write")
       endif
       write(31,*) dat%fg_map(23000,map_n,param%fg_ref_loc(1),1)
       close(31)
       
       title = trim(param%outdir) // 'pixel_23000_beta_s_' // trim(tqu(map_n)) // '.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(32,file=title, status="old",position="append", action="write")
       else
          open(32,file=title, status="new", action="write")
       endif
       write(32,*) comp%beta_s(23000,map_n)
       close(32)
       
       title = trim(param%outdir) // 'total_chisq_' // trim(tqu(map_n)) // '.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(33,file=title, status="old",position="append", action="write")
       else
          open(33,file=title, status="new", action="write")
       endif
       call compute_chisq(param,dat,comp,map_n)
       write(33,*) dat%chisq
       close(33)
       
       do i = 1, param%ntemp
          title = trim(param%outdir) //  trim(param%temp_label(i)) // '_' //trim(tqu(map_n)) // '_amplitudes.dat'
          inquire(file=title,exist=exist)
          if (exist) then
             open(34,file=title, status="old", &
                  position="append", action="write")
          else
             open(34,file=title, status="new", action="write")
          endif
          write(34,'(10(E17.8))') dat%temp_amps(:,map_n,i)
          close(34)
       end do
       
    else if (trim(param%mode) == 'hi_fit') then
       title = trim(param%outdir) // 'HI_amplitudes.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(35,file=title,status="old",position="append",action="write") 
       else
          open(35,file=title,status="new",action="write")
       end if
       write(35,'(10(E17.8))') comp%HI_amps
       close(35)
       
       title = trim(param%outdir)//'HI_chisq.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(36,file=title,status="old",position="append",action="write") 
       else
          open(36,file=title,status="new",action="write")
       end if
       call compute_chisq(param,dat,comp,1)
       write(36,'(E17.8)') dat%chisq
       close(36)

       ! title = trim(param%outdir)//'band_gains.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(37,file=title,status="old",position="append",action="write") 
       ! else
       !    open(37,file=title,status="new",action="write")
       !    write(37,"(3x,8(A16))") param%dat_label
       ! end if
       ! write(37,"(8(E16.4))") dat%gain
       ! close(37)

       ! title = trim(param%outdir)//'band_offsets.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(38,file=title,status="old",position="append",action="write") 
       ! else
       !    open(38,file=title,status="new",action="write")
       !    write(38,"(3x,8(A16))") param%dat_label
       ! end if
       ! write(38,"(8(E16.4))") dat%offset
       ! close(38)
       
    end if
    
  end subroutine write_data

end module dang_data_mod
