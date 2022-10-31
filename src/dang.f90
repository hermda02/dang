program dang
  use healpix_types
  use pix_tools
  use fitstools
  use udgrade_nr
  use dang_util_mod
  use dang_param_mod
  use dang_bp_mod
  use dang_linalg_mod
  use dang_data_mod
  use dang_component_mod
  use dang_sample_mod
  use dang_cg_mod
  use dang_swap_mod
  implicit none
  
  !------------------------------------------------------------------------------------------------------
  ! This program was designed to fit the Planck (NPIPE) 353 GHz dust map to LFI bands to set constraints|
  ! on the level of polarized emission from Anomalous Microwave Emission in the LFI bands. It has grown |
  ! into a more full Gibbs sampler, thus inpsiring the new name:                                        |
  !                                                                                                     |  
  !\                          Daniel's Amazing New Gibbs sampler (DANG)                                /|
  ! \                                                                                                 / |
  !  \                                    Daniel Herman 2020                                         /  |
  !   \                                                                                             /   |  
  !-----------------------------------------------------------------------------------------------------|  
  
  !-----------------------------------------------------------------------------------------------------|  
  ! What we want here: tale a joint fit of CMB, synchrotron, and dust emission. In order to do this     |
  ! effectively, this will be a mulit-frequency Gibbs sampler. Using the BeyondPlanck LFI maps, along   |
  ! side the WMAP data, we iteratively fit CMB (per pixel), synchrotron (per pixel), and dust (global). |
  ! Once this has been done, estimates to the level of AME polarization will be made using the lowest   |
  ! frequency bands used here.                                                                          |
  !-----------------------------------------------------------------------------------------------------|  
  integer(i4b)                          :: i, j, k, l, m, n
  integer(i4b)                          :: bp_iter, namps
  real(dp), allocatable, dimension(:,:) :: map, rms, true_synch

  ! Object Orient
  type(dang_params) :: dpar
  type(dang_data)   :: ddata
  type(dang_comps)  :: dcomps

  call init_mpi()
  call read_param_file(dpar)
  call init_bp_mod(dpar)

  !----------------------------------------------------------------------------------------------------------
  ! General paramters
  if (trim(dpar%mode) == 'comp_sep') then
     i = getsize_fits(dpar%mask_file, nside=nside, ordering=ordering, nmaps=nmaps)
  else if (trim(dpar%mode) == 'hi_fit') then
     i = getsize_fits(dpar%mask_file,    nside=nside, ordering=ordering, nmaps=nmaps)
  else
     write(*,*) 'Unrecognized operational mode. Do better!'
  end if
  ddata%npix = nside2npix(nside) 
  npix       = ddata%npix
  nbands     = dpar%numinc
  nfgs       = dpar%ncomp+dpar%ntemp
  ncomp      = dpar%ncomp
  ncg_groups = dpar%ncggroup
  nsample    = dpar%nsample
  ml_mode    = trim(dpar%ml_mode)
  npixpar    = 0.d0
  nglobalpar = 0.d0
  nump       = 0
  namps      = 0
  nlheader   = size(header)
  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------

  call RANDOM_SEED()

  ! Either do component separation, or carry out the HI fit
  !------------------------------------------------------------------------------------------------
  if (trim(dpar%mode) == 'comp_sep') then
     ! Initialize ddata and components
     !----------------------------------------------------------------------------------------------------------
     call ddata%init_data_maps(dpar)
     call ddata%read_data_maps(dpar)
     call ddata%convert_maps(dpar)
     ! do j = 1, nbands
     !    ! Check to see if any maps need to be dust corrected
     !    if (dpar%dust_corr(j) .and. .not. dpar%bp_map(j)) then
     !       call dust_correct_band(ddata,dpar,dcomps,j)
     !    end if
     ! end do
     write(*,*) ''
     call initialize_components(dpar)
     call initialize_cg_groups(dpar)
     call ddata%update_sky_model
     call write_maps(dpar,ddata)
     write(*,*) '---------------------------'
     write(*,*) ' Starting main Gibbs Chain '
     write(*,*) '---------------------------'
     call comp_sep
  
  !------------------------------------------------------------------------------------------------
  else if (trim(dpar%mode) == 'hi_fit') then
     
     ! Initialize ddata and components
     call ddata%init_data_maps(dpar)
     call ddata%read_data_maps(dpar)
     call ddata%convert_maps(dpar)
     call initialize_components(dpar)
     call initialize_cg_groups(dpar)
     call write_maps(dpar,ddata)
     write(*,*) '---------------------------'
     write(*,*) ' Starting main Gibbs Chain '
     write(*,*) '---------------------------'
     call hi_fit   
  end if

contains
  
  !----------------------------------------------------------------|
  ! Functions and subroutines                                      |
  !----------------------------------------------------------------|
  
  subroutine comp_sep

   !--------------------------------------------------------------|
   !                   Calculation portion                        |               
   !--------------------------------------------------------------|
   
   do iter = 1, dpar%ngibbs
      !--------------------- BP SWAP CHUNK -----------------------|
      ! -- Swap in a different BeyondPlanck map each iteration -- |
      !-----------------------------------------------------------|
      if (dpar%bp_swap) then
         call swap_bp_maps(ddata,dpar)
         write(*,*) ''
         bp_iter = bp_iter + 1
         call convert_bp_maps(ddata, dpar)
         write(*,*) ''
         ! ! Check to see if any swapped maps need to be dust corrected                               
         ! do j = 1, nbands
         !    if (dpar%bp_map(j)) then
         !       if (dpar%dust_corr(j)) then
         !          call dust_correct_band(ddata,dpar,dcomps,j,iter)
         !       end if
         !    end if
         ! end do
         write(*,*) ''
      end if

      ! ------------------------------------------------------------------------------------------
      ! Sample each CG group for amplitudes
      ! ------------------------------------------------------------------------------------------
      call sample_cg_groups(dpar,ddata)
      call ddata%update_sky_model
      call write_stats_to_term(ddata,dpar,iter)

      ! ------------------------------------------------------------------------------------------
      ! Sample each spectral parameter
      ! ------------------------------------------------------------------------------------------
      call sample_spectral_parameters(ddata)
      call ddata%update_sky_model
      call write_stats_to_term(ddata,dpar,iter)
      
      ! ------------------------------------------------------------------------------------------
      ! Sample band calibrators
      ! ------------------------------------------------------------------------------------------
      if (iter > 1) then
         call sample_calibrators(ddata)
         call ddata%update_sky_model
      end if

      ! ------------------------------------------------------------------------------------------
      ! Write out the data
      ! ------------------------------------------------------------------------------------------
      do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
         ! if (rank == master) then
         if (mod(iter,dpar%iter_out) .EQ. 0) then
            call write_data(dpar,ddata,k)
         end if
      end do
      if (mod(iter,dpar%iter_out) .EQ. 0) then
         call write_maps(dpar,ddata)
      end if
      write(*,*) ''
      ! ------------------------------------------------------------------------------------------
   end do
   call mpi_finalize(ierr)
 end subroutine comp_sep

 ! ------------------------------------------------------------------------------------------
 ! Specifically for the hi_fitting mode
 ! ------------------------------------------------------------------------------------------ 
 subroutine hi_fit

   do i = 1, ncomp
      if (component_list(i)%p%type == 'hi_fit') then
         call ddata%mask_hi(dpar,component_list(i)%p)
      end if
   end do

   do iter = 1, dpar%ngibbs
            
      ! ------------------------------------------------------------------------------------------
      ! Sample each CG group for amplitudes
      ! ------------------------------------------------------------------------------------------
      call sample_cg_groups(dpar,ddata)
      call ddata%update_sky_model
      call write_stats_to_term(ddata,dpar,iter)
      
      ! ------------------------------------------------------------------------------------------
      ! Sample each spectral parameter
      ! ------------------------------------------------------------------------------------------
      call sample_spectral_parameters(ddata)
      call ddata%update_sky_model
      call write_stats_to_term(ddata,dpar,iter)

      ! ------------------------------------------------------------------------------------------
      ! Write out the data
      ! ------------------------------------------------------------------------------------------
      do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
         if (rank == master) then
            call write_data(dpar,ddata,k)
         end if
      end do
      if (mod(iter,dpar%iter_out) .EQ. 0) then
         call write_maps(dpar,ddata)
      end if
      ! ------------------------------------------------------------------------------------------

   end do
 end subroutine hi_fit
end program dang
