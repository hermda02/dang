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
  i = getsize_fits(dpar%mask_file, nside=nside, ordering=ordering, nmaps=nmaps)
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
  verbosity  = 1
  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------

  call RANDOM_SEED()

  ! Initialize ddata and components
  !----------------------------------------------------------------------------------------------------------
  call initialize_components(dpar)
  call ddata%initialize_data_module(dpar)
  call initialize_cg_groups(dpar)
  if (mask_hi) call ddata%mask_hi_threshold(dpar)
  call ddata%update_sky_model
  call write_maps(ddata,dpar)     
  do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
     call ddata%write_data(dpar, k)
  end do
  write(*,*) '---------------------------'
  write(*,*) ' Starting main Gibbs Chain '
  write(*,*) '---------------------------'

  !--------------------------------------------------------------|
  !                   Computation portion                        |               
  !--------------------------------------------------------------|
  do iter = 1, dpar%ngibbs
      t0 = mpi_wtime()
     !--------------------- BP SWAP CHUNK -----------------------|
     ! -- Swap in a different BeyondPlanck map each iteration -- |
     !-----------------------------------------------------------|
     if (dpar%cg_swap .and. iter .ne. 1) then
        call swap_cg_maps(ddata, dpar)
        write(*,*) ''
        call convert_cg_maps(ddata)
        write(*,*) ''
     end if
     ! ------------------------------------------------------------------------------------------
     ! Sample each CG group for amplitudes
     ! ------------------------------------------------------------------------------------------
     call sample_cg_groups(dpar,ddata)
     if (iter > 1) then
        ! ------------------------------------------------------------------------------------------
        ! Sample each spectral parameter
        ! ------------------------------------------------------------------------------------------
        call sample_spectral_parameters(dpar,ddata)
        ! ------------------------------------------------------------------------------------------
        ! Sample band calibrators
        ! ------------------------------------------------------------------------------------------
        call sample_calibrators(ddata)
     end if
     
     ! ------------------------------------------------------------------------------------------
     ! Write out the data
     ! ------------------------------------------------------------------------------------------
     do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
        call ddata%write_data(dpar, k)
     end do
     if (mod(iter,dpar%iter_out) .EQ. 0) then
        call ddata%write_maps(dpar)
     end if
     write(*,*) ''
     t00 = mpi_wtime()
     write(*,fmt='(a,f12.5)') 'Total CPU time = ',t00-t0
     ! ------------------------------------------------------------------------------------------
  end do
  call mpi_finalize(ierr)

end program dang
