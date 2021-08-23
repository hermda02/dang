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
  integer(i4b)                          :: bp_iter
  real(dp), allocatable, dimension(:,:) :: map, rms, true_synch
      
  ! Object Orient
  type(dang_params) :: dpar
  type(dang_data)   :: ddata
  type(dang_comps)  :: dcomps
  ! type(dang_bp)     :: dbp 

  ! allocate(dbp(dpar%numinc))

  call init_mpi()
  call read_param_file(dpar)
  call init_bp_mod(dpar)

  ! write(*,*) dpar%fg_nu_ref(1)
  ! do i = 1, dpar%numinc
  !    write(*,*) bp(i)%nu_c
  ! end do
  ! stop

  !----------------------------------------------------------------------------------------------------------
  ! General paramters
  if (trim(dpar%mode) == 'comp_sep') then
     i = getsize_fits(dpar%temp_file(1), nside=nside, ordering=ordering, nmaps=nmaps)
  else if (trim(dpar%mode) == 'hi_fit') then
     i = getsize_fits(dpar%mask_file,    nside=nside, ordering=ordering, nmaps=nmaps)
  else
     write(*,*) 'Unrecognized operational mode. Do better!'
  end if
  ddata%npix = nside2npix(nside) 
  npix       = ddata%npix
  nbands     = dpar%numinc
  nfgs       = dpar%ncomp+dpar%ntemp
  npixpar    = 0.d0
  nglobalpar = 0.d0
  nump       = 0
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

     ! Initialize component
     call init_synch(dcomps,dpar,npix,nmaps)
     call init_dust(dcomps,dpar,npix,nmaps)

     !----------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------
     ! Read maps
     call ddata%read_data_maps(dpar)
     
     call convert_maps(ddata,dpar)
     do j = 1, nbands
        ! Check to see if any maps need to be dust corrected
        if (dpar%dust_corr(j) .and. .not. dpar%bp_map(j)) then
           call dust_correct_band(ddata,dpar,dcomps,j)
        end if
     end do
     write(*,*) ''
     call comp_sep
  
  !------------------------------------------------------------------------------------------------
  else if (trim(dpar%mode) == 'hi_fit') then
     ! Initialize ddata and components
     !----------------------------------------------------------------------------------------------------------
     call ddata%init_data_maps(dpar)
     call init_hi_fit(dcomps, dpar, npix)

     call read_data_maps(ddata,dpar)

     do i = 0, npix-1
        if (dcomps%HI(i,1) > dpar%thresh) then
           ddata%masks(i,1) = missval
        else if (ddata%masks(i,1) == missval) then
           ddata%masks(i,1) = missval
        else if (ddata%rms_map(i,1,1) == 0.d0) then
           ddata%masks(i,1) = missval
        else
           ddata%masks(i,1) = 1.d0
        end if
     end do
     nump = 0
     do i = 0, npix-1
        do j = 1, nmaps
           if (ddata%masks(i,j) == 0.d0 .or. ddata%masks(i,j) == missval) then
              ddata%masks(i,j) = missval
           else 
              nump = nump + 1
           end if
        end do
     end do
     
     call hi_fit 
  
  end if

contains
  
  !----------------------------------------------------------------|
  ! Functions and subroutines                                      |
  !----------------------------------------------------------------|
  
  subroutine comp_sep

    ! Count up the degrees of freedom for chisq nromalization
    do n = 1, dpar%ncomp
       ! Count up foregrounds in the joint sampler
       if (ANY(dpar%joint_comp == trim(dpar%fg_label(n))) .and. dpar%joint_sample) then
          if (dpar%joint_pol) then
             npixpar = npixpar + 2*nump
          else
             npixpar = npixpar + nump
          end if
       ! Count up foregrounds not included in the joint sampler
       else if (dpar%fg_samp_amp(n)) then
             npixpar = npixpar + size(dpar%pol_type)*nump
       end if
       ! Count up spectral index parameters, either fullsky or per pixel
       if (dpar%fg_samp_spec(n,1)) then
          if (index(dpar%fg_ind_region(n,1),'pix') /= 0) then
             ! Skip for now because the contribution is not trivial - esp with a mask
             ! if (dpar%fg_spec_joint(n,1)) then

          else if (index(dpar%fg_ind_region(n,1),'full') /= 0) then
             if (dpar%fg_spec_joint(n,1)) then
                ! Sampling fullsky jointly in Q and U
                nglobalpar = nglobalpar + 1
             else 
                ! Sampling fullsky separately in Q and U
                nglobalpar = nglobalpar + 2
             end if
          end if
       end if
    end do
    ! 
    do n = 1, dpar%ntemp
       if (ANY(dpar%joint_comp == trim(dpar%temp_label(n)))) then
          nglobalpar = nglobalpar + dpar%temp_nfit(n)
       else if (.not. ANY(dpar%joint_comp == trim(dpar%temp_label(n))) .or. .not. (dpar%joint_sample)) then
          nglobalpar = nglobalpar + size(dpar%pol_type)*dpar%temp_nfit(n)
       end if
   end do

   ! write(*,fmt='(a,i6)') 'npixpar: ', npixpar
   ! write(*,fmt='(a,i6)') 'nglobalpar: ', nglobalpar

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
          call convert_maps_bp(ddata, dpar)
          write(*,*) ''
          ! Check to see if any swapped maps need to be dust corrected                               
          do j = 1, nbands
             if (dpar%bp_map(j)) then
                if (dpar%dust_corr(j)) then
                   call dust_correct_band(ddata,dpar,dcomps,j,iter)
                end if
             end if
          end do
          write(*,*) ''
       end if       
       ! --------------------------------------------------------------
       ! Joint sampling section - checks for templates and foregrounds
       ! in the joint sampling component list
       ! --------------------------------------------------------------
       if (dpar%joint_sample) then
          if (dpar%joint_pol) then
             call sample_joint_amp(dpar,ddata,dcomps,2,trim(dpar%solver))
          else
             do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
                call sample_joint_amp(dpar,ddata,dcomps,k,trim(dpar%solver))
             end do
          end if             
          do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
             do m = 1, size(dpar%joint_comp)
                ! Extrapolate foreround solutions
                do n = 1, dpar%ncomp
                   if (trim(dpar%joint_comp(m)) == trim(dpar%fg_label(n))) then
                      call extrapolate_foreground(dpar,ddata,dcomps,n,k)
                   end if
                end do
                ! Extrapolate template solutions
                do n = 1, dpar%ntemp
                   if (trim(dpar%joint_comp(m)) == trim(dpar%temp_label(n))) then
                      call extrapolate_template(dpar,ddata,dcomps,n,k)
                   end if
                end do
             end do
          end do
          ! How good is the fit and what are the parameters looking like?
          call write_stats_to_term(ddata,dpar,dcomps,iter)
      end if

       ! ------------------------------------------------------------------------------------------
       ! Sample amplitudes
       ! ------------------------------------------------------------------------------------------
       do n = 1, dpar%ncomp
          if (dpar%fg_samp_amp(n)) then
             write(*,*) "Sample "//trim(dpar%fg_label(n))//" amplitudes."
             do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
                ddata%fg_map(:,k,dpar%fg_ref_loc(n),n) = sample_fg_amp(dpar,ddata,dcomps,n,k)
                call extrapolate_foreground(dpar,ddata,dcomps,n,k)
             end do
          end if
       end do
       do n = 1, dpar%ntemp
          if (.not. ANY(dpar%joint_comp == trim(dpar%temp_label(n))) .or. .not. (dpar%joint_sample)) then
             if (dpar%temp_sample(n)) then
                write(*,*) "Sample "//trim(dpar%fg_label(n))//" template amplitudes."
                do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
                   call template_fit(dpar,ddata,dcomps,k,n)
                   call extrapolate_template(dpar,ddata,dcomps,n,k)
                end do
             end if
          end if
       end do

       ! ------------------------------------------------------------------------------------------
       ! Sample spectral parameters
       ! ------------------------------------------------------------------------------------------
       do n = 1, dpar%ncomp
          if (dpar%fg_samp_spec(n,1)) then
             if (dpar%fg_spec_joint(n,1)) then
                write(*,*) "Sample "//trim(dpar%fg_label(n))//" beta jointly."
                write(*,*) "---------------------------------"
                ! call sample_new_index(dpar,ddata,dcomps,n,-1)
                call sample_index(dpar,ddata,dcomps,n,-1)
             else
                do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
                write(*,*) "Sample "//trim(dpar%fg_label(n))//" beta for "//trim(tqu(k))//"."
                write(*,*) "---------------------------------"
                   call sample_index(dpar,ddata,dcomps,n,k)
                end do
             end if
             do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
                call extrapolate_foreground(dpar,ddata,dcomps,n,k)
             end do
             ! How good is the fit and what are the parameters looking like?
             call write_stats_to_term(ddata,dpar,dcomps,iter)
          end if
       end do
       ! ------------------------------------------------------------------------------------------
       ddata%res_map = ddata%sig_map
       do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
          do j = 1, nfgs
             ddata%res_map(:,k,:)  = ddata%res_map(:,k,:) - ddata%fg_map(:,k,:,j)
          end do
          if (rank == master) then
             call write_data(dpar,ddata,dcomps,k)
          end if
       end do
       !-------------------------------------------------------------------------------------------

       ! How good is the fit and what are the parameters looking like?
       call write_stats_to_term(ddata,dpar,dcomps,iter)

       if (mod(iter,dpar%iter_out) .EQ. 0) then
          call write_maps(dpar,ddata,dcomps)
       end if
    end do
    call mpi_finalize(ierr)
  end subroutine comp_sep

  ! ------------------------------------------------------------------------------------------
  ! Specifically for the hi_fitting mode
  ! ------------------------------------------------------------------------------------------ 
  subroutine hi_fit

    npixpar    = nump
    nglobalpar = nbands

    do iter = 1, dpar%ngibbs

       if (iter > 1) then
          do j = 1, nbands
             if (dpar%fit_gain(j)) then
                call sample_band_gain(dpar, ddata, dcomps, 1, j, 1, 1)
             end if
             if (dpar%fit_offs(j)) then
                call sample_band_offset(dpar, ddata, dcomps, 1, j, 1)
             end if
          end do
       end if

       write(*,*) 'Fit templates'
       call template_fit(dpar, ddata, dcomps, 1)
       call write_stats_to_term(ddata,dpar,dcomps,iter)
       write(*,*) 'Sample Td'
       call sample_HI_T(dpar, ddata, dcomps, 1)

       do j = 1, nbands
          do i = 0, npix-1
             if (ddata%masks(i,1) == missval .or. ddata%masks(i,1) == 0.d0) cycle
             ddata%res_map(i,1,j) = (ddata%sig_map(i,1,j)-ddata%offset(j))/ddata%gain(j) &
                  - dcomps%HI_amps(j)*dcomps%HI(i,1)*planck(dpar%band_nu(j)*1d9,dcomps%T_d(i,1))
          end do
       end do

       ! How good is the fit and what are the parameters looking like?
       call write_stats_to_term(ddata,dpar,dcomps,iter)

       if (mod(iter,dpar%iter_out) .EQ. 0) then
          call write_maps(dpar,ddata,dcomps)
       end if
       call write_data(dpar,ddata,dcomps,1)
       write(*,*) ''
    end do
  end subroutine hi_fit
end program dang
