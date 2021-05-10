program dang
  use healpix_types
  use pix_tools
  use fitstools
  use udgrade_nr
  use dang_util_mod
  use dang_param_mod
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
  real(dp), allocatable, dimension(:,:) :: map, rms, mask, true_synch

  real(dp), allocatable, dimension(:)   :: chi_pix

  real(dp)                   :: result
  character(len=6)           :: pixel_str
  
  ! Object Orient
  type(params)       :: par
  type(data)         :: dang_data
  type(component)    :: comp
  
  call init_mpi()
  call read_param_file(par)
  
  tqu(1)            = 'T'
  tqu(2)            = 'Q'
  tqu(3)            = 'U'
  
  !----------------------------------------------------------------------------------------------------------
  ! General paramters
  if (trim(par%mode) == 'comp_sep') then
     i              = getsize_fits(par%temp_file(1), nside=nside, ordering=ordering, nmaps=nmaps)
  else if (trim(par%mode) == 'hi_fit') then
     i                 = getsize_fits(par%mask_file, nside=nside, ordering=ordering, nmaps=nmaps)
  else
     write(*,*) 'Unrecognized operational mode. Do better!'
  end if
  dang_data%npix    = nside2npix(nside) 
  npix              = dang_data%npix
  nbands            = par%numinc
  nfgs              = par%ncomp+par%ntemp
  npixpar           = 0.d0
  nglobalpar        = 0.d0
  nump              = 0
  nlheader          = size(header)
  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------
  ! Array allocation - dummy maps used for loading in data
  allocate(mask(0:npix-1,1))
  allocate(map(0:npix-1,nmaps))
  allocate(rms(0:npix-1,nmaps))
  allocate(true_synch(0:npix-1,nmaps))
  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------
  ! Either do component separation, or carry out the HI fit

  call RANDOM_SEED()

  ! Initialize data and components
  if (trim(par%mode) == 'comp_sep') then
     ! Joint Sampler Info
     !----------------------------------------------------------------------------------------------------------
     call init_fg_map(dang_data,npix,nmaps,nbands,nfgs)
     call init_data_maps(dang_data,npix,nmaps,nbands)
     call init_mask_maps(dang_data,npix,nmaps)
     call init_template(dang_data,npix,nmaps,par%ntemp,nbands)
     call init_synch(comp,par,npix,nmaps)
     call init_dust(comp,par,npix,nmaps)

     dang_data%fg_map    = 0.0
     dang_data%temp_amps = 0.0
     !----------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------
     ! Read maps
     do i = 1, par%ntemp
        ! Read in templates
        call read_bintab(par%temp_file(i), dang_data%temps(:,:,i), npix,nmaps,nullval,anynull,header=header)
        ! Normalize templates
        do k = 1, nmaps
           dang_data%temp_norm(k,i) = 1.d0!maxval(dang_data%temps(:,k,i))
           dang_data%temps(:,k,i)   = dang_data%temps(:,k,i)/dang_data%temp_norm(k,i)
        end do
     end do
  
     write(*,*) 'Reading MASKFILE'
     call read_bintab(par%mask_file,dang_data%masks,npix,1,nullval,anynull,header=header)
     do i = 0, npix-1
        do j = 1, nmaps
           if (dang_data%masks(i,j) == 0.d0 .or. dang_data%masks(i,j) == missval) then
              dang_data%masks(i,j) = missval
           else 
              nump = nump + 1
           end if
        end do
     end do
     write(*,*) ''
  
     ! Read maps in
     do j = 1, nbands
        ! If not supposed to be swapped BP maps, load that map
        if (.not. par%bp_map(j)) then
           call read_bintab(trim(par%datadir) // trim(par%dat_noisefile(j)), &
                rms,npix,nmaps,nullval,anynull,header=header)
           dang_data%rms_map(:,:,j) = rms
           call read_bintab(trim(par%datadir) // trim(par%dat_mapfile(j)), &
                map,npix,nmaps,nullval,anynull,header=header)
           dang_data%sig_map(:,:,j) = map
        end if
     end do

     call read_bintab('data/test_data/new_sims/synch/synch_030_n0064_rj.fits',map,npix,nmaps,nullval,anynull,header=header)

     true_synch = map
     
     dang_data%fg_map(:,:,par%fg_ref_loc(1),1) = true_synch
     
     deallocate(map,rms)
     
     call convert_maps(par,dang_data)
  
     do j = 1, nbands
        ! Check to see if any maps need to be dust corrected
        if (par%dust_corr(j)) then
           call dust_correct_band(dang_data,par,comp,j)
        end if
     end do
     write(*,*) ''


     ! All here and below (to call comp sep) is just for testing purposes!!!
     !----------------------------------------------------------------------------------------------------------
     ! ! 1.00 sim
     ! dang_data%temp_amps(1,:,1) = 0.41957897*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(2,:,1) = 0.17704154*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(3,:,1) = 0.09161571*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(4,:,1) = 0.12402716*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(5,:,1) = 0.20367266*dang_data%temp_norm(:,1)

     ! ! 0.50 sim
     ! dang_data%temp_amps(1,:,1) = 0.20019717*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(2,:,1) = 0.0709024475*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(3,:,1) = 0.0136504706*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(4,:,1) = 7.42349435d-4*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(5,:,1) = 3.57833056d-4*dang_data%temp_norm(:,1)
 
     ! dang_data%temp_amps(1,:,1) = 0.196d-2*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(2,:,1) = 0.475d-2*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(3,:,1) = 0.000d0*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(4,:,1) = 0.283d-2*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(5,:,1) = -0.104d-01*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(6,:,1) = 0.000d0*dang_data%temp_norm(:,1)
     ! dang_data%temp_amps(7,:,1) = 0.110d-1*dang_data%temp_norm(:,1)

     ! Extrapolating solved for bands for ze cleaning
     do k = par%pol_type(1), par%pol_type(size(par%pol_type))
        call extrapolate_template(par,dang_data,comp,1,k)
        call extrapolate_foreground(par,dang_data,comp,1,k)
     end do

     !----------------------------------------------------------------------------------------------------------

     call comp_sep
  
  else if (trim(par%mode) == 'hi_fit') then
     call init_data_maps(dang_data,npix,nmaps,nbands)
     call init_mask_maps(dang_data,npix,nmaps)
     call init_template(dang_data,npix,nmaps,par%ntemp,nbands)
     call init_hi_fit(comp, par, npix)

     write(*,*) 'Reading MASKFILE'
     call read_bintab(par%mask_file,dang_data%masks,npix,1,nullval,anynull,header=header)
     do i = 0, npix-1
        if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
           dang_data%masks(i,1) = missval
        else 
           nump = nump + 1
        end if
     end do
     write(*,*) ''

     ! Read maps in
     do j = 1, nbands
        call read_bintab(trim(par%datadir) // trim(par%dat_noisefile(j)), &
             rms,npix,nmaps,nullval,anynull,header=header)
        dang_data%rms_map(:,:,j) = rms
        call read_bintab(trim(par%datadir) // trim(par%dat_mapfile(j)), &
             map,npix,nmaps,nullval,anynull,header=header)
        dang_data%sig_map(:,:,j) = map
        ! Initialize gain and offset values from parameter file
        dang_data%gain(j)   = par%init_gain(j)
        dang_data%offset(j) = 0.d0 
     end do

     do i = 0, npix-1
        if (comp%HI(i,1) > par%thresh) then
           dang_data%masks(i,1) = missval
        else if (dang_data%rms_map(i,1,1) == 0.d0) then
           dang_data%masks(i,1) = missval
        else
           dang_data%masks(i,1) = 1.d0
        end if
     end do

     par%fit_offs(:) = .false.
     
     deallocate(map,rms)

     call hi_fit 
  
  end if

contains
  
  !----------------------------------------------------------------|
  ! Functions and subroutines                                      |
  !----------------------------------------------------------------|
  
  subroutine comp_sep

    ! Count up the degrees of freedom for chisq nromalization
    do n = 1, par%ncomp
       ! Count up foregrounds in the joint sampler
       if (ANY(par%joint_comp == trim(par%fg_label(n))) .and. par%joint_sample) then
          if (par%joint_pol) then
             npixpar = npixpar + 2*nump
          else
             npixpar = npixpar + nump
          end if
       ! Count up foregrounds not included in the joint sampler
       else if (par%fg_samp_amp(n)) then
             npixpar = npixpar + size(par%pol_type)*nump
       end if
       ! Count up spectral index parameters, either fullsky or per pixel
       if (par%fg_samp_inc(n,1)) then
          if (index(par%fg_ind_region(n,1),'pix') /= 0) then
             ! Skip for now because the contribution is not trivial - esp with a mask
             ! if (par%fg_spec_joint(n,1)) then

          else if (index(par%fg_ind_region(n,1),'full') /= 0) then
             if (par%fg_spec_joint(n,1)) then
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
    do n = 1, par%ntemp
       if (ANY(par%joint_comp == trim(par%temp_label(n)))) then
          nglobalpar = nglobalpar + par%temp_nfit(n)
       else if (.not. ANY(par%joint_comp == trim(par%temp_label(n))) .or. .not. (par%joint_sample)) then
          nglobalpar = nglobalpar + size(par%pol_type)*par%temp_nfit(n)
       end if
   end do

   ! write(*,fmt='(a,i6)') 'npixpar: ', npixpar
   ! write(*,fmt='(a,i6)') 'nglobalpar: ', nglobalpar

    !--------------------------------------------------------------|
    !                   Calculation portion                        |               
    !--------------------------------------------------------------|

    do iter = 1, par%ngibbs

       !--------------------- BP SWAP CHUNK -----------------------|
       ! -- Swap in a different BeyondPlanck map each iteration -- |
       !-----------------------------------------------------------|
       if (par%bp_swap) then
          call swap_bp_maps(dang_data,par)
          write(*,*) ''
          bp_iter = bp_iter + 1
          call convert_maps_bp(par, dang_data)
          write(*,*) ''
          ! Check to see if any swapped maps need to be dust corrected                               
          do j = 1, nbands
             if ( par%bp_map(j)) then
                if (par%dust_corr(j)) then
                   call dust_correct_band(dang_data,par,comp,j)
                end if
             end if
          end do
          write(*,*) ''
       end if       
       ! --------------------------------------------------------------
       ! Joint sampling section - checks for templates and foregrounds
       ! in the joint sampling component list
       ! --------------------------------------------------------------
       if (par%joint_sample) then
          if (par%joint_pol) then
             call sample_joint_amp(par,dang_data,comp,2,trim(par%solver))
          else
             do k = par%pol_type(1), par%pol_type(size(par%pol_type))
                call sample_joint_amp(par,dang_data,comp,k,trim(par%solver))
             end do
          end if             
          do k = par%pol_type(1), par%pol_type(size(par%pol_type))
             do m = 1, size(par%joint_comp)
                ! Extrapolate foreround solutions
                do n = 1, par%ncomp
                   if (trim(par%joint_comp(m)) == trim(par%fg_label(n))) then
                      call extrapolate_foreground(par,dang_data,comp,n,k)
                   end if
                end do
                ! Extrapolate template solutions
                do n = 1, par%ntemp
                   if (trim(par%joint_comp(m)) == trim(par%temp_label(n))) then
                      call extrapolate_template(par,dang_data,comp,n,k)
                   end if
                end do
             end do
          end do
          ! How good is the fit and what are the parameters looking like?
          call write_stats_to_term(par,dang_data,comp,iter)
      end if

       ! ------------------------------------------------------------------------------------------
       ! Sample amplitudes
       ! ------------------------------------------------------------------------------------------
       do n = 1, par%ncomp
          if (par%fg_samp_amp(n)) then
             write(*,*) "Sample "//trim(par%fg_label(n))//" amplitudes."
             do k = par%pol_type(1), par%pol_type(size(par%pol_type))
                dang_data%fg_map(:,k,par%fg_ref_loc(n),n) = sample_fg_amp(par,dang_data,comp,n,k)
                call extrapolate_foreground(par,dang_data,comp,n,k)
             end do
          end if
       end do
       do n = 1, par%ntemp
          if (.not. ANY(par%joint_comp == trim(par%temp_label(n))) .or. .not. (par%joint_sample)) then
             if (par%temp_sample(n)) then
                write(*,*) "Sample "//trim(par%fg_label(n))//" template amplitudes."
                do k = par%pol_type(1), par%pol_type(size(par%pol_type))
                   call template_fit(par,dang_data,comp,k,n)
                   call extrapolate_template(par,dang_data,comp,n,k)
                end do
             end if
          end if
       end do

       ! ------------------------------------------------------------------------------------------
       ! Sample spectral parameters
       ! ------------------------------------------------------------------------------------------
       do n = 1, par%ncomp
          if (par%fg_samp_inc(n,1)) then
             if (par%fg_spec_joint(n,1)) then
                write(*,*) "Sample "//trim(par%fg_label(n))//" beta jointly."
                write(*,*) "---------------------------------"
                call sample_index(par,dang_data,comp,n,-1)
             else
                do k = par%pol_type(1), par%pol_type(size(par%pol_type))
                write(*,*) "Sample "//trim(par%fg_label(n))//" beta for "//trim(tqu(k))//"."
                write(*,*) "---------------------------------"
                   call sample_index(par,dang_data,comp,n,k)
                end do
             end if
             do k = par%pol_type(1), par%pol_type(size(par%pol_type))
                call extrapolate_foreground(par,dang_data,comp,n,k)
             end do
             ! How good is the fit and what are the parameters looking like?
             call write_stats_to_term(par,dang_data,comp,iter)
          end if
       end do
       ! ------------------------------------------------------------------------------------------
       dang_data%res_map = dang_data%sig_map
       do k = par%pol_type(1), par%pol_type(size(par%pol_type))
          do j = 1, nfgs
             dang_data%res_map(:,k,:)  = dang_data%res_map(:,k,:) - dang_data%fg_map(:,k,:,j)
          end do
          if (rank == master) then
             call write_data(par,dang_data,comp,k)
          end if
       end do

       !-------------------------------------------------------------------------------------------
       ! write(iter_str, '(i0.5)') iter
       ! title = trim(par%outdir) // 'samp_minus_true_' // trim(iter_str) // '.dat'
       ! inquire(file=title,exist=exist)
       ! if (exist) then
       !    open(50,file=title, status="old",position="append", action="write")
       ! else
       !    open(50,file=title, status="new", action="write")
       ! endif
       ! do i = 0, npix-1
       !    if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
       !       cycle
       !    else
       !       write(50,fmt='(3(f16.8))') dang_data%fg_map(i,:,par%fg_ref_loc(1),1) - true_synch(i,:)
       !    end if
       ! end do
       ! close(50)

       !-------------------------------------------------------------------------------------------

       ! How good is the fit and what are the parameters looking like?
       call write_stats_to_term(par,dang_data,comp,iter)

       if (mod(iter,par%iter_out) .EQ. 0) then
          call write_maps(par,dang_data,comp)
       end if
    end do
    call mpi_finalize(ierr)
  end subroutine comp_sep

  ! ------------------------------------------------------------------------------------------
  ! Specifically for the hi_fitting mode
  ! ------------------------------------------------------------------------------------------ 
  subroutine hi_fit
    do iter = 1, par%ngibbs

       if (iter > 1) then
          write(*,*) 'Calc HI gain offset'
          do j = 1, nbands
             if (par%fit_gain(j)) then
                call sample_band_gain(par, dang_data, comp, 1, j, 1, 1)
             end if
             if (par%fit_offs(j)) then
                call sample_band_offset(par, dang_data, comp, 1, j, 1)
             end if
             ! if (par%fit_gain(j) .and. par%fit_offs(j)) then
             !    call calc_HI_gain_offset(par, dang_data, comp, 1, 1, j)
             ! end if
          end do
          ! write(*,"(12x,8(A16))") par%dat_label
          ! write(*,"(a8,8(E16.4))") 'Gain: ',dang_data%gain
          ! write(*,"(a8,8(E16.4))") 'Offset: ',dang_data%offset
       end if

       write(*,*) 'template_fit'
       call template_fit(par, dang_data, comp, 1)

       write(*,fmt='(i6, a, E10.3, a, e10.3)')&
            iter, " - chisq: " , dang_data%chisq, " - T_d: ",&
            mask_avg(comp%T_d(:,1),dang_data%masks(:,1))
       write(*,fmt='(a)') '---------------------------------------------'

       call sample_HI_T(par, dang_data, comp, 1)

       do j = 1, nbands
          do i = 0, npix-1
             if (dang_data%masks(i,1) == missval .or. dang_data%masks(i,1) == 0.d0) cycle
             dang_data%res_map(i,1,j) = (dang_data%sig_map(i,1,j)-dang_data%offset(j))/dang_data%gain(j) &
                  - comp%HI_amps(j)*comp%HI(i,1)*planck(par%dat_nu(j)*1d9,comp%T_d(i,1))
          end do
       end do

       call compute_chisq(par,dang_data,comp,1)

       if (rank == master) then
          if (mod(iter, 1) == 0 .or. iter == 1) then
             if (nbands .lt. 10) then
                write(*,fmt='(i6, a, E10.3, a, e10.3, a, 10e10.3)')&
                     iter, " - chisq: " , dang_data%chisq, " - T_d: ",&
                     mask_avg(comp%T_d(:,1),dang_data%masks(:,1)), ' - A_HI: ', comp%HI_amps
                write(*,fmt='(a)') '---------------------------------------------'
             else
                write(*,fmt='(i6, a, E10.3, a, e10.3)')&
                     iter, " - chisq: " , dang_data%chisq, " - T_d: ",&
                     mask_avg(comp%T_d(:,1),dang_data%masks(:,1))
                write(*,fmt='(a)') '---------------------------------------------'
             end if
          end if
          if (mod(iter,par%iter_out) .EQ. 0) then
             call write_maps(par,dang_data,comp)
          end if
          call write_data(par,dang_data,comp,1)
       end if
       write(*,*) ''
    end do
  end subroutine hi_fit
end program dang
