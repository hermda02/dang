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
  real(dp), allocatable, dimension(:,:) :: map, rms, mask
  
  ! Object Orient
  type(params)       :: par
  type(data)         :: dang_data
  type(component)    :: comp
  
  call init_mpi()
  call read_param_file(par)
  
  allocate(comp%joint(2))
  comp%joint(1) = 'synch'
  comp%joint(2) = 'template01'
  
  !----------------------------------------------------------------------------------------------------------
  ! General paramters
  tqu(1)            = 'T'
  tqu(2)            = 'Q'
  tqu(3)            = 'U'
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
  npar              = 3
  nump              = 0
  nlheader          = size(header)
  !----------------------------------------------------------------------------------------------------------
  ! Array allocation - dummy maps used for loading in data
  allocate(mask(0:npix-1,1))
  allocate(map(0:npix-1,nmaps))
  allocate(rms(0:npix-1,nmaps))
  !----------------------------------------------------------------------------------------------------------

  call RANDOM_SEED()

  ! Either do component separation, or carry out the HI fit
  !------------------------------------------------------------------------------------------------
  if (trim(par%mode) == 'comp_sep') then
     ! Initialize data and components
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
           call read_bintab(trim(par%datadir) // trim(par%band_noisefile(j)), &
                rms,npix,nmaps,nullval,anynull,header=header)
           dang_data%rms_map(:,:,j) = rms
           call read_bintab(trim(par%datadir) // trim(par%band_mapfile(j)), &
                map,npix,nmaps,nullval,anynull,header=header)
           dang_data%sig_map(:,:,j) = map
        end if
     end do
     
     deallocate(map,rms)
     
     call convert_maps(par,dang_data)
  
     do j = 1, nbands
        ! Check to see if any maps need to be dust corrected
        if (par%dust_corr(j)) then
           call dust_correct_band(dang_data,par,comp,j)
        end if
     end do
     write(*,*) ''

     call comp_sep
  
  !------------------------------------------------------------------------------------------------
  else if (trim(par%mode) == 'hi_fit') then
     ! Initialize data and components
     !----------------------------------------------------------------------------------------------------------
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

     k = 1 
     ! Read maps in
     do j = 1, par%numband
        if (par%band_inc(j)) then
           call read_bintab(trim(par%datadir) // trim(par%band_noisefile(k)), &
                rms,npix,nmaps,nullval,anynull,header=header)
           dang_data%rms_map(:,:,k) = rms
           call read_bintab(trim(par%datadir) // trim(par%band_mapfile(k)), &
                map,npix,nmaps,nullval,anynull,header=header)
           dang_data%sig_map(:,:,k) = map
           ! Initialize gain and offset values from parameter file
           dang_data%gain(k)   = par%init_gain(k)
           dang_data%offset(k) = 0.d0 
           k = k + 1
        end if
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
       if (par%joint_sample) then
          call sample_joint_amp(par,dang_data,comp,2,trim(par%solver))
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
          do k = par%pol_type(1), par%pol_type(size(par%pol_type))
             call compute_chisq(k,chisq,par%mode)
             if (rank == master) then
                if (mod(iter, 1) == 0 .or. iter == 1) then
                   write(*,fmt='(i6, a, E10.3, a, f7.3, a, f8.4, a, 10e10.3)')&
                        iter, " - chisq: " , chisq, " - A_s: ",&
                        dang_data%fg_map(23000,k,par%fg_ref_loc(1),1),  " - beta_s: ",&
                        mask_avg(comp%beta_s(:,k),dang_data%masks(:,1)), ' - A_d: ', &
                        dang_data%temp_amps(:,k,1)/dang_data%temp_norm(k,1)
                   write(*,fmt='(a)') '---------------------------------------------'
                end if
             end if
          end do
       end if

       ! ------------------------------------------------------------------------------------------
       ! Sample amplitudes
       ! ------------------------------------------------------------------------------------------
       do n = 1, par%ncomp
          if (par%fg_samp_amp(n)) then
             do k = par%pol_type(1), par%pol_type(size(par%pol_type))
                dang_data%fg_map(:,k,par%fg_ref_loc(n),n) = sample_spec_amp(par,dang_data,comp,n,k)
                call extrapolate_foreground(par,dang_data,comp,n,k)
             end do
          end if
       end do
       ! do n = 1, par%ntemp
       !    if (.not. ANY(par%joint_comp == trim(par%temp_label(n)))) then
       !       if (par%fg_amp_joint(n)) then
       !          call fit_template()
       !       end if
       !    end if
       !    call extrapolate_template(par,dang_data,comp,n,k)
       ! end do

       ! ------------------------------------------------------------------------------------------
       ! Sample spectral parameters
       ! ------------------------------------------------------------------------------------------
       do n = 1, par%ncomp
          if (par%fg_samp_spec(n,1)) then
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
             do k = par%pol_type(1), par%pol_type(size(par%pol_type))
                call compute_chisq(k,chisq,par%mode)
                if (rank == master) then
                   if (mod(iter, 1) == 0 .or. iter == 1) then
                      write(*,fmt='(i6, a, E10.3, a, f7.3, a, f8.4, a, 10e10.3)')&
                           iter, " - chisq: " , chisq, " - A_s: ",&
                           dang_data%fg_map(23000,k,par%fg_ref_loc(1),1),  " - beta_s: ",&
                           mask_avg(comp%beta_s(:,k),dang_data%masks(:,1)), ' - A_d: ', &
                           dang_data%temp_amps(:,k,1)/dang_data%temp_norm(k,1)
                      write(*,fmt='(a)') '---------------------------------------------'
                   end if
                end if
             end do
          end if
       end do
       ! ------------------------------------------------------------------------------------------
       dang_data%res_map = dang_data%sig_map
       do k = par%pol_type(1), par%pol_type(size(par%pol_type))
          do j = 1, nfgs
             dang_data%res_map(:,k,:)  = dang_data%res_map(:,k,:) - dang_data%fg_map(:,k,:,j)
          end do
          
          call compute_chisq(k,chisq,par%mode)
          
          if (rank == master) then
             if (mod(iter, 1) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, E10.3, a, f7.3, a, f8.4, a, 10e10.3)')&
                     iter, " - chisq: " , chisq, " - A_s: ",&
                     dang_data%fg_map(23000,k,par%fg_ref_loc(1),1),  " - beta_s: ",&
                     mask_avg(comp%beta_s(:,k),dang_data%masks(:,1)), ' - A_d: ', &
                     dang_data%temp_amps(:,k,1)/dang_data%temp_norm(k,1)
                write(*,fmt='(a)') '---------------------------------------------'
             end if
             call write_data(par%mode)
          end if
       end do
       if (mod(iter,par%iter_out) .EQ. 0) then
          call write_maps(par%mode)
       end if
       write(*,*) ''
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
          end do
       end if

       write(*,*) 'template_fit'
       call template_fit(par, dang_data, comp, 1)

       write(*,fmt='(i6, a, E10.3, a, e10.3)')&
            iter, " - chisq: " , chisq, " - T_d: ",&
            mask_avg(comp%T_d(:,1),dang_data%masks(:,1))
       write(*,fmt='(a)') '---------------------------------------------'

       call sample_HI_T(par, dang_data, comp, 1)

       do j = 1, nbands
          do i = 0, npix-1
             if (dang_data%masks(i,1) == missval .or. dang_data%masks(i,1) == 0.d0) cycle
             dang_data%res_map(i,1,j) = (dang_data%sig_map(i,1,j)-dang_data%offset(j))/dang_data%gain(j) &
                  - comp%HI_amps(j)*comp%HI(i,1)*planck(par%band_nu(j)*1d9,comp%T_d(i,1))
          end do
       end do

       call compute_chisq(1,chisq,par%mode)

       if (rank == master) then
          if (mod(iter, 1) == 0 .or. iter == 1) then
             if (nbands .lt. 10) then
                write(*,fmt='(i6, a, E10.3, a, e10.3, a, 10e10.3)')&
                     iter, " - chisq: " , chisq, " - T_d: ",&
                     mask_avg(comp%T_d(:,1),dang_data%masks(:,1)), ' - A_HI: ', comp%HI_amps
                write(*,fmt='(a)') '---------------------------------------------'
             else
                write(*,fmt='(i6, a, E10.3, a, e10.3)')&
                     iter, " - chisq: " , chisq, " - T_d: ",&
                     mask_avg(comp%T_d(:,1),dang_data%masks(:,1))
                write(*,fmt='(a)') '---------------------------------------------'
             end if
          end if
          if (mod(iter,par%iter_out) .EQ. 0) then
             call write_maps(par%mode)
          end if
          call write_data(par%mode)
       end if
       write(*,*) ''
    end do
  end subroutine hi_fit
  
  subroutine write_maps(mode)
    implicit none
    
    character(len=16), intent(in)       :: mode
    real(dp), dimension(0:npix-1,nmaps) :: map
    real(dp)                            :: s, signal
    integer(i4b)                        :: n, mn
    character(len=128)                  :: title

    write(*,*) 'Output data maps'
    
    if (trim(mode) == 'comp_sep') then
       
       write(iter_str, '(i0.5)') iter
       if (par%output_fg .eqv. .true.) then
          do j = 1, nbands
             do n = 1, par%ntemp
                title = trim(par%outdir) // trim(par%band_label(j)) //'_'// trim(par%temp_label(n)) //&
                     '_k' // trim(iter_str) // '.fits'
                map(:,:)   = dang_data%fg_map(:,:,j,n+par%ncomp)
                do i = 0, npix-1
                   if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
                      map(i,:) = missval
                   end if
                end do
                call write_result_map(trim(title), nside, ordering, header, map)
             end do
             do n = 1, par%ncomp
                title = trim(par%outdir) // trim(par%band_label(j)) //'_'// trim(par%fg_label(n)) //&
                     '_amplitude_k' // trim(iter_str) // '.fits'
                map(:,:)   = dang_data%fg_map(:,:,j,n)
                do i = 0, npix-1
                   if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
                      map(i,1) = missval
                   end if
                end do
                call write_result_map(trim(title), nside, ordering, header, map)
             end do
          end do
       else 
          do n = 1, par%ncomp
             title = trim(par%outdir) // trim(par%band_label(par%fg_ref_loc(n))) //'_'// trim(par%fg_label(n)) //&
                  '_amplitude_k' // trim(iter_str) // '.fits'
             map(:,:)   = dang_data%fg_map(:,:,par%fg_ref_loc(n),n)
             do i = 0, npix-1
                if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
                   map(i,1) = missval
                end if
             end do
             call write_result_map(trim(title), nside, ordering, header, map)
          end do
       end if
       do j = 1, nbands
          title = trim(par%outdir) // trim(par%band_label(j)) // '_residual_k' // trim(iter_str) // '.fits'
          map(:,:)   = dang_data%res_map(:,:,j)
          do i = 0, npix-1
             if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
                map(i,1) = missval
             end if
          end do
          call write_result_map(trim(title), nside, ordering, header, map)
       end do
       title = trim(par%outdir) // 'synch_beta_k' // trim(iter_str) // '.fits'
       map(:,:)   = comp%beta_s(:,:)
       do i = 0, npix-1
          if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
             map(i,1) = missval
          end if
       end do
       call write_result_map(trim(title), nside, ordering, header, map)
       do mn = 1, nmaps
          dang_data%chi_map(:,mn) = 0.d0
          do i = 0, npix-1
             do j = 1, nbands
                s      = 0.d0
                do l = 1, nfgs
                   signal = dang_data%fg_map(i,mn,j,l)
                   s      = s + signal
                end do
                dang_data%chi_map(i,mn) = dang_data%chi_map(i,mn) + dang_data%masks(i,1)*&
                     (dang_data%sig_map(i,mn,j) - s)**2.d0/dang_data%rms_map(i,mn,j)**2.d0
             end do
          end do
       end do
       dang_data%chi_map(:,:) = dang_data%chi_map(:,:)/(nbands+nfgs)
       title = trim(par%outdir) // 'chisq_k'// trim(iter_str) // '.fits'
       map(:,:)   = dang_data%chi_map(:,:)
       do i = 0, npix-1
          if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
             map(i,1) = missval
          end if
       end do
       call write_result_map(trim(title), nside, ordering, header, map)
    else if (trim(mode) == 'hi_fit') then

       write(iter_str, '(i0.5)') iter
       n = 1
       do j = 1, par%numband
          if (par%band_inc(j)) then
             title = trim(par%outdir) // trim(par%band_label(j)) //'_hi_amplitude_k'// trim(iter_str) // '.fits'
             map(:,1)   = comp%HI_amps(n)*comp%HI(:,1)
             do i = 0, npix-1
                if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
                   map(i,1) = missval
                end if
             end do
             call write_bintab(map,npix,1, header, nlheader, trim(title))
             n = n + 1
          end if
       end do

       n = 1
       do j = 1, par%numband
          if (par%band_inc(j)) then
             title = trim(par%outdir) // trim(par%band_label(j)) // '_residual_k' // trim(iter_str) // '.fits'
             map(:,:)   = dang_data%res_map(:,:,n)
             do i = 0, npix-1
                if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
                   map(i,1) = missval
                end if
             end do
             call write_bintab(map,npix,1, header, nlheader, trim(title))
             n = n + 1
          end if
       end do
       title = trim(par%outdir) // 'T_d_k'// trim(iter_str) // '.fits'
       map(:,1)   = comp%T_d(:,1)
       do i = 0, npix-1
          if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
             map(i,1) = missval
          end if
       end do
       call write_bintab(map,npix,1, header, nlheader, trim(title))
       dang_data%chi_map = 0.d0
       do i = 0, npix-1
          do j = 1, nbands
             s =  comp%HI_amps(j)*comp%HI(i,1)*planck(par%band_nu(j)*1d9,comp%T_d(i,1))
             dang_data%chi_map(i,1) = dang_data%chi_map(i,1) + dang_data%masks(i,1)*&
                  (dang_data%sig_map(i,1,j) - s)**2.d0/dang_data%rms_map(i,1,j)**2.d0
          end do
       end do
       dang_data%chi_map(:,1) = dang_data%chi_map(:,1)/(nbands+nfgs)
       title = trim(par%outdir) // 'chisq_k' // trim(iter_str) // '.fits'
       map(:,1)   = dang_data%chi_map(:,1)
       do i = 0, npix-1
          if (dang_data%masks(i,1) == 0.d0 .or. dang_data%masks(i,1) == missval) then
             map(i,1) = missval
          end if
       end do
       call write_bintab(map,npix,1, header, nlheader, trim(title))
          
    end if
    
  end subroutine write_maps
  
  subroutine write_data(mode)
    implicit none
    character(len=16), intent(in) :: mode
    character(len=2)              :: temp_n
    character(len=128)            :: title

    if (trim(mode) == 'comp_sep') then
       
       title = trim(par%outdir) // 'pixel_23000_A_d_' // trim(tqu(k)) // '.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(30,file=title, status="old",position="append", action="write")
       else
          open(30,file=title, status="new", action="write")
          write(30,*) 
       endif
       write(30,*) dang_data%fg_map(23000,k,par%fg_ref_loc(1),2)
       close(30)
       
       title = trim(par%outdir) // 'pixel_23000_A_s_' // trim(tqu(k)) // '.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(31,file=title, status="old",position="append", action="write")
       else
          open(31,file=title, status="new", action="write")
       endif
       write(31,*) dang_data%fg_map(23000,k,par%fg_ref_loc(1),1)
       close(31)
       
       title = trim(par%outdir) // 'pixel_23000_beta_s_' // trim(tqu(k)) // '.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(32,file=title, status="old",position="append", action="write")
       else
          open(32,file=title, status="new", action="write")
       endif
       write(32,*) comp%beta_s(23000,k)
       close(32)
       
       title = trim(par%outdir) // 'total_chisq_' // trim(tqu(k)) // '.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(33,file=title, status="old",position="append", action="write")
       else
          open(33,file=title, status="new", action="write")
       endif
       call compute_chisq(k,chisq,par%mode)
       write(33,*) chisq
       close(33)
       
       do i = 1, par%ntemp
          title = trim(par%outdir) //  trim(par%temp_label(i)) // '_' //trim(tqu(k)) // '_amplitudes.dat'
          inquire(file=title,exist=exist)
          if (exist) then
             open(34,file=title, status="old", &
                  position="append", action="write")
          else
             open(34,file=title, status="new", action="write")
          endif
          write(34,'(10(E17.8))') dang_data%temp_amps(:,k,i)
          close(34)
       end do
       
    else if (trim(mode) == 'hi_fit') then
       title = trim(par%outdir) // 'HI_amplitudes.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(35,file=title,status="old",position="append",action="write") 
       else
          open(35,file=title,status="new",action="write")
       end if
       write(35,'(10(E17.8))') comp%HI_amps
       close(35)
       
       title = trim(par%outdir)//'HI_chisq.dat'
       inquire(file=title,exist=exist)
       if (exist) then
          open(36,file=title,status="old",position="append",action="write") 
       else
          open(36,file=title,status="new",action="write")
       end if
       call compute_chisq(1,chisq,mode)
       write(36,'(E17.8)') chisq
       close(36)
    end if
  end subroutine write_data
  
  subroutine compute_chisq(map_n,chisq,mode)
    use healpix_types
    implicit none
    integer(i4b),                                 intent(in)    :: map_n
    character(len=16),                            intent(in)    :: mode
    real(dp),                                     intent(inout) :: chisq
    real(dp)                                                    :: s, signal
    integer(i4b)                                                :: i,j,w
    
    if (trim(mode) == 'comp_sep') then
       chisq = 0.d0
       do i = 0, npix-1
          if (dang_data%masks(i,1) == missval .or. dang_data%masks(i,1) == 0.d0) cycle
          do j = 1, nbands
             s = 0.d0
             do w = 1, nfgs
                signal = dang_data%fg_map(i,map_n,j,w)
                s = s + signal
             end do
             chisq = chisq + (((dang_data%sig_map(i,map_n,j) - s)**2))/(dang_data%rms_map(i,map_n,j)**2)
          end do
       end do
       chisq = chisq/(nump*(nbands-npar))
       
    else if (trim(mode) == 'hi_fit') then
       chisq = 0.d0
       do i = 0, npix-1
          if (dang_data%masks(i,1) == missval .or. dang_data%masks(i,1) == 0.d0) cycle
          do j = 1, nbands    
             s = 0.0
             s = dang_data%gain(j)*comp%HI_amps(j)*comp%HI(i,1)*planck(par%band_nu(j)*1d9,comp%T_d(i,1))+dang_data%offset(j)
             chisq = chisq + (dang_data%sig_map(i,map_n,j)-s)**2.d0/(dang_data%rms_map(i,map_n,j)**2.d0)
          end do
       end do
       chisq = chisq/(nump*(nbands-2))
    end if
    
  end subroutine compute_chisq
  
end program dang
