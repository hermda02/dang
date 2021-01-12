
program dang
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
    use utility_mod
    use dang_param_mod
    use linalg_mod
    use dang_data_mod
    use dang_component_mod
    use sample_mod
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
      
    integer(i4b)       :: i, j, k, l, m
    integer(i4b)       :: iterations, n2fit, chain_num
    integer(i4b)       :: output_iter, bp_iter
    logical(lgt)       :: test, output_fg
  
    character(len=128) :: template_file_01, template_file_02, mask_file, arg1
    character(len=128) :: template_file_03, template_file_04
    character(len=128) :: mapfile, title, direct
  
    real(dp), allocatable, dimension(:,:,:,:)    :: fg_map
    real(dp), allocatable, dimension(:,:,:)      :: maps, rmss, model, res
    real(dp), allocatable, dimension(:,:,:)      :: synchonly
    real(dp), allocatable, dimension(:,:)        :: map, rms
    real(dp), allocatable, dimension(:,:)        :: mask, HI
    real(dp), allocatable, dimension(:)          :: temp_norm_01, temp_norm_02
    real(dp)                                     :: chisq, T_d_mean
    character(len=80), allocatable, dimension(:) :: joint_poltype
    character(len=10)                            :: solver
    character(len=5)                             :: iter_str
    real(dp), allocatable, dimension(:,:)        :: mat_l, mat_u

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
    i                 = getsize_fits(par%temp_file(1), nside=nside, ordering=ordering, nmaps=nmaps)
    dang_data%npix    = nside2npix(nside) 
    npix              = dang_data%npix
    nbands            = par%numband
    nfgs              = par%ncomp+par%ntemp
    nlheader          = size(header)
    nmaps             = nmaps
    iterations        = par%nsample              ! # of iterations in the samplers
    if (par%bp_swap) then
       niter          = (par%bp_max - par%bp_burnin)*par%num_chains
       bp_iter        = par%bp_burnin
       chain_num      = 1
    else
       niter          = par%ngibbs               ! # of MC-MC iterations
    end if
    output_iter       = par%iter_out             ! Output maps every <- # of iterations
    output_fg         = par%output_fg            ! Option for outputting foregrounds for all bands
    direct            = par%outdir               ! Output directory name
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Array allocation
    allocate(mask(0:npix-1,1), synchonly(0:npix-1,nmaps,nbands))
    allocate(HI(0:npix-1,nmaps),joint_poltype(2))
    allocate(map(0:npix-1,nmaps))
    allocate(rms(0:npix-1,nmaps))
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Joint Sampler Info
    solver = par%solver

    joint_poltype(1) = 'Q'
    joint_poltype(2) = 'U'
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Read maps

    call init_fg_map(dang_data,npix,nmaps,nbands,nfgs)
    call init_data_maps(dang_data,npix,nmaps,nbands)
    call init_mask_maps(dang_data,npix,nmaps)
    call init_template(dang_data,npix,nmaps,par%ntemp)
    call init_temp_amps(dang_data,nbands,nmaps,par%ntemp)
    call init_synch(comp,par,npix,nmaps)
    call init_dust(comp,par,npix,nmaps)

    do i = 1, par%ntemp
       call read_bintab(par%temp_file(i), dang_data%temps(:,:,i), npix,nmaps,nullval,anynull,header=header)
    end do

    write(*,*) 'Reading MASKFILE'
    call read_bintab(par%mask_file,dang_data%masks,npix,1,nullval,anynull,header=header)
    do i = 0, npix-1
       do j = 1, nmaps
          if (dang_data%masks(i,j) == 0.d0) dang_data%masks(i,j) = missval
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

    deallocate(map,rms)

    call convert_maps(par)

    do j = 1, nbands
       ! Check to see if any maps need to be dust corrected
       if (par%dust_corr(j)) then
          call dust_correct_band(dang_data,par,comp,j)
       end if
    end do
    write(*,*) ''
    !----------------------------------------------------------------------------------------------------------
    ! Calculation portion
    !----------------------------------------------------------------------------------------------------------

    dang_data%fg_map    = 0.0
    dang_data%temp_amps = 0.0

    do iter = 1, niter
       
       ! ------------ BP SWAP CHUNK ------------------------------------------------------------------
       if (par%bp_swap) then
          if (bp_iter > par%bp_max) then
             write(*,*) ''
             write(*,fmt='(a)') 'Switching to chain '// trim(par%bp_chain_list(chain_num+1))
             write(*,*) ''
             chain_num = chain_num + 1
             bp_iter   = par%bp_burnin
          end if
          call swap_bp_maps(dang_data,par,bp_iter,par%bp_chain_list(chain_num))
          write(*,*) ''
          bp_iter = bp_iter + 1
          call convert_maps_bp(par)
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
       ! -------------------------------------------------------------------------------------------------------------------

       ! -------------------------------------------------------------------------------------------------------------------
       if (par%joint_sample) then
          call sample_joint_amp(par,dang_data,comp,2,trim(solver),joint_poltype)
          do k = par%pol_type(1), par%pol_type(size(par%pol_type))

             ! Extrapolating A_synch to bands
             if (ANY(par%joint_comp=='synch')) then
                do i = 0, npix-1
                   do j = 1, nbands
                      dang_data%fg_map(i,k,j,1) = dang_data%fg_map(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,comp,1,par%dat_nu(j),i,k)
                   end do
                end do
             end if
                
             ! Extrapolating A_dust to bands
             !if (ANY(par%joint_comp=='dust')) then
             !   write(*,*) 'Extrapolating A_dust to bands'
             !   do i = 0, npix-1
             !      do j = 1, nbands
             !         dang_data%fg_map(i,k,j,3) = dang_data%fg_map(i,k,par%fg_ref_loc(2),3)*compute_spectrum(par,comp,2,par%dat_nu(j),i,k)
             !      end do
             !   end do
             !end if
                
             ! Applying dust templates to make dust maps
             if (ANY(par%joint_comp=='template01')) then
                !write(*,*) 'Extrapolating temps to bands'
                do i = 0, npix-1
                   do j = 1, nbands
                      dang_data%fg_map(i,k,j,2) = dang_data%temp_amps(j,k,1)*dang_data%temps(i,k,1)
                   end do
                end do
             end if
             if (ANY(par%joint_comp=='template02')) then
                !write(*,*) 'Extrapolating temps to bands'
                do i = 0, npix-1
                   do j = 1, nbands
                      dang_data%fg_map(i,k,j,3) = dang_data%temp_amps(j,k,2)*dang_data%temps(i,k,2)
                   end do
                end do
             end if
             if (ANY(par%joint_comp=='template03')) then
                !write(*,*) 'Extrapolating temps to bands'
                do i = 0, npix-1
                   do j = 1, nbands
                      dang_data%fg_map(i,k,j,4) = dang_data%temp_amps(j,k,3)*dang_data%temps(i,k,3)
                   end do
                end do
             end if
             if (ANY(par%joint_comp=='template04')) then
                !write(*,*) 'Extrapolating temps to bands'
                do i = 0, npix-1
                   do j = 1, nbands
                      dang_data%fg_map(i,k,j,5) = dang_data%temp_amps(j,k,4)*dang_data%temps(i,k,4)
                   end do
                end do
             end if
          end do
       end if

       ! -------------------------------------------------------------------------------------------------------------------
       ! For sampling foreground amplitudes individually
       ! -------------------------------------------------------------------------------------------------------------------
       if (par%fg_samp_amp(1)) then
          do k = par%pol_type(1), par%pol_type(size(par%pol_type))
             dang_data%fg_map(:,k,par%fg_ref_loc(1),1) =  sample_spec_amp(par,comp,dang_data%sig_map,dang_data%rms_map,1,k)
             do i = 0, npix-1
                do j = 1, nbands
                   dang_data%fg_map(i,k,j,1) = dang_data%fg_map(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,comp,1,par%dat_nu(j),i,k)
                end do
             end do
          end do
          do k = par%pol_type(1), par%pol_type(size(par%pol_type))
             call compute_chisq(k,chisq,par%mode)
             if (rank == master) then
                if (mod(iter, 1) == 0 .or. iter == 1) then
                   write(*,fmt='(i6, a, E10.3, a, f7.3, a, a, 10e10.3)')&
                        iter, " - chisq: " , chisq, " - A_s: ", dang_data%fg_map(23000,k,par%fg_ref_loc(1),1),& 
                        " Pol_type = " // trim(tqu(k)), ' - A_d ', dang_data%temp_amps(:,k,1)
                   write(*,fmt='(a)') '---------------------------------------------'
                end if
             end if
          end do
       end if
       ! -------------------------------------------------------------------------------------------------------------------
       ! Jointly sample synchrotron beta
       ! -------------------------------------------------------------------------------------------------------------------
       if (par%fg_samp_inc(1,1)) then
          write(*,*) 'Jointly sample synch beta'
          do i = 1, par%ntemp
             synchonly(:,:,:) = dang_data%sig_map(:,:,:)-dang_data%fg_map(:,:,:,i+1)
          end do
          call sample_index(par,comp,dang_data,synchonly,par%fg_samp_nside(1,1),1,-1)
          do i = 0, npix-1
             comp%beta_s(i,:) = comp%beta_s(i,:)!*dang_data%masks(i,1)
          end do
          do i = 0, npix-1
             do j = 1, nbands
                do k = par%pol_type(1), par%pol_type(size(par%pol_type))
                   dang_data%fg_map(i,k,j,1) = dang_data%fg_map(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,comp,1,par%dat_nu(j),i,k)
                end do
             end do
          end do
       end if
 
       ! -------------------------------------------------------------------------------------------------------------------
       
       dang_data%res_map = dang_data%sig_map
       do k = par%pol_type(1), par%pol_type(size(par%pol_type))
          write(*,*) k
          do j = 1, nfgs
             dang_data%res_map(:,k,:)  = dang_data%res_map(:,k,:) - dang_data%fg_map(:,k,:,j)
          end do
          
          call compute_chisq(k,chisq,par%mode)
          
          if (rank == master) then
             if (mod(iter, 1) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, E10.3, a, f7.3, a, f8.4, a, 10e10.3)')&
                     iter, " - chisq: " , chisq, " - A_s: ",&
                     dang_data%fg_map(23000,k,par%fg_ref_loc(1),1),  " - beta_s: ",&
                     mask_avg(comp%beta_s(:,k),dang_data%masks(:,1)), ' - A_d: ', dang_data%temp_amps(:,k,1)
                write(*,fmt='(a)') '---------------------------------------------'
             end if
             if (mod(iter,output_iter) .EQ. 0) then
                call write_maps(k,par%mode)
             end if
             call write_data(par%mode)
          end if
       end do
       write(*,*) ''
    end do
    call mpi_finalize(ierr)
       
contains

    !----------------------------------------------------------------------------------------------------------
    ! Functions and subroutines
    !----------------------------------------------------------------------------------------------------------


    subroutine write_maps(nm,mode)
      implicit none

      integer(i4b), intent(in)          :: nm
      character(len=16), intent(in)     :: mode
      real(dp), dimension(0:npix-1,1)   :: map
      real(dp)                          :: s, signal
      integer(i4b)                      :: n

      write(*,*) 'Output data maps'

      if (trim(mode) == 'comp_sep') then
         
         write(iter_str, '(i0.5)') iter
         if (output_fg .eqv. .true.) then
            do j = 1, nbands
               do i = 1, par%ntemp
                  title = trim(direct) // trim(par%dat_label(j)) //'_'// trim(par%temp_label(i)) //&
                       '_'// trim(tqu(nm)) // '_' // trim(iter_str) // '.fits'
                  map(:,1)   = dang_data%fg_map(:,nm,j,i+1)
                  do n = 0, npix-1
                     if (dang_data%masks(n,1) == 0.d0 .or. dang_data%masks(n,1) == missval) then
                        map(n,1) = missval
                     end if
                  end do
                  call write_bintab(map,npix,1, header, nlheader, trim(title))
               end do
               title = trim(direct) // trim(par%dat_label(j)) // '_synch_amplitude_' //  trim(tqu(nm)) &
                    // '_' // trim(iter_str) // '.fits'
               map(:,1)   = dang_data%fg_map(:,nm,j,1)
               do n = 0, npix-1
                  if (dang_data%masks(n,1) == 0.d0 .or. dang_data%masks(n,1) == missval) then
                     map(n,1) = missval
                  end if
               end do
               call write_bintab(map,npix,1, header, nlheader, trim(title))
            end do
         else 
            title = trim(direct) // trim(par%dat_label(par%fg_ref_loc(1))) // '_synch_amplitude_' //  trim(tqu(nm)) &
                 // '_' // trim(iter_str) // '.fits'
            map(:,1)   = dang_data%fg_map(:,nm,par%fg_ref_loc(1),1)
            do n = 0, npix-1
               if (dang_data%masks(n,1) == 0.d0 .or. dang_data%masks(n,1) == missval) then
                  map(n,1) = missval
               end if
            end do
            call write_bintab(map,npix,1, header, nlheader, trim(title))
         end if
         do j = 1, nbands
            title = trim(direct) // trim(par%dat_label(j)) // '_residual_' // trim(tqu(nm)) & 
                 // '_' // trim(iter_str) // '.fits'
            map(:,1)   = dang_data%res_map(:,nm,j)
            do n = 0, npix-1
               if (dang_data%masks(n,1) == 0.d0 .or. dang_data%masks(n,1) == missval) then
                  map(n,1) = missval
               end if
            end do
            call write_bintab(map,npix,1, header, nlheader, trim(title))
         end do
         title = trim(direct) // 'synch_beta_' // trim(tqu(nm)) // '_' // trim(iter_str) // '.fits'
         map(:,1)   = comp%beta_s(:,nm)
         do n = 0, npix-1
            if (dang_data%masks(n,1) == 0.d0 .or. dang_data%masks(n,1) == missval) then
               map(n,1) = missval
            end if
         end do
         call write_bintab(map,npix,1, header, nlheader, trim(title))
         dang_data%chi_map = 0.d0
         do i = 0, npix-1
            do j = 1, nbands
               s      = 0.d0
               do l = 1, nfgs
                  signal = dang_data%fg_map(i,nm,j,l)
                  s      = s + signal
               end do
               dang_data%chi_map(i,nm) = dang_data%chi_map(i,nm) + dang_data%masks(i,1)*(dang_data%sig_map(i,nm,j) - s)**2.d0/dang_data%rms_map(i,nm,j)**2.d0
            end do
         end do
         dang_data%chi_map(:,nm) = dang_data%chi_map(:,nm)/(nbands+nfgs)
         title = trim(direct) // 'chisq_' // trim(tqu(nm)) // '_' // trim(iter_str) // '.fits'
         map(:,1)   = dang_data%chi_map(:,nm)
         do n = 0, npix-1
            if (dang_data%masks(n,1) == 0.d0 .or. dang_data%masks(n,1) == missval) then
               map(n,1) = missval
            end if
         end do
         call write_bintab(map,npix,1, header, nlheader, trim(title))
      end if

    end subroutine write_maps

    subroutine write_data(mode)
        implicit none
        character(len=16), intent(in) :: mode
        character(len=2)              :: temp_n

        if (trim(mode) == 'comp_sep') then

            title = trim(direct) // 'pixel_23000_A_d_' // trim(tqu(k)) // '.dat'
            inquire(file=title,exist=exist)
            if (exist) then
                open(30,file=title, status="old",position="append", action="write")
            else
                open(30,file=title, status="new", action="write")
                write(30,*) 
            endif
            write(30,*) dang_data%fg_map(23000,k,par%fg_ref_loc(1),2)
            close(30)

            title = trim(direct) // 'pixel_23000_A_s_' // trim(tqu(k)) // '.dat'
            inquire(file=title,exist=exist)
            if (exist) then
                open(31,file=title, status="old",position="append", action="write")
            else
                open(31,file=title, status="new", action="write")
            endif
            write(31,*) dang_data%fg_map(23000,k,par%fg_ref_loc(1),1)
            close(31)

            title = trim(direct) // 'pixel_23000_beta_s_' // trim(tqu(k)) // '.dat'
            inquire(file=title,exist=exist)
            if (exist) then
                open(32,file=title, status="old",position="append", action="write")
            else
                open(32,file=title, status="new", action="write")
            endif
            write(32,*) comp%beta_s(23000,k)
            close(32)

            title = trim(direct) // 'total_chisq_' // trim(tqu(k)) // '.dat'
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
               title = trim(direct) //  trim(par%temp_label(i)) // '_' //trim(tqu(k)) // '_amplitudes.dat'
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

        else if (trim(mode) == 'HI_fit') then
            title = trim(direct)//'HI_amplitudes.dat'
            inquire(file=title,exist=exist)
            if (exist) then
                open(35,file=title,status="old",position="append",action="write") 
            else
                open(35,file=title,status="new",action="write")
            end if
            write(35,'(6(E17.8))')
            close(35)

            title = trim(direct)//'HI_chisq.dat'
            inquire(file=title,exist=exist)
            if (exist) then
                open(36,file=title,status="old",position="append",action="write") 
            else
                open(36,file=title,status="new",action="write")
            end if
            call compute_chisq(k,chisq,par%mode)
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
           if (dang_data%masks(i,1) == missval) cycle
           do j = 1, nbands
              s = 0.d0
              do w = 1, nfgs
                 signal = dang_data%fg_map(i,map_n,j,w)
                 s = s + signal
              end do
              chisq = chisq + (((dang_data%sig_map(i,map_n,j) - s)**2))/(dang_data%rms_map(i,map_n,j)**2)!*mask(i,1)
           end do
        end do
        chisq = chisq/(npix+nbands+nfgs)

    else if (trim(mode) == 'HI_fit') then
        chisq = 0.d0
        do i = 0, npix-1
            if (HI(i,1) > par%thresh) then
                cycle
            else
                do j = 1, nbands    
                    ! s = temp01_amps(j)*HI(i,1)*planck(par%dat_nu(j)*1.d9,T_d(i,1))
                    ! chisq = chisq + (dang_data%sig_map(i,map_n,j)-s)**2.d0/(dang_data%rms_map(i,map_n,j)**2.d0)*mask(i,1)
                    chisq = chisq + (res(i,map_n,j))**2.d0/(dang_data%rms_map(i,map_n,j)**2.d0)
                end do
            end if
        end do
        chisq = chisq/(n2fit+1)
    end if

    end subroutine compute_chisq

    subroutine convert_maps(self)
      implicit none
      type(params), intent(inout) :: self
      real(dp)                    :: cmb_to_rj, y
      
      do j = 1, nbands
         if (.not. par%bp_map(j)) then
            if (trim(self%dat_unit(j)) == 'uK_RJ') then
               cycle
            else if (trim(self%dat_unit(j)) == 'uK_cmb') then
               write(*,*) 'Putting band ', trim(self%dat_label(j)), ' from uK_cmb to uK_RJ.'
               y           = h*(self%dat_nu(j)*1.0d9) / (k_B*T_CMB)
               cmb_to_rj   = (y**2.d0*exp(y))/(exp(y)-1)**2.d0
               dang_data%sig_map(:,:,j) = cmb_to_rj*dang_data%sig_map(:,:,j)
            else if (trim(self%dat_unit(j)) == 'MJy/sr') then
               write(*,*) 'Unit conversion not for MJy/sr not added!'
               stop
            else
               write(*,*) 'Not a unit, dumbass!'
               stop
            end if
         end if
      end do
      
    end subroutine convert_maps

    subroutine convert_maps_bp(self)
      implicit none
      type(params), intent(inout) :: self
      real(dp)                    :: cmb_to_rj, y
      real(dp)                    :: mjy_to_rj
      
      do j = 1, nbands
         if (par%bp_map(j)) then
            if (trim(self%dat_unit(j)) == 'uK_RJ') then
               cycle
            else if (trim(self%dat_unit(j)) == 'uK_cmb') then
               write(*,*) 'Putting band ', trim(self%dat_label(j)), ' from uK_cmb to uK_RJ.'
               y           = h*(self%dat_nu(j)*1.0d9) / (k_B*T_CMB)
               cmb_to_rj   = (y**2.d0*exp(y))/(exp(y)-1)**2.d0
               dang_data%sig_map(:,:,j) = cmb_to_rj*dang_data%sig_map(:,:,j)
            else if (trim(self%dat_unit(j)) == 'MJy/sr') then
               y           = h*(self%dat_nu(j)*1.0d9) / (k_B*T_CMB)
               mjy_to_rj   = 1.0/(1e14*(2.0*h*(self%dat_nu(j)*1.0d9)**3.d0/&
                    &(c**2.d0*(exp(y)-1)))*(exp(y)/(exp(y)-1.d0))*(h*(self%dat_nu(j)*1.0d9))/(k_B*T_CMB**2.d0))
               dang_data%sig_map(:,:,j) = mjy_to_rj*dang_data%sig_map(:,:,j)
               write(*,*) 'Unit conversion not for MJy/sr not added!'
               stop
            else
               write(*,*) 'Not a unit, dumbass!'
               stop
            end if
         end if
      end do
      
    end subroutine convert_maps_bp

    function mask_avg(array,mask)
      real(dp), dimension(:), intent(in) :: array
      real(dp), dimension(:), intent(in) :: mask
      real(dp)                           :: sum, mask_avg
      integer(i4b)                       :: mask_sum

      sum = 0.d0
      mask_sum = 0

      do i = 1, npix
         if (mask(i) == missval) then
            cycle
         else
            sum      = sum + array(i)
            mask_sum = mask_sum + 1
         end if
      end do

      mask_avg = sum/mask_sum
     
    end function mask_avg
  end program dang
