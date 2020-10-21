program dang
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
    use utility_mod
    use param_mod
    use linalg_mod
    use data_mod
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
      
    integer(i4b)       :: i, j, k, l, iter, npix, nside, nmaps, ordering, m
    integer(i4b)       :: nlheader, niter, nfgs, iterations, n2fit
    integer(i4b)       :: output_iter
    real(dp)           :: nullval
    real(dp)           :: missval = -1.6375d30
    logical(lgt)       :: anynull, test, exist, output_fg
  
    character(len=128) :: template_file_01, template_file_02, mask_file, arg1
    character(len=128) :: mapfile, title, direct
  
    real(dp), allocatable, dimension(:,:,:,:)    :: fg_map
    real(dp), allocatable, dimension(:,:,:)      :: maps, rmss, model, res
    real(dp), allocatable, dimension(:,:,:)      :: nodust
    real(dp), allocatable, dimension(:,:)        :: template_01, template_02, map, rms
    real(dp), allocatable, dimension(:,:)        :: beta_s, T_d, beta_d, chi_map, mask, HI
    real(dp), allocatable, dimension(:)          :: temp01_amps, temp02_amps, temp_norm_01, temp_norm_02
    real(dp)                                     :: chisq, T_d_mean
    character(len=80), allocatable, dimension(:) :: joint_poltype
    character(len=10)                            :: solver
    character(len=5)                             :: iter_str
    real(dp), allocatable, dimension(:,:)        :: mat_l, mat_u

    ! Object Orient
    type(params)       :: par

    call init_mpi()
    call read_param_file(par)

    !----------------------------------------------------------------------------------------------------------
    ! General paramters
    template_file_01  = par%temp_file(1)
    tqu(1)            = 'T'
    tqu(2)            = 'Q'
    tqu(3)            = 'U'
    i                 = getsize_fits(template_file_01, nside=nside, ordering=ordering, nmaps=nmaps)
    npix              = nside2npix(nside) 
    nbands            = par%numband
    nfgs              = par%ncomp+par%ntemp
    nlheader          = size(header)
    nmaps             = nmaps
    niter             = par%ngibbs               ! # of MC-MC iterations
    iterations        = par%nsample              ! # of iterations in the samplers
    output_iter       = par%iter_out             ! Output maps every <- # of iterations
    output_fg         = par%output_fg            ! Option for outputting foregrounds for all bands
    direct            = par%outdir               ! Output directory name
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Array allocation
    allocate(mask(0:npix-1,1), nodust(0:npix-1,nmaps,nbands))
    allocate(temp_norm_01(3), temp_norm_02(3))
    allocate(HI(0:npix-1,nmaps))
    allocate(map(0:npix-1,nmaps))
    allocate(rms(0:npix-1,nmaps))
    allocate(joint_poltype(2))
    !----------------------------------------------------------------------------------------------------------
    beta_s     = -3.10d0    ! Synchrotron beta initial guess
    beta_d     =  1.60d0    ! Dust beta initial guess
    !----------------------------------------------------------------------------------------------------------
    ! Read maps

    do j = 1, nbands
        call read_bintab(trim(par%datadir) // trim(par%dat_noisefile(j)), &
        rms,npix,nmaps,nullval,anynull,header=header)
        rmss(:,:,j) = rms
        call read_bintab(trim(par%datadir) // trim(par%dat_mapfile(j)), &
        map,npix,nmaps,nullval,anynull,header=header)
        maps(:,:,j) = map
    end do

    call convert_maps
    
    deallocate(map,rms)

    write(*,*) 'Reading TEMPLATE01'
    call read_bintab(template_file_01,template_01,npix,nmaps,nullval,anynull,header=header)
    write(*,*) 'Reading MASKFILE'
    call read_bintab(par%mask_file,mask,npix,1,nullval,anynull,header=header)
!    call read_bintab(template_file_02,template_02,npix,nmaps,nullval,anynull,header=header)

    !----------------------------------------------------------------------------------------------------------
    ! Normalize template to avoid large values in the matrix equation
    do k = 1, nmaps
       temp_norm_01(k)  = maxval(template_01(:,k))
       template_01(:,k) = template_01(:,k)/temp_norm_01(k)
    end do

    do j = 1, nbands
       do k = 1, nmaps
          do i = 0, npix-1
             if (maps(i,k,j) == missval) maps(i,k,j) = 0.d0
          end do
       end do
    end do


    !----------------------------------------------------------------------------------------------------------
    ! Joint Sampler Info

    allocate(joint_comps(2))
    solver = par%solver
    joint_comps(1) = 'synch'
    joint_comps(2) = 'template01'

    joint_poltype(1) = 'Q'
    joint_poltype(2) = 'U'
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Calculation portion
    !----------------------------------------------------------------------------------------------------------

    dang_data%fg_map    = 0.0
    dang_data%temp_amps = 0.0

    do iter = 1, niter
       ! do k = par%pol_type(1), par%pol_type(size(par%pol_type))
        
       ! -------------------------------------------------------------------------------------------------------------------
       if (par%joint_sample) then
          call sample_joint_amp(par,dang_data,comp,2,trim(solver),joint_poltype)
          ! Extrapolating A_synch to bands
          do k = par%pol_type(1), par%pol_type(size(par%pol_type))
             if (ANY(comp%joint=='synch')) then
                do i = 0, npix-1
                   do j = 1, nbands
                      dang_data%fg_map(i,k,j,1) = dang_data%fg_map(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,comp,1,par%dat_nu(j),i,k)
                   end do
                end if
                
                ! Extrapolating A_dust to bands
                if (ANY(joint_comps=='dust')) then
                   do i = 0, npix-1
                      do j = 1, nbands
                         ! if (mask(i,k) == 0.d0) then
                            ! fg_map(i,k,j,3) = missval
                         fg_map(i,k,j,3) = fg_map(i,k,par%fg_ref_loc(2),3)*compute_spectrum(par,2,par%dat_nu(j),i,k)
                      end do
                   end do
                end if
                  
                ! Applying dust templates to make dust maps
                if (ANY(joint_comps=='template01')) then
                   do i = 0, npix-1
                      do j = 1, nbands
                         !if (mask(i,k) == 0.d0 .or. mask(i,k) == missval) then
                         !   fg_map(i,k,j,2) = missval
                         !else
                         fg_map(i,k,j,2) = temp01_amps(j)*template_01(i,k)
                         !end if
                      end do
                   end do
                end if
             end if
             ! -------------------------------------------------------------------------------------------------------------------
             ! -------------------------------------------------------------------------------------------------------------------
             if (par%fg_samp_amp(1)) then
                fg_map(:,k,par%fg_ref_loc(1),1) =  sample_spec_amp(par,maps,rmss,1,k)
                do i = 0, npix-1
                   do j = 1, nbands
                      fg_map(i,k,j,1) = fg_map(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,1,par%dat_nu(j),i,k)
                   end do
                end do
             end if
             
             ! Applying dust templates to make dust maps
             if (ANY(comp%joint=='template01')) then
                do i = 0, npix-1
                   do j = 1, nbands
                      !if (mask(i,k) == 0.d0  .or. mask(i,k) == missval) then
                      !   fg_map(i,k,j,1) = missval
                      !else
                      fg_map(i,k,j,1) = fg_map(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,1,par%dat_nu(j),i,k)
                      !end if
                   end do
                end do
             end if
          end do
       end if
       ! -------------------------------------------------------------------------------------------------------------------
       ! -------------------------------------------------------------------------------------------------------------------
       if (par%fg_samp_amp(1)) then
          dang_data%fg_map(:,k,par%fg_ref_loc(1),1) =  sample_spec_amp(par,comp,dang_data%sig_map,dang_data%rms_map,1,k)
          do i = 0, npix-1
             do j = 1, nbands
                dang_data%fg_map(i,k,j,1) = dang_data%fg_map(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,comp,1,par%dat_nu(j),i,k)
             end do
          end do
       end if
       ! -------------------------------------------------------------------------------------------------------------------
       do k = par%pol_type(1), par%pol_type(size(par%pol_type))
          call compute_chisq(k,chisq,par%mode)

          if (rank == master) then
             if (mod(iter, 1) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, E10.3, a, f7.3, a, a, 6e10.3)')&
                     iter, " - chisq: " , chisq, " - A_s: ", dang_data%fg_map(100,k,par%fg_ref_loc(1),1),& 
                     " Pol_type = " // trim(tqu(k)), ' - A_d: ', dang_data%temp_amps(:,k,1)/temp_norm_01(k)
                write(*,fmt='(a)') '---------------------------------------------'
             end if
          end do
          
       else if (trim(par%mode) == 'HI_fit') then
          write(*,*) 'Right on'
          
          call read_bintab(par%hi_file,HI,npix,1,nullval,anynull,header=header)
          T_d = par%HI_Td_init
          mask  = 1.d0
          n2fit = 0
            
          do i = 0, npix-1
             if (HI(i,k) > par%thresh) then
                mask(i,k) = 0.d0
             else
                n2fit = n2fit + 1
             end if
          end do

          temp01_amps = 0.d0
          
          do j =1, nbands
             temp01_amps(j) = temp_fit(maps(:,k,j),HI(:,1),rmss(:,k,j),mask(:,1),par%dat_nu(j))
          end do

          write(*,*) 'Intial template amplitudes:'
          write(*,fmt='(6(E12.4))') temp01_amps
          
          do iter = 1, niter
             call sample_HI_T(par,k)
             
             T_d_mean = 0.d0
             do i = 0, npix-1
                if (HI(i,1) > par%thresh) then
                   cycle
                else
                   T_d_mean = T_d_mean + T_d(i,1)
                end if
             end do
             T_d_mean = T_d_mean/n2fit
             
             write(*,*) "Mean T_d: ", T_d_mean

             do j =1, nbands
                temp01_amps(j) = temp_fit(maps(:,k,j),HI(:,1),rmss(:,k,j),mask(:,1),par%dat_nu(j))
             end do

             do i = 0, npix-1
                do j = 1, nbands
                   res(i,k,j) = maps(i,k,j) - temp01_amps(j)*HI(i,k)*planck(par%dat_nu(j)*1.d9,T_d(i,1))
                end do
             end do
             
             call compute_chisq(k,chisq,par%mode)
             
             if (mod(iter, 1) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, E10.3, a, f7.3, a, 6e12.4)')&
                     iter, " - chisq: " , chisq, " - mean T_d: ", T_d_mean, ' - Amps: ', temp01_amps
             end if
             
             if (mod(iter,output_iter) .EQ. 0) then
                call write_maps(k,par%mode)
             end if
             
             call write_data(par%mode)
             
          end do
       else   
          write(*,*) 'Unrecognized SOLVER_MODE (comp_sep, HI_fit)'
          stop
       end if
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
      
      if (trim(mode) == 'comp_sep') then
         
         write(iter_str, '(i0.5)') iter
         if (output_fg .eqv. .true.) then
            do j = 1, nbands
               title = trim(direct) // trim(par%dat_label(j)) // '_dust_fit_'// trim(tqu(nm)) & 
                    // '_' // trim(iter_str) // '.fits'
               map(:,1)   = dang_data%temp_amps(j,nm,1)*dang_data%temps(:,nm,1)
               call write_bintab(map,npix,1, header, nlheader, trim(title))
               title = trim(direct) // trim(par%dat_label(j)) // '_synch_amplitude_' //  trim(tqu(nm)) &
                    // '_' // trim(iter_str) // '.fits'
               map(:,1)   = dang_data%fg_map(:,nm,j,1)
               call write_bintab(map,npix,1, header, nlheader, trim(title))
            end do
         else 
            title = trim(direct) // trim(par%dat_label(par%fg_ref_loc(1))) // '_synch_amplitude_' //  trim(tqu(nm)) &
                 // '_' // trim(iter_str) // '.fits'
            map(:,1)   = dang_data%fg_map(:,nm,par%fg_ref_loc(1),1)
            call write_bintab(map,npix,1, header, nlheader, trim(title))
         end if
         do j = 1, nbands
            title = trim(direct) // trim(par%dat_label(j)) // '_residual_' // trim(tqu(nm)) & 
                 // '_' // trim(iter_str) // '.fits'
            map(:,1)   = res(:,nm,j)
            call write_bintab(map,npix,1, header, nlheader, trim(title))
         end do
         title = trim(direct) // 'synch_beta_' // trim(tqu(nm)) // '_' // trim(iter_str) // '.fits'
         map(:,1)   = comp%beta_s(:,nm)
         call write_bintab(map,npix,1, header, nlheader, trim(title))
         dang_data%chi_map = 0.d0
         do i = 0, npix-1
            do j = 1, nbands
               dang_data%chi_map(i,nm) = dang_data%chi_map(i,nm) + (dang_data%sig_map(i,nm,j) - dang_data%fg_map(i,nm,j,1) - dang_data%fg_map(i,nm,j,2))**2.d0/dang_data%rms_map(i,nm,j)**2.d0
            end do
         end do
         dang_data%chi_map(:,nm) = dang_data%chi_map(:,nm)/(nbands+nfgs)
         title = trim(direct) // 'chisq_' // trim(tqu(nm)) // '_' // trim(iter_str) // '.fits'
         map(:,1)   = dang_data%chi_map(:,nm)
         call write_bintab(map,npix,1, header, nlheader, trim(title))
         
      else if (trim(mode) == 'HI_fit') then
         write(iter_str, '(i0.5)') iter
         do j = 1, nbands
            title = trim(direct) // trim(par%dat_label(j)) // '_residual_' // trim(iter_str) // '.fits'
            do i=0,npix-1
               if (HI(i,1) > par%thresh) then
                  map(i,1) = missval
               else
                  map(i,1) = res(i,nm,j)
               end if
            end do
            
            call write_bintab(map,npix,1,header,nlheader,trim(title))
            
            title = trim(direct) //  trim(par%dat_label(j)) // '_model_' // trim(iter_str) // '.fits'
            do i = 0, npix-1
               if (HI(i,1) > par%thresh) then
                  map(i,1) = missval
               else
                  map(i,1) = dang_data%temp_amps(j,nm,1)*HI(i,1)*planck(par%dat_nu(j)*1.d9,comp%T_d(i,1))
               end if
            end do
            call write_bintab(map,npix,1,header,nlheader,trim(title))
         end do
         title = trim(direct) // 'T_d_' // trim(iter_str) // '.fits'
         call write_bintab(comp%T_d,npix,1,header,nlheader,trim(title))
         dang_data%chi_map = 0.d0
         do i = 0, npix-1
            !if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) then
            !   chi_map(i,nm) = missval
            !   cycle
            !else
               do j = 1, nbands
                  if (HI(i,1) > par%thresh) then
                     dang_data%chi_map(i,nm) = missval
                  else
                     dang_data%chi_map(i,nm) = dang_data%chi_map(i,nm) + &
                          (dang_data%sig_map(i,nm,j) - &
                          dang_data%temp_amps(j,nm,1)*HI(i,1)*planck(par%dat_nu(j)*1.d9,comp%T_d(i,1)))**2.d0/dang_data%rms_map(i,nm,j)**2.d0
                  end if
               end do
            !end if
         end do
         title = trim(direct) // 'chisq_' // trim(iter_str) // '.fits'
         call write_bintab(map,npix,1,header,nlheader,trim(title))
      end if

    end subroutine write_maps

    subroutine write_data(mode)
        implicit none
        character(len=16), intent(in) :: mode

        if (trim(mode) == 'comp_sep') then

            title = trim(direct) // 'pixel_100_A_d_' // trim(tqu(k)) // '.dat'
            inquire(file=title,exist=exist)
            if (exist) then
                open(30,file=title, status="old",position="append", action="write")
            else
                open(30,file=title, status="new", action="write")
                write(30,*) 
            endif
            write(30,*) dang_data%fg_map(100,k,par%fg_ref_loc(1),2)
            close(30)

            title = trim(direct) // 'pixel_100_A_s_' // trim(tqu(k)) // '.dat'
            inquire(file=title,exist=exist)
            if (exist) then
                open(31,file=title, status="old",position="append", action="write")
            else
                open(31,file=title, status="new", action="write")
            endif
            write(31,*) dang_data%fg_map(100,k,par%fg_ref_loc(1),1)
            close(31)

            title = trim(direct) // 'pixel_100_beta_s_' // trim(tqu(k)) // '.dat'
            inquire(file=title,exist=exist)
            if (exist) then
                open(32,file=title, status="old",position="append", action="write")
            else
                open(32,file=title, status="new", action="write")
            endif
            write(32,*) comp%beta_s(100,k)
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

            inquire(file=trim(direct) // 'dust_' // trim(tqu(k)) // '_amplitudes.dat',exist=exist)
            if (exist) then
                open(34,file = trim(direct) // 'dust_' // trim(tqu(k)) // '_amplitudes.dat', status="old", &
                            position="append", action="write")
            else
                open(34,file = trim(direct) // 'dust_' // trim(tqu(k)) // '_amplitudes.dat', status="new", action="write")
            endif
            write(34,'(10(E17.8))') dang_data%temp_amps(:,k,1)/temp_norm_01(k)
            close(34)

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
      integer(i4b)                                                :: i,j
  
      if (trim(mode) == 'comp_sep') then

        chisq = 0.d0
        do i = 0, npix-1
            do j = 1, nbands
                s = 0.d0
                signal = dang_data%fg_map(i,map_n,j,1) + dang_data%fg_map(i,map_n,j,2)
                s = s + signal
                chisq = chisq + (((dang_data%sig_map(i,map_n,j) - s)**2))/(dang_data%rms_map(i,map_n,j)**2)*mask(i,1)
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
      end do
      
    end subroutine


  end program dang
