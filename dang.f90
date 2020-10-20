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
      
    integer(i4b)       :: i, j, k, l, ordering, m
    integer(i4b)       :: nlheader, iterations, n2fit
    integer(i4b)       :: output_iter
    real(dp)           :: nullval
    real(dp)           :: missval = -1.6375d30
    logical(lgt)       :: anynull, test, output_fg
  
    character(len=128) :: template_file_01, template_file_02, mask_file, arg1
    character(len=128) :: mapfile, title, direct
  
    real(dp), allocatable, dimension(:,:,:,:)    :: fg_map
    real(dp), allocatable, dimension(:,:,:)      :: maps, rmss, model, res
    real(dp), allocatable, dimension(:,:,:)      :: nodust
    real(dp), allocatable, dimension(:,:)        :: template_01, template_02, map, rms
    real(dp), allocatable, dimension(:,:)        :: beta_s, T_d, beta_d, chi_map, mask, HI
    real(dp), allocatable, dimension(:,:)        :: temp01_amps, temp02_amps
    real(dp), allocatable, dimension(:)          :: temp_norm_01, temp_norm_02
    real(dp)                                     :: chisq, T_d_mean
    character(len=80), dimension(180)            :: header
    character(len=80), allocatable, dimension(:) :: joint_comps
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

    !----------------------------------------------------------------------------------------------------------
    ! General paramters
    template_file_01  = par%temp_file(1)
    tqu(1)            = 'T'
    tqu(2)            = 'Q'
    tqu(3)            = 'U'
    i                 = getsize_fits(template_file_01, nside=nside, ordering=ordering, nmaps=nmaps)
    dang_data%npix    = nside2npix(nside) 
    npix              = dang_data%npix
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
    allocate(template_01(0:npix-1,nmaps), template_02(0:npix-1,nmaps))
    !allocate(maps(0:npix-1,nmaps,nbands), rmss(0:npix-1,nmaps,nbands), nodust(0:npix-1,nmaps,nbands))
    allocate(mask(0:npix-1,1))!, res(0:npix-1,nmaps,nbands), chi_map(0:npix-1,nmaps))
    !allocate(fg_map(0:npix-1,nmaps,nbands,nfgs))
    allocate(T_d(0:npix-1,nmaps), beta_d(0:npix-1,nmaps), beta_s(0:npix-1,nmaps))
    allocate(temp01_amps(nbands,nmaps), temp02_amps(nbands,nmaps), temp_norm_01(3), temp_norm_02(3))
    allocate(HI(0:npix-1,nmaps))
    allocate(map(0:npix-1,nmaps))
    allocate(rms(0:npix-1,nmaps))
    !----------------------------------------------------------------------------------------------------------
    beta_s     = -3.10d0    ! Synchrotron beta initial guess
    beta_d     =  1.60d0    ! Dust beta initial guess
    !----------------------------------------------------------------------------------------------------------
    ! Read maps

    call init_fg_map(dang_data,npix,nmaps,nbands,nfgs)
    call init_data_maps(dang_data,npix,nmaps,nbands)

    do j = 1, nbands
        call read_bintab(trim(par%datadir) // trim(par%dat_noisefile(j)), &
        rms,npix,nmaps,nullval,anynull,header=header)
        dang_data%rms_map(:,:,j) = rms
        call read_bintab(trim(par%datadir) // trim(par%dat_mapfile(j)), &
        map,npix,nmaps,nullval,anynull,header=header)
        dang_data%sig_map(:,:,j) = map
    end do

    call convert_maps(par)
    
    deallocate(map,rms)

    write(*,*) 'Reading TEMPLATE01'
    call read_bintab(template_file_01,template_01,npix,nmaps,nullval,anynull,header=header)
    write(*,*) 'Reading MASKFILE'
    call read_bintab(par%mask_file,mask,npix,1,nullval,anynull,header=header)
    ! call read_bintab(template_file_02,template_02,npix,nmaps,nullval,anynull,header=header)
    write(*,*) ''

    !----------------------------------------------------------------------------------------------------------
    ! Normalize template to avoid large values in the matrix equation
    do k = 1, nmaps
       temp_norm_01(k)  = maxval(template_01(:,k))
       template_01(:,k) = template_01(:,k)/temp_norm_01(k)
    end do

    do j = 1, nbands
       do k = 1, nmaps
          do i = 0, npix-1
             if (dang_data%sig_map(i,k,j) == missval) dang_data%sig_map(i,k,j) = 0.d0
          end do
       end do
    end do

    !----------------------------------------------------------------------------------------------------------
    ! Joint Sampler Info

    !allocate(joint_comps(2))
    solver = par%solver
    comp%joint(1) = 'synch'
    comp%joint(2) = 'template01'

    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Calculation portion
    !----------------------------------------------------------------------------------------------------------

    dang_data%fg_map = 0.0
    temp01_amps = 0.0
    do iter = 1, niter
       do k = par%pol_type(1), par%pol_type(size(par%pol_type))
        
          ! -------------------------------------------------------------------------------------------------------------------
          if (par%joint_sample) then
             call sample_joint_amp(par,dang_data,comp,k,trim(solver))
             ! Extrapolating A_synch to bands
             if (ANY(joint_comps=='synch')) then
                do i = 0, npix-1
                   do j = 1, nbands
                      dang_data%fg_map(i,k,j,1) = dang_data%fg_map(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,comp,1,par%dat_nu(j),i,k)
                   end do
                end do
             end if
             
             ! Extrapolating A_dust to bands
             if (ANY(joint_comps=='dust')) then
                do i = 0, npix-1
                   do j = 1, nbands
                      dang_data%fg_map(i,k,j,3) = dang_data%fg_map(i,k,par%fg_ref_loc(2),3)*compute_spectrum(par,comp,2,par%dat_nu(j),i,k)
                   end do
                end do
             end if
                  
             ! Applying dust templates to make dust maps
             if (ANY(joint_comps=='template01')) then
                do i = 0, npix-1
                   do j = 1, nbands
                      dang_data%fg_map(i,k,j,2) = temp01_amps(j,k)*template_01(i,k)
                   end do
                end do
             end if
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
          call compute_chisq(k,chisq,par%mode)
       
          if (rank == master) then
             if (mod(iter, 1) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, E10.3, a, f7.3, a, a, 6e10.3)')&
                     iter, " - chisq: " , chisq, " - A_s: ", dang_data%fg_map(100,k,par%fg_ref_loc(1),1),& 
                     " Pol_type = " // trim(tqu(k)), ' - A_d: ', temp01_amps(:,k)/temp_norm_01(k)
                write(*,fmt='(a)') '---------------------------------------------'
             end if
          end if
       end do
       ! -------------------------------------------------------------------------------------------------------------------
       ! Jointly sample synchrotron beta
       ! -------------------------------------------------------------------------------------------------------------------
       if (par%fg_samp_inc(1,1)) then
          nodust(:,:,:) = dang_data%sig_map(:,:,:)-dang_data%fg_map(:,:,:,2)
          call sample_index(par,comp,nodust,par%fg_samp_nside(1,1),1,-1)
          do i = 0, npix-1
             do j = 1, nbands
                dang_data%fg_map(i,k,j,1) = dang_data%fg_map(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,comp,1,par%dat_nu(j),i,k)
             end do
          end do
       end if
       ! -------------------------------------------------------------------------------------------------------------------
             
       res = dang_data%sig_map
       do k = par%pol_type(1), par%pol_type(size(par%pol_type))
          do j = 1, nfgs
             res(:,k,:)  = res(:,k,:) - dang_data%fg_map(:,k,:,j)
          end do
       
          call compute_chisq(k,chisq,par%mode)
       
          if (rank == master) then
             if (mod(iter, 1) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, E10.3, a, f7.3, a, f8.4, a, 6e10.3)')&
                     iter, " - chisq: " , chisq, " - A_s: ",&
                     dang_data%fg_map(100,k,par%fg_ref_loc(1),1),  " - beta_s: ",&
                     sum(beta_s(:,k))/npix, ' - A_d: ', temp01_amps(:,k)/temp_norm_01(k)
                write(*,fmt='(a)') '---------------------------------------------'
             end if
             if (mod(iter,output_iter) .EQ. 0) then
                call write_maps(k,par%mode)
             end if
             call write_data(par%mode)
          end if
       end do
    end do
    call mpi_finalize(ierr)
  
  contains

    !----------------------------------------------------------------------------------------------------------
    ! Functions and subroutines
    !----------------------------------------------------------------------------------------------------------
   
  

    ! This architecture of this function has not been verified yet
    subroutine sample_HI_T(self, map_n) !(self, band, nside2, map_n)
        implicit none
  
        class(params)                                          :: self
        integer(i4b), intent(in)                               :: map_n!, nside2
        integer(i4b)                                           :: nside1, npix2, nside2
        ! real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: band
        real(dp), dimension(0:npix-1,nmaps,nbands)             :: cov
        real(dp), dimension(0:npix-1,nmaps)                    :: te
        real(dp), dimension(0:npix-1)                          :: te_sample
        real(dp), allocatable, dimension(:,:,:)                :: maps_low, cov_low
        real(dp), allocatable, dimension(:,:)                  :: T_low
        real(dp), allocatable, dimension(:)                    :: sample_T_low
        real(dp), dimension(nbands)                            :: signal, tmp
        real(dp), dimension(2)                                 :: x
        real(dp)                                               :: a, b, c, num, sam, t, p, sol, naccept_t_d

        te      = T_d
        cov     = dang_data%rms_map*dang_data%rms_map

        nside2 = 8

        nside1 = npix2nside(npix)
        if (nside1 == nside2) then
            npix2 = npix
        else
            npix2 = nside2npix(nside2)
        end if
        allocate(maps_low(0:npix2-1,nmaps,nbands))
        allocate(T_low(0:npix2-1,nmaps),cov_low(0:npix2-1,nmaps,nbands))
        allocate(sample_T_low(0:npix2-1))

        ! if (nside1 /= nside2) then 
        !     if (ordering == 1) then
        !         call udgrade_ring(te,nside1,T_low,nside2)
        !     else
        !         call udgrade_nest(te,nside1,T_low,nside2)
        !     end if
        !     do j = 1, nbands
        !         if (ordering == 1) then
        !             call udgrade_ring(map2fit(:,:,j),nside1,band_low(:,:,j),nside2)
        !             call convert_nest2ring(nside1,map2fit(:,:,j))
        !             call udgrade_ring(cov(:,:,j),nside1,cov_low(:,:,j),nside2)
        !             call convert_nest2ring(nside1,dang_data%rms_map(:,:,j))
        !         else
        !             call udgrade_nest(map2fit(:,:,j),nside1,band_low(:,:,j),nside2)
        !             call udgrade_nest(dang_data%rms_map(:,:,j),nside1,cov_low(:,:,j),nside2)
        !         end if
        !     end do
        !     cov_low = sqrt(cov_low / (npix/npix2))
        ! else
        do j = 1, nbands
            maps_low(:,:,j)   = dang_data%sig_map(:,:,j)
            cov_low(:,:,j)    = cov(:,:,j)
        end do
        T_low = te
        ! end if

        x(1) = 1.d0
        do i = 0, npix2-1
            if (mask(i,1) == 0.d0  .or. mask(i,1) == missval) then
                sample_T_low(i) = missval
                cycle
            else
                a   = 0.d0
                sol = T_low(i,map_n)

                ! Chi-square from the most recent Gibbs chain update
                do j = 1, nbands
                    a = a + (((temp01_amps(j,map_n)* HI(i,1)*planck(par%dat_nu(j)*1.d9,sol)) &
                        - maps_low(i,map_n,j))**2.d0)/cov_low(i,map_n,j)
                end do
                c   = a

                do l = 1, iterations
                    ! Begin sampling from the prior
                    ! t = rand_normal(sol,par%HI_Td_std)
                    t = rand_normal(par%HI_Td_mean,par%HI_Td_std)
                    b = 0.d0
                    do j = 1, nbands
                        tmp(j) = temp01_amps(j,map_n)*HI(i,1)*planck(par%dat_nu(j)*1.d9,t)
                        b      = b + ((tmp(j)-maps_low(i,map_n,j))**2.d0)/cov_low(i,map_n,j)
                    end do
                    b = b

                    if (b < c .and. t .lt. 35.d0 .and. t .gt. 10.d0) then
                        sam = t
                        c   = b
                    else
                        x(2) = exp(0.5d0*(c-b))
                        p = minval(x)
                        call RANDOM_NUMBER(num)
                        if (num < p) then
                            if (t .lt. 35.d0 .and. t .gt. 10.d0) then
                                sam = t
                                c   = b
                            end if
                        end if
                    end if
                end do
                sol             = sam
                sample_T_low(i) = sol
            end if
        end do
        if (nside1 /= nside2) then
            if (ordering == 1) then
                call udgrade_ring(sample_T_low, nside2, te_sample, nside1)
                call convert_nest2ring(nside2, sample_T_low)
            else
                call udgrade_nest(sample_T_low, nside2, te_sample, nside1)
            end if
        else
            te_sample =  sample_T_low
        end if
        T_d(:,1) = te_sample

        deallocate(maps_low)
        deallocate(T_low)
        deallocate(cov_low)
        deallocate(sample_T_low)

    end subroutine sample_HI_T
    ! ------------------------------------------------------------

    subroutine sample_joint_amp(para, dat, compo, map_n, method)
        !------------------------------------------------------------------------
        ! Solving the matrix equation Ab = c                                    |
        !                                                                       |
        ! For this purpose, we solve the equation:                              |
        !         sum_nu ((T_nu)^T N^-1 T_nu)amp = sum_nu ((T_nu)^T N^-1 d_nu)  |
        !                                                                       |
        ! where T_nu and a are defined below, and d_nu is the data at band nu.  |
        !                                                                       |
        ! For the matrix equation, A = sum_nu ((T_nu)^T N_nu^-1 T_nu), b = amp, |
        !                      and c =  sum_nu ((T_nu)^T d_nu).                 !                
        ! If I have this correct, A should be n x n where n=npix+nband-1, which |
        ! checks out since a is npix+nband-1.                                   |
        !------------------------------------------------------------------------

        implicit none
        type(params)                              :: para
        type(data),                intent(inout)  :: dat
        type(component)                           :: compo
        integer(i4b),              intent(in)     :: map_n
        character(len=*),          intent(in)     :: method
        real(dp), allocatable, dimension(:,:)     :: A, val
        integer(i4b), allocatable, dimension(:,:) :: col_ptr, row_ind
        real(dp), allocatable, dimension(:,:)     :: mat_l, mat_u, unc_a_s
        real(dp), allocatable, dimension(:)       :: b, c, d, rand, samp, unc_a_d
        character(len=256)                        :: title
        integer(i4b)                              :: x, y, z, nfit1, nfit2, w, l, m, n
        integer(i4b)                              :: vi, ci, ri, co, nnz, nnz_a
        integer(i4b)                              :: info

        real(dp)                                  :: q, t6, t7

        if (rank == master) then
           write(*,fmt='(a)') 'Starting joint sampling for synch and dust_template.'
           t1 = mpi_wtime()
        end if

        ! Load which components to jointly fit for arary allocation
        ! vv These will not change based off of components
        x = dat%npix
        y = 0
        z = nbands

        do i = 1, size(compo%joint)
            if (compo%joint(i) == 'synch') then
                y = y + dat%npix
            else if (compo%joint(i) == 'dust') then
                y = y + dat%npix
            else if (compo%joint(i) == 'template01') then
                ! Must have at most nbands -1 templates to avoid degeneracy
                ! Count how many bands are not being fit
                nfit1 = 0
                do j = 1, nbands
                     if (para%temp_corr(1,j)) then
                        nfit1 = nfit1 + 1
                    end if
                end do
                y = y + nfit1
            else if (compo%joint(i) == 'template02') then
                ! Must have at most nbands -1 templates to avoid degeneracy
                ! Count how many bands are not being fit
                nfit2 = 0
                do j = 1, nbands
                    if (para%temp_corr(2,j)) then
                        nfit2 = nfit2 + 1
                    end if
                end do
                y = y + nfit2
            end if
        end do

        !------------------------------------------------------------------------------|
        ! Since A is sparse (mostly 0s), we'll try to save time and memory by          |
        ! putting all of the values, with 'pointers' into a smaller array              |
        !                                                                              |
        ! We will use the CSC (Compressed Sparse Column) scheme to save time here.     |
        ! CSC is used instead of CSR due to the nature of the sparse matrix calc.      |
        !                                                                              |
        ! There are three arrays that will be used here, namely col_ptr, row_ind, and  |
        ! val. col_ptr is size y+1 (int), row_ind is nnz (num non-zero entries) (int), |
        ! and val is also size nnz (dp).                                               |
        !                                                                              |
        ! so for each band, we have nnz = 2*npix                                       |
        !------------------------------------------------------------------------------|

        nnz = x + x 

        allocate(A(y,y))
        allocate(b(y),c(y),d(y))
        allocate(mat_l(y,y),mat_u(y,y))
        allocate(rand(y),samp(y))
        allocate(col_ptr(y+1,z),row_ind(nnz,z),val(nnz,z))

        ! write(*,*) 'Initialize'
        ! Initialize arrays
        A(:,:)            = 0.d0
        b(:)              = 0.d0
        c(:)              = 0.d0
        d(:)              = 0.d0
        rand(:)           = 0.d0
        samp(:)           = 0.d0
        mat_l(:,:)        = 0.d0

        ! write(*,*) 'Fill Template Matrix'
        ! Fill template matrix
 
        l  = 1
        do j = 1, z
           vi = 1
           ci = 1
           ri = 1
           co = 1
           col_ptr(ci,j) = 1
           ci = ci + 1
!           if (mask(i,1) == 0.d0  .or. mask(i,1) == missval) cycle
           do m = 1,x 
              val(vi,j)     = compute_spectrum(para,compo,1,para%dat_nu(j),m-1,map_n)/dat%rms_map(m-1,map_n,j)
              row_ind(ri,j) = m
              co            = co + 1
              vi            = vi + 1
              ri            = ri + 1
              col_ptr(ci,j) = co
              ci            = ci + 1
           end do
           if (para%temp_corr(1,j)) then
              do i = 1, x
                 val(vi,j)     = template_01(i-1,map_n)/dat%rms_map(i-1,map_n,j)
                 co            = co + 1
                 row_ind(ri,j) = i
                 vi            = vi + 1
                 ri            = ri + 1
              end do
              col_ptr(x:ci+l-1,j)  = col_ptr(x,j)
              col_ptr(ci+l-1,j)    = co
              col_ptr(ci+l-1:,j)   = col_ptr(ci+l-1,j)
              l = l + 1
              ! write(*,*) 'its fine, true', j, l, ci+l-1
           else
              if (ci+l-1 > y+1) then
                 ! write(*,*) 'out of bounds'
                 cycle
              else
                 col_ptr(m:ci+l-1,j)  = col_ptr(m,j)
                 col_ptr(ci+l-1,j)    = co
                 col_ptr(ci+l-1:,j)   = col_ptr(ci+l-1,j)
                 ! write(*,*) 'its fine, false', j, l, ci+l-1
              end if
           end if
        end do

        ! write(*,*) 'Compute RHS of matrix eqn.'
        ! Computing the LHS and RHS of the linear equation
        ! RHS
        w = 0 
        do m = 1, size(compo%joint)
            if (compo%joint(m) == 'synch') then
               do j=1, z
                  do i=1, x
!                     if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
                     c(i) = c(i) + 1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n)
                  end do
               end do
               w = w + x
            else if (compo%joint(m) == 'dust') then
               do j=1, z
                  do i=1, x
!                     if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
                     c(i) = c(i) + 1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*compute_spectrum(para,compo,2,para%dat_nu(j),i-1,map_n)
                  end do
               end do
               w = w + x
            end if
        end do
        do m = 1, size(compo%joint)
            if (compo%joint(m) == 'template01') then
                l = 1
                do j = 1, z
                    if (para%temp_corr(1,j)) then
                        do i = 1, x
                           ! if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
                           c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                                    template_01(i-1,map_n)
                        end do
                        l = l + 1
                    end if
                end do
                w = w + l  
            else if (compo%joint(m) == 'template02') then
                l = 1
                do j = 1, z
                    if (para%temp_corr(2,j)) then
                        do i = 1, x
                           ! if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
                           c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                                    template_02(i-1,map_n)
                        end do
                        l = l + 1
                    end if
                end do
                w = w + l
            end if
        end do

        ! if (rank == master) write(*,*) 'Compute LHS of matrix eqn.'
        !LHS

        t2 = mpi_wtime()

        do j = 1, z
           ! write(*,*) j
           A(:,:) = A(:,:) + compute_ATA_CSC(val(:,j),row_ind(:,j),col_ptr(:,j))
           ! write(*,*) "finished band ", j
        end do
        nnz_a = count(A/=0)

        t3 = mpi_wtime()

        write(*,fmt='(a,E12.4,a)') 'Sparse matrix multiply: ', t3-t2, 's.'

        ! Computation
        if (trim(method) == 'cholesky') then
           mat_u(:,:)        = 0.d0
           if (rank == master) write(*,fmt='(a)') 'Joint sampling using Cholesky Decomposition.'
           call cholesky_decomp(A,mat_l)
           mat_u  = transpose(mat_l)
           call forward_sub(mat_l,d,c)
           call backward_sub(mat_u,b,d)
        else if (trim(method) == 'cg') then
           if (rank == master) write(*,*) 'Joint sampling using CG.'
           call compute_cg(A,b,c,y,nnz_a,para%cg_iter,para%cg_converge)
        else if (trim(method) == 'lu') then
           write(*,*) 'Currently deprecated -- replace with LAPACK'
           stop
           mat_u(:,:)        = 0.d0
           if (rank == master) write(*,*) 'Joint sampling using LU Decomp'
           call LUDecomp(A,mat_l,mat_u,y)
           call forward_sub(mat_l,d,c)
           call backward_sub(mat_u,b,d)
        end if
        ! Draw a sample by cholesky decompsing A^-1, and multiplying 
        ! the subsequent lower triangular by a vector of random numbers

        if (rank == master) write(*,fmt='(a,i6)') 'Draw a sample for iteration ', iter
        do i = 1, y
           rand(i) = rand_normal(0.d0,1.d0)
        end do

        t2 = mpi_wtime()
        call cholesky_decomp(A,mat_l)
        t3 = mpi_wtime()
        write(*,fmt='(a,E12.4,a)') 'Cholesky completed in ', t3-t2, 's.'
        call forward_sub(mat_l,d,rand)
       
        b = b + d

        ! Output amplitudes to the appropriate variables
        w = 0
        do m = 1, size(compo%joint)
           if (compo%joint(m) == 'synch') then
              do i = 1, x
                 dat%fg_map(i-1,map_n,para%fg_ref_loc(1),1) = b(w+i)
              end do
              w = w + x
           else if (compo%joint(m) == 'dust') then
              do i = 1, x
                 dat%fg_map(i-1,map_n,para%fg_ref_loc(2),3) = b(w+i)
              end do
              w = w + x
           end if
        end do
        do m = 1, size(compo%joint)
           if (compo%joint(m) == 'template01') then
              l = 1
              do while (l .lt. (nfit1))
                 do j= 1, z
                    if (para%temp_corr(1,j)) then
                       temp01_amps(j,map_n) = b(w+l)
                       l = l + 1
                    else
                       temp01_amps(j,map_n) = 0.d0
                    end if
                 end do
              end do
           else if (compo%joint(m) == 'template02') then
              l = 1
              do while (l .lt. (nfit2))
                 do j= 1, z
                    if (para%temp_corr(2,j)) then
                       temp02_amps(j,map_n) = b(w+l)
                       l = l + 1
                    end if
                 end do
              end do
           end if
        end do

        if (rank == master) then
           t3 = mpi_wtime()
           write(*,fmt='(a,f10.3,a)') 'Joint Sampler completed in ', t3-t1, 's.'
        end if

        if (para%output_unc .and. iter == niter) then
           allocate(unc_a_d(nfit1))
           allocate(unc_a_s(0:dat%npix-1,1))

           write(*,*) 'Dust amplitude uncertainties: '
        
           call invert_matrix_dp(A,.true.)

           write(*,*) 'Done inverting.'
           do j = 1, x
              unc_a_s(j-1,1) = sqrt(A(j,j))
           end do
           do j = 1, nfit1
              unc_a_d(j) = sqrt(A(x+j,x+j))
           end do

           inquire(file=trim(para%outdir) // 'dust_' // trim(tqu(k)) // '_uncertainties.dat',exist=exist)
           if (exist) then
              open(40,file = trim(para%outdir) // 'dust_' // trim(tqu(k)) // '_uncertainties.dat', status="old", &
                   position="append", action="write")
           else
              open(40,file = trim(para%outdir) // 'dust_' // trim(tqu(k)) // '_uncertainties.dat', status="new", action="write")
           endif
           write(40,'(6(E17.8))') unc_a_d
           close(40)

           title = trim(para%outdir) // 'a_synch_uncertainty_'// trim(tqu(k)) // '.fits'
           call write_bintab(unc_a_s,dat%npix,1, header, nlheader, trim(title))

        end if

        ! Sure to deallocate all arrays here to free up memory
        deallocate(A)
        deallocate(b)
        deallocate(c)
        deallocate(d)
        deallocate(mat_l)
        deallocate(mat_u)
        deallocate(rand)

    end subroutine sample_joint_amp


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
               map(:,1)   = temp01_amps(j,nm)*template_01(:,nm)
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
         map(:,1)   = beta_s(:,nm)
         call write_bintab(map,npix,1, header, nlheader, trim(title))
         chi_map = 0.d0
         do i = 0, npix-1
            do j = 1, nbands
               chi_map(i,nm) = chi_map(i,nm) + (dang_data%sig_map(i,nm,j) - dang_data%fg_map(i,nm,j,1) - dang_data%fg_map(i,nm,j,2))**2.d0/dang_data%rms_map(i,nm,j)**2.d0
            end do
         end do
         chi_map(:,nm) = chi_map(:,nm)/(nbands+nfgs)
         title = trim(direct) // 'chisq_' // trim(tqu(nm)) // '_' // trim(iter_str) // '.fits'
         map(:,1)   = chi_map(:,nm)
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
                  map(i,1) = temp01_amps(j,nm)*HI(i,1)*planck(par%dat_nu(j)*1.d9,T_d(i,1))
               end if
            end do
            call write_bintab(map,npix,1,header,nlheader,trim(title))
         end do
         title = trim(direct) // 'T_d_' // trim(iter_str) // '.fits'
         call write_bintab(T_d,npix,1,header,nlheader,trim(title))
         chi_map = 0.d0
         do i = 0, npix-1
            !if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) then
            !   chi_map(i,nm) = missval
            !   cycle
            !else
               do j = 1, nbands
                  if (HI(i,1) > par%thresh) then
                     chi_map(i,nm) = missval
                  else
                     chi_map(i,nm) = chi_map(i,nm) + (dang_data%sig_map(i,nm,j) - &
                          temp01_amps(j,nm)*HI(i,1)*planck(par%dat_nu(j)*1.d9,T_d(i,1)))**2.d0/dang_data%rms_map(i,nm,j)**2.d0
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
            write(32,*) beta_s(100,k)
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
            write(34,'(10(E17.8))') temp01_amps/temp_norm_01(k)
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
