program dust_fit
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
    use init_mod
    use linalg_mod
    use foreground_mod
    use data_mod
    implicit none
  
    !------------------------------------------------------------------------------------------------------
    ! Daniel Herman 2020                                                                                  |
    !                                                                                                     |  
    ! This program fits the Planck (NPIPE) 353 GHz dust map to LFI bands to set constraints on the level  |
    ! of polarized emission from Anomalous Microwave Emission in the LFI bands.                           |
    !                                                                                                     |  
    !-----------------------------------------------------------------------------------------------------|  
  
    !-----------------------------------------------------------------------------------------------------|  
    ! What we want here: we are taking a joint fit of CMB, synchrotron, and dust emission. In order to do |
    ! effectively, this will be a mulit-frequency Gibbs sampler. Using the BeyondPlanck LFI maps, along   |
    ! side the WMAP data, we iteratively fit CMB (per pixel), synchrotron (per pixel), and dust (global). |
    ! Once this has been done, estimates to the level of AME polarization will be made using the lowest   |
    ! frequency bands used here.                                                                          |
    !-----------------------------------------------------------------------------------------------------|  
      
    integer(i4b)       :: i, j, k, l, iter, npix, nside, nmaps, ordering
    integer(i4b)       :: beta_samp_nside, nlheader, niter, nfgs, iterations
    integer(i4b)       :: output_iter, like_iter, m
    real(dp)           :: nullval
    real(dp)           :: missval = -1.6375d30
    logical(lgt)       :: anynull, double_precision, test, exist, output_fg
  
    character(len=128) :: template_file_01, template_file_02, mask_file, arg1
    character(len=128) :: mapfile, title, direct
  
    real(dp), allocatable, dimension(:,:,:,:)    :: fg_amp
    real(dp), allocatable, dimension(:,:,:)      :: maps, rmss, model, res
    real(dp), allocatable, dimension(:,:,:)      :: synch_map, dust_map, cmb_map, nodust
    real(dp), allocatable, dimension(:,:)        :: template_01, template_02, map, rms
    real(dp), allocatable, dimension(:,:)        :: beta_s, T_d, beta_d, chi_map, mask, HI
    real(dp), allocatable, dimension(:)          :: nuz, chi, tump, accept, prob, par, dust_amps
    real(dp)                                     :: nu_ref_s, beta_s_std, beta_s_mu
    real(dp)                                     :: nu_ref_d, beta_d_std, beta_d_mu, T_d_mu
    real(dp)                                     :: chisq, temp_norm_01, temp_norm_02
    character(len=80), dimension(180)            :: header
    character(len=80), dimension(3)              :: tqu
    character(len=80), allocatable, dimension(:) :: joint_comps
    character(len=5)                             :: iter_str
    logical(lgt), allocatable, dimension(:)      :: j_corr01, j_corr02

    real(dp), allocatable, dimension(:,:)        :: mat_test, mat_l, mat_u
    real(dp), allocatable, dimension(:)          :: x, b, d

    ! Object Orient
    type(fg_comp)      :: frgrnd(3)
    type(band)         :: bands(6) 

    !----------------------------------------------------------------------------------------------------------
    ! General paramters
    template_file_01  = 'data/test_data/npipe6v20_353_map_Q_n0004.fits'
    template_file_02  = 'data/temp_synch_030_n0004.fits'
    mask_file         = 'data/mask_fullsky_n0004.fits'
    tqu(1)            = 'T'
    tqu(2)            = 'Q'
    tqu(3)            = 'U'
    i                 = getsize_fits(template_file_01, nside=nside, ordering=ordering, nmaps=nmaps)
    npix              = nside2npix(nside) 
    nbands            = 6
    nfgs              = 2
    nlheader          = size(header)
    nmaps             = 1

    niter             = 1       ! # of MC-MC iterations
    iterations        = 100        ! # of iterations in the samplers
    output_iter       = 1        ! Output maps every <- # of iterations
    like_iter         = 1000       ! Output likelihood test every <- # of iterations
    nu_ref_s          = 45.0d0     ! Synchrotron reference frequency
    nu_ref_d          = 353.d0     ! Dust reference frequency
    beta_s_mu         = -3.10d0    ! \beta_synch Gaussian prior mean
    beta_s_std        = 0.1d0      ! \beta_synch Gaussian prior std
    beta_samp_nside   = 4          ! \beta_synch nside sampling
    output_fg         = .true.     ! Option for outputting foregrounds for all bands
    test              = .true.     ! Testing Metropolis-Hasting apparatus
    !----------------------------------------------------------------------------------------------------------

    call getarg(1,arg1)
    direct = arg1

    !----------------------------------------------------------------------------------------------------------
    ! Array allocation
    allocate(template_01(0:npix-1,nmaps), template_02(0:npix-1,nmaps), j_corr01(nbands), j_corr02(nbands))
    allocate(maps(0:npix-1,nmaps,nbands), rmss(0:npix-1,nmaps,nbands), nodust(0:npix-1,nmaps,nbands))
    allocate(mask(0:npix-1,1), res(0:npix-1,nmaps,nbands), chi_map(0:npix-1,nmaps))
    allocate(fg_amp(0:npix-1,nmaps,nbands,nfgs), beta_s(0:npix-1,nmaps))
    allocate(T_d(0:npix-1,nmaps), beta_d(0:npix-1,nmaps), HI(0:npix-1,nmaps))
    allocate(cmb_map(0:npix-1,nmaps,nbands), dust_map(0:npix-1,nmaps,nbands), synch_map(0:npix-1,nmaps,nbands))
    allocate(dust_amps(nbands))
    allocate(map(0:npix-1,nmaps))
    allocate(rms(0:npix-1,nmaps))
    !----------------------------------------------------------------------------------------------------------
    beta_s     = -3.10d0    ! Synchrotron beta initial guess
    beta_d     = 1.60d0     ! Dust beta initial guess

    ! Load band info
    call bands(1)%get_map('norm_pol_020_n0004.fits')
    call bands(2)%get_map('norm_pol_045_n0004.fits')
    call bands(3)%get_map('norm_pol_070_n0004.fits')
    call bands(4)%get_map('norm_pol_100_n0004.fits')
    call bands(5)%get_map('norm_pol_200_n0004.fits')
    call bands(6)%get_map('norm_pol_353_n0004.fits')

    call bands(1)%get_rms('norm_pol_020_rms_n0004.fits')
    call bands(2)%get_rms('norm_pol_045_rms_n0004.fits')
    call bands(3)%get_rms('norm_pol_070_rms_n0004.fits')
    call bands(4)%get_rms('norm_pol_100_rms_n0004.fits')
    call bands(5)%get_rms('norm_pol_200_rms_n0004.fits')
    call bands(6)%get_rms('norm_pol_353_rms_n0004.fits')

    call bands(1)%get_nu(20.d0)
    call bands(2)%get_nu(45.d0)
    call bands(3)%get_nu(70.d0)
    call bands(4)%get_nu(100.d0)
    call bands(5)%get_nu(200.d0)
    call bands(6)%get_nu(353.d0)

    call bands(1)%get_label('norm_pol_020')
    call bands(2)%get_label('norm_pol_045')
    call bands(3)%get_label('norm_pol_070')
    call bands(4)%get_label('norm_pol_100')
    call bands(5)%get_label('norm_pol_200')
    call bands(6)%get_label('norm_pol_353')
    
    ! Load foreground info
    frgrnd(1)%type       = 'synch'
    frgrnd(1)%nu_ref     = 45.d0
    frgrnd(1)%gauss_mean = -3.10d0
    frgrnd(1)%gauss_std  = 0.1d0
    frgrnd(1)%uni_min    = -4.5d0
    frgrnd(1)%uni_max    = -1.5d0
    frgrnd(1)%loc        = minloc(abs(bands%nu-nu_ref_s),1)

    !----------------------------------------------------------------------------------------------------------
    ! Read maps

    do j = 1, nbands
        call read_bintab('data/test_data/norm_pol/' // trim(bands(j)%rms_map), &
        rms,npix,nmaps,nullval,anynull,header=header)
        rmss(:,:,j) = rms
        call read_bintab('data/test_data/norm_pol/' // bands(j)%sig_map, &
        map,npix,nmaps,nullval,anynull,header=header)
        maps(:,:,j) = map
    end do

    deallocate(map,rms)

    call read_bintab(template_file_01,template_01,npix,nmaps,nullval,anynull,header=header)
    call read_bintab(mask_file,mask,npix,1,nullval,anynull,header=header)
    call read_bintab(template_file_02,template_02,npix,nmaps,nullval,anynull,header=header)

    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Choose which bands for fitting templates
    j_corr01(1) = .false.
    j_corr01(2) = .true.
    j_corr01(3) = .true.
    j_corr01(4) = .true.
    j_corr01(5) = .true.
    j_corr01(6) = .true.
    j_corr02(1) = .false.
    j_corr02(2) = .true.
    j_corr02(3) = .true.
    j_corr02(4) = .true.
    j_corr02(5) = .true.
    j_corr02(6) = .true.
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Metropolis-Hastings testing things
    allocate(chi(iterations*niter+1), tump(iterations*niter+1), accept(iterations*niter+1), prob(iterations*niter+1))
    !----------------------------------------------------------------------------------------------------------
    ! Normalize template to avoid large values in the matrix equation
    temp_norm_01 = maxval(template_01)
    template_01  = template_01/temp_norm_01
    temp_norm_02 = maxval(template_02)
    template_02  = template_02/temp_norm_02
    !----------------------------------------------------------------------------------------------------------
    ! Calculation portion
    !----------------------------------------------------------------------------------------------------------

    allocate(joint_comps(2))

    do k = 1, nmaps
        
        ! if (k == 1) then
        !     write(*,*) 'Sampling Temperature'
        !     write(*,*) '-----------------------'
        if (k == 1) then
            write(*,*) 'Stokes Q'
            write(*,*) '-----------------------'
        else if (k == 2) then
            write(*,*) 'Stokes U'
            write(*,*) '-----------------------'
        end if 

        do iter = 1, niter

            joint_comps(1) = 'synch'
            joint_comps(2) = 'template01'
    
            call sample_joint_amp(npix,k,'lu')  ! Method possibilities are 'cg', 'lu', and 'cholesky'
            dust_amps = fg_amp(0,k,:,2)

            ! -------------------------------------------------------------------------------------------------------------------
            ! Extrapolating A_synch to bands
            if (ANY(joint_comps=='synch')) then
                do i = 0, npix-1
                    frgrnd(1)%p(1) = beta_s(i,k)
                    do j = 1, nbands
                        fg_amp(i,k,j,1) = fg_amp(i,k,frgrnd(1)%loc,1)*compute_spectrum(frgrnd(1),bands(j)%nu)
                    end do
                end do
                synch_map(:,k,:)        = fg_amp(:,k,:,1)
            end if

            ! Extrapolating A_dust to bands
            if (ANY(joint_comps=='dust')) then
                write(*,*) 'works, dust'
                do i = 0, npix-1
                    frgrnd(3)%p(1) = beta_d(i,k)
                    frgrnd(3)%p(2) = T_d(i,k)
                    do j = 1, nbands
                        fg_amp(i,k,j,3) = fg_amp(i,k,frgrnd(3)%loc,3)*compute_spectrum(frgrnd(3),bands(j)%nu)
                    end do
                end do
                dust_map(:,k,:)         = fg_amp(:,k,:,3)
            end if

            ! Applying dust templates to make dust maps
            if (ANY(joint_comps=='template01')) then
                do j = 1, nbands
                    dust_map(:,k,j) = fg_amp(:,k,j,2)*template_01(:,k)
                end do
            end if

            nodust = maps-dust_map
            ! -------------------------------------------------------------------------------------------------------------------
            !call sample_index(nodust,'synch',beta_samp_nside,k)
            !do i = 0, npix-1
            !   par(1) = beta_s(i,k)
            !   do j = 1, nbands
            !       fg_amp(i,k,j,1) = fg_amp(i,k,frgrnd(1)%loc,1)*compute_spectrum('synch',nuz(j),par)
            !   end do
            !end do
            !synch_map(:,k,:)        = fg_amp(:,k,:,1)
            ! -------------------------------------------------------------------------------------------------------------------

            res       = maps - synch_map - dust_map

            call compute_chisq(fg_amp,k)

            if (mod(iter, 1) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, f10.3, a, f7.3, a, f8.4, a, 6e10.3)')&
                 iter, " - chisq: " , chisq, " - A_s: ",&
                 fg_amp(100,k,frgrnd(1)%loc,1),  " - beta_s: ",&
                 sum(beta_s(:,k))/npix, ' - A_d: ', dust_amps/temp_norm_01
            end if

            call write_data
            if (mod(iter,output_iter) .EQ. 0) then
                call write_maps(k)
            end if
        end do    
    end do
  
  contains

    !----------------------------------------------------------------------------------------------------------
    ! Functions and subroutines
    !----------------------------------------------------------------------------------------------------------
 
    function rand_normal(mean,stdev) result(c)
         double precision :: mean,stdev,c,temp(2),theta,r
         if (stdev <= 0.0d0) then
            write(*,*) "Standard Deviation must be positive."
         else
            call RANDOM_NUMBER(temp)
            r=(-2.0d0*log(temp(1)))**0.5
            theta = 2.0d0*PI*temp(2)
            c= mean+stdev*r*sin(theta)
      end if
    end function

    function compute_spectrum(self, freq)! result(spectrum)
        implicit none
        class(fg_comp)                 :: self
        real(dp),           intent(in) :: freq
        real(dp)                       :: y, rj_cmb ! Conversion factor from K_{RJ} -> K_{CMB}
        real(dp)                       :: z, compute_spectrum
        real(dp)                       :: spectrum

        y = h*(freq*1.0d9) / (k_B*T_CMB)
        rj_cmb = ((exp(y)-1)**2.d0)/(y**2.d0*exp(y))

        if (trim(self%type) == 'synch') then
            compute_spectrum = (freq/self%nu_ref)**self%p(1) !*rj_cmb
        else if (trim(self%type) == 'dust') then
            z = h / (k_B*self%p(2))
            compute_spectrum = (exp(z*self%nu_ref*1d9)-1.d0) / (exp(z*freq*1d9)-1.d0) * (freq/self%nu_ref)**(self%p(1)+1.d0)!*rj_cmb
        end if

    end function compute_spectrum

    function planck(fre,T)
        implicit none
        real(dp), intent(in)  :: fre
        real(dp), intent(in)  :: T
        real(dp)              :: planck
        ! Output in units of [W sr^-1 m^-2 Hz^-1]
        planck  = ((2.d0*h*fre**3.d0)/(c**2.d0))*(1.d0/(exp((h*fre)/(k_B*T))-1))
    end function planck

    function temp_fit(data,template,noise,t_mask)
        implicit none
    
        real(dp), dimension(0:npix-1), intent(in) :: data, template, noise, t_mask
        real(dp), dimension(0:npix-1)             :: cov
        real(dp)                                  :: temp_fit, norm, sam, p
        real(dp)                                  :: amp, old, sum1, sum2, num
        real(dp)                                  :: chi, chi_0, chi_00
    
        cov = noise*2

        ! Uncertainty, used for sampling.
        norm  = sum((template(:)**2.d0)/cov(:))/sum(t_mask)
        
        sum1 = 0.d0
        sum2 = 0.d0
    
        do i=0,npix-1
            sum1   = sum1 + (data(i)*template(i))/cov(i)*t_mask(i)
            sum2   = sum2 + template(i)**2.d0/cov(i)*t_mask(i)
        end do

        ! Don't allow negative amplitudes, if negative, set to 0.d0.
        if (sum1 < 0.d0) then
            amp = 0.d0
        else
            amp   = sum1/sum2 + rand_normal(0.d0,1.d0)/sqrt(norm)
        end if
        temp_fit  = amp
  
    end function temp_fit
  
    function sample_spec_amp(data,type,noise)
        !------------------------------------------------------------------------
        ! Samples spectral amplitude (per pixel), following the spectrum of foreground 'type'. Returns a full map of amplitudes.
        !------------------------------------------------------------------------
        implicit none
  
        real(dp), dimension(0:npix-1,nbands), intent(in)       :: data, noise
        character(len=*),                     intent(in)       :: type
        real(dp)                                               :: sum1, sum2, spec
        real(dp)                                               :: chi, chi_0, chi_00, p
        real(dp)                                               :: amp, num, t, sam
        real(dp), dimension(2)                                 :: pars
        real(dp), dimension(nbands)                            :: tmp
        real(dp), dimension(0:npix-1)                          :: norm
        real(dp), dimension(0:npix-1,nbands)                   :: cov
        real(dp), dimension(0:npix-1)                          :: sample_spec_amp

        cov = noise*2

        do i = 0, npix-1
            if (trim(type) == 'synch') then 
                pars(1) = beta_s(i,k)
            else if (trim(type) == 'dust') then
                pars(1) = beta_d(i,k)
                pars(2) = T_d(i,k)
            end if
            sum1 = 0.0d0
            sum2 = 0.0d0
            do j = 1, nbands
                spec          = compute_spectrum(frgrnd(1),nuz(j))
                sum1           = sum1 + (data(i,j)*spec)/cov(i,j)
                sum2           = sum2 + (spec)**2.d0/cov(i,j)
                norm(i)        = norm(i) + ((spec)**2.d0)/cov(i,j)
            end do
            norm(i)            = norm(i)!/nbands
            amp                = sum1/sum2
            sample_spec_amp(i) = amp + rand_normal(0.d0,1.d0)/sqrt(norm(i))
        end do

    end function sample_spec_amp
  
    subroutine sample_index(data, type, nside2, map_n)
        implicit none
  
        integer(i4b),                               intent(in) :: map_n, nside2
        character(len=*),                           intent(in) :: type
        real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: data
        integer(i4b)                                           :: nside1, npix2
        real(dp), dimension(0:npix-1,nmaps,nbands)             :: map2fit, cov
        real(dp), dimension(0:npix-1,nmaps)                    :: indx
        real(dp), dimension(0:npix-1)                          :: indx_sample
        real(dp), allocatable, dimension(:,:,:)                :: data_low, fg_amp_low, rms_low
        real(dp), allocatable, dimension(:,:)                  :: indx_low
        real(dp), allocatable, dimension(:)                    :: indx_sample_low
        real(dp), dimension(nbands)                            :: signal, tmp
        real(dp), dimension(2)                                 :: x, pars, prior
        real(dp)                                               :: a, b, c, num, sam, t, p, sol
        real(dp)                                               :: mu, sigma, d1, d2

        real(dp)                                               :: naccept   
        logical                                                :: exist

        !------------------------------------------------------------------------
        ! Spectral index sampler, using the Metropolis-Hastings approach.
        !------------------------------------------------------------------------

        map2fit = data
        cov     = rmss*rmss

        !------------------------------------------------------------------------
        ! Load priors for the appropriate spectrum
        !------------------------------------------------------------------------
        if (trim(type) == 'synch') then 
            prior(1) = beta_s_mu
            prior(2) = beta_s_std
            indx     = beta_s
        else if (trim(type) == 'dust') then 
            prior(1) = beta_d_mu
            prior(2) = beta_d_std
            indx     = beta_d
        end if


        if (mod(iter,output_iter) .EQ. 0) then
           write(*,*) 'Sampling ' // trim(type) // ' beta at nside', beta_samp_nside
        end if

        !------------------------------------------------------------------------
        ! Check to see if the data nside is the same as the sampling nside
        ! If not equal, downgrade the data before sampling
        !------------------------------------------------------------------------
        nside1 = npix2nside(npix)
        if (nside1 == nside2) then
            npix2 = npix
        else
            npix2 = nside2npix(nside2)
        end if
        allocate(data_low(0:npix2-1,nmaps,nbands),fg_amp_low(0:npix2-1,nmaps,nbands))
        allocate(indx_low(0:npix2-1,nmaps),rms_low(0:npix2-1,nmaps,nbands))
        allocate(indx_sample_low(0:npix2-1))

        if (nside1 /= nside2) then 
            if (ordering == 1) then
                call udgrade_ring(indx,nside1,indx_low,nside2)
            else
                call udgrade_nest(indx,nside1,indx_low,nside2)
            end if
            do j = 1, nbands
                if (ordering == 1) then
                    call udgrade_ring(map2fit(:,:,j),nside1,data_low(:,:,j),nside2)
                    call convert_nest2ring(nside1,map2fit(:,:,j))
                    call udgrade_ring(fg_amp(:,:,j,1),nside1,fg_amp_low(:,:,j),nside2)
                    call convert_nest2ring(nside1,fg_amp(:,:,j,1))
                    call udgrade_ring(cov(:,:,j),nside1,rms_low(:,:,j),nside2)
                    call convert_nest2ring(nside1,rmss(:,:,j))
                else
                    call udgrade_nest(map2fit(:,:,j),nside1,data_low(:,:,j),nside2)
                    call udgrade_nest(fg_amp(:,:,j,1),nside1,fg_amp_low(:,:,j),nside2)
                    call udgrade_nest(rmss(:,:,j),nside1,rms_low(:,:,j),nside2)
                end if
            end do
            rms_low = sqrt(rms_low / (npix/npix2))
        else
            do j = 1, nbands
                data_low(:,:,j)   = data(:,:,j)
                fg_amp_low(:,:,j) = fg_amp(:,:,j,1)
                rms_low(:,:,j)    = rmss(:,:,j)
            end do
            indx_low = indx
        end if

        x(1) = 1.d0           

        !------------------------------------------------------------------------
        ! Sampling portion. Determine the log-likelihood, and accept based off of
        ! the improvement in the fit.
        !------------------------------------------------------------------------
        do i = 0, npix2-1
            a       = 0.d0
            sol     = indx_low(i,map_n)
            pars(1) = indx_low(i,map_n)

            ! Chi-square from the most recent Gibbs chain update
            do j = 1, nbands
                a = a + (((fg_amp_low(i,map_n,frgrnd(1)%loc) * compute_spectrum(frgrnd(1),nuz(j))) &
                      - data_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
            end do
            c = a
            if (test .eqv. .true.) then
                if (i == 350) then
                    d1        = a
                    prob(1)   = 1.d0
                    chi(1)    = a
                    tump(1)   = sol
                    accept(1) = 0.d0
                    naccept   = 0
                end if
            end if

            do l = 1, iterations

                ! Sampling from the prior
                t       = rand_normal(prior(1), prior(2))
                pars(1) = t
                b       = 0.d0

                do j = 1, nbands
                    tmp(j) = fg_amp_low(i,map_n,frgrnd(1)%loc)*compute_spectrum(frgrnd(1),nuz(j))
                    b      = b + ((tmp(j)-data_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
                end do
                b = b

                if (b < c .and. t .lt. -2.5 .and. t .gt. -3.5) then
                    sam = t
                    c   = b
                    if (test .eqv. .true.) then
                        if (i == 350) then
                            naccept  = naccept + 1
                            prob(l+1)= 1
                        end if
                    end if
                else
                    x(2) = exp(0.5d0*(c-b))
                    p = minval(x)
                    call RANDOM_NUMBER(num)
                    if (num < p) then
                        if (t .lt. -2.5 .and. t .gt. -3.5) then
                            sam = t
                            c   = b
                            if (test .eqv. .true.) then
                                if (i == 350) then
                                    naccept = naccept + 1
                                end if
                            end if
                        end if
                    end if
                    ! Metropolis testing apparatus
                    !-----------------------------
                    if (test .eqv. .true.) then
                        if (i == 350) then
                            prob(l+1)   = p
                        end if
                    end if
                end if
                if (test .eqv. .true.) then
                    if (i == 350) then 
                        chi(l+1)    = c
                        tump(l+1)   = sam
                        accept(l+1) = naccept/l
                    end if
                end if
                !-----------------------------
            end do
            if (c < a) then
                sol = sam
            else
                sol = indx_low(i,map_n)
            end if

            ! Metropolis testing apparatus
            !-----------------------------
            if (test .eqv. .true.) then
                if (i == 350) then
                    inquire(file=trim(direct) // trim(tqu(k)) // '_prob.dat',exist=exist)
                    if (exist) then
                        open(40,file = trim(direct) // trim(tqu(k)) // '_prob.dat',&
                                       status="old", position="append", action="write")
                    else
                        open(40,file = trim(direct) // trim(tqu(k)) // '_prob.dat',&
                                       status="new", action="write")
                    endif
    
                    inquire(file=trim(direct) // trim(tqu(k)) // '_chi.dat',exist=exist)
                    if (exist) then
                        open(41,file = trim(direct) // trim(tqu(k)) // '_chi.dat',&
                                       status="old", position="append", action="write")
                    else
                        open(41,file = trim(direct) // trim(tqu(k)) // '_chi.dat',&
                                       status="new", action="write")
                    endif
    
                    inquire(file=trim(direct) // trim(tqu(k)) // '_temps.dat',exist=exist)
                    if (exist) then
                        open(42,file = trim(direct) // trim(tqu(k)) // '_temps.dat',&
                                       status="old", position="append", action="write")
                    else
                        open(42,file = trim(direct) // trim(tqu(k)) // '_temps.dat',&
                                       status="new", action="write")
                    endif
    
                    inquire(file=trim(direct) // trim(tqu(k)) // '_accept.dat',exist=exist)
                    if (exist) then
                        open(43,file = trim(direct) // trim(tqu(k)) // '_accept.dat',&
                                       status="old", position="append", action="write")
                    else
                        open(43,file = trim(direct) // trim(tqu(k)) // '_accept.dat',&
                                       status="new", action="write")
                    endif
    
                    do l = 1, iterations+1
                        write(40,*) prob(l)
                        write(41,*) chi(l)
                        write(42,*) tump(l)
                        write(43,*) accept(l)
                    end do
                    close(40)
                    close(41)
                    close(42)
                    close(43)
                end if
            end if
            !-----------------------------
            if (i == 350) then
               d2 = c
            end if
            indx_sample_low(i) = sol
        end do

        if (nside1 /= nside2) then
            if (ordering == 1) then
                call udgrade_ring(indx_sample_low,nside2,indx_sample,nside1)
                call convert_nest2ring(nside2,indx_sample_low)
            else
                call udgrade_nest(indx_sample_low,nside2,indx_sample,nside1)
            end if
        else
            indx_sample = indx_sample_low
        end if

        if (trim(type) == 'synch') then 
            beta_s(:,k) = indx_sample
        else if (trim(type) == 'dust') then 
            beta_d(:,k) = indx_sample
        end if
        deallocate(data_low)
        deallocate(fg_amp_low)
        deallocate(indx_low)
        deallocate(rms_low)
        deallocate(indx_sample_low)

    end subroutine sample_index


    ! This architecture of this function has not been modified yet
    function sample_T(band, npix, map_n, sigma, T_map, nside2)
        implicit none
  
        integer(i4b), intent(in)                               :: npix, map_n, nside2
        integer(i4b)                                           :: nside1, npix2
        real(dp), intent(in)                                   :: sigma
        real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: band
        real(dp), dimension(0:npix-1,nmaps), intent(in)        :: T_map
        real(dp), dimension(0:npix-1)                          :: sample_T
        real(dp), dimension(0:npix-1,nmaps,nbands)             :: map2fit, cov
        real(dp), dimension(0:npix-1,nmaps)                    :: te
        real(dp), allocatable, dimension(:,:,:)                :: band_low, fg_amp_low, rms_low
        real(dp), allocatable, dimension(:,:)                  :: T_low
        real(dp), allocatable, dimension(:)                    ::  sample_T_low
        real(dp), dimension(nbands)                            :: signal, tmp
        real(dp), dimension(2)                                 :: x
        real(dp)                                               :: a, b, c, num, sam, t, p, sol, naccept_t_d

        te      = T_map
        map2fit = band
        cov     = rmss*rmss

        nside1 = npix2nside(npix)
        if (nside1 == nside2) then
            npix2 = npix
        else
            npix2 = nside2npix(nside2)
        end if
        allocate(band_low(0:npix2-1,nmaps,nbands),fg_amp_low(0:npix2-1,nmaps,nbands))
        allocate(T_low(0:npix2-1,nmaps),rms_low(0:npix2-1,nmaps,nbands))
        allocate(sample_T_low(0:npix2-1))

        if (nside1 /= nside2) then 
            if (ordering == 1) then
                call udgrade_ring(te,nside1,T_low,nside2)
            else
                call udgrade_nest(te,nside1,T_low,nside2)
            end if
            do j = 1, nbands
                if (ordering == 1) then
                    call udgrade_ring(map2fit(:,:,j),nside1,band_low(:,:,j),nside2)
                    call convert_nest2ring(nside1,map2fit(:,:,j))
                    call udgrade_ring(fg_amp(:,:,j,1),nside1,fg_amp_low(:,:,j),nside2)
                    call convert_nest2ring(nside1,fg_amp(:,:,j,1))
                    call udgrade_ring(cov(:,:,j),nside1,rms_low(:,:,j),nside2)
                    call convert_nest2ring(nside1,rmss(:,:,j))
                else
                    call udgrade_nest(map2fit(:,:,j),nside1,band_low(:,:,j),nside2)
                    call udgrade_nest(fg_amp(:,:,j,1),nside1,fg_amp_low(:,:,j),nside2)
                    call udgrade_nest(rmss(:,:,j),nside1,rms_low(:,:,j),nside2)
                end if
            end do
            rms_low = sqrt(rms_low / (npix/npix2))
        else
            do j = 1, nbands
                band_low(:,:,j)   = band(:,:,j)
                fg_amp_low(:,:,j) = fg_amp(:,:,j,1)
                rms_low(:,:,j)    = rmss(:,:,j)
            end do
            T_low = T_map
        end if

        x(1) = 1.d0
        do i = 0, npix2-1
            if (mask(i,1) == 0.d0) then
                sample_T_low(i) = missval
               cycle
            else
                a   = 0.d0
                sol = T_low(i,map_n)

                ! Chi-square from the most recent Gibbs chain update
                do j = 1, nbands
                    a = a + (((fg_amp_low(i,map_n,j) * HI(i,1)*planck(nuz(j)*1.d9,sol)) &
                        - band_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2
                end do
                c   = a
                sam = T_d_mu
                if (test .eqv. .true. .and. iter .eq. 1) then
                    if (i == 350) then
                        prob(1)     = 1.d0
                        chi(1)      = a
                        tump(1)     = sol
                        accept(1)   = 0.d0
                        naccept_t_d = 0
                    end if
                end if

                do l = 1, iterations
                    ! Begin sampling from the prior
                    t = rand_normal(sam, sigma)
                    b = 0.d0
                    do j = 1, nbands
                        tmp(j) = fg_amp_low(i,map_n,j)*HI(i,1)*planck(nuz(j)*1.d9,t)
                        b      = b + ((tmp(j)-band_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
                    end do
                    b = b

                    if (b < c .and. t .lt. 35.d0 .and. t .gt. 10.d0) then
                        sam = t
                        c   = b
                        if (test .eqv. .true.) then
                            if (i == 350) then
                                naccept_t_d  = naccept_t_d + 1
                                prob((iter-1)*100+l+1)= 1
                            end if
                        end if
                    else
                        x(2) = exp(0.5d0*(c-b))
                        p = minval(x)
                        call RANDOM_NUMBER(num)
                        if (num < p) then
                            if (t .lt. 35.d0 .and. t .gt. 10.d0) then
                                sam = t
                                c   = b
                                if (test .eqv. .true.) then
                                    if (i == 350) then
                                        naccept_t_d = naccept_t_d + 1
                                    end if
                                end if
                            end if
                        end if
                        ! Metropolis testing apparatus
                        !-----------------------------
                        if (test .eqv. .true.) then
                            if (i == 350) then
                                prob((iter-1)*100+l+1)   = p
                            end if
                        end if
                    end if
                    if (test .eqv. .true.) then
                        if (i == 350) then 
                            chi((iter-1)*100+l+1)    = c
                            tump((iter-1)*100+l+1)   = sam
                            accept((iter-1)*100+l+1) = naccept_t_d/((iter-1)*100+l+1)
                        end if
                    end if
                    !-----------------------------
                end do
                if (c < a) then
                    sol = sam
                else
                    sol = T_low(i,map_n)
                end if
            end if
            sample_T_low(i) = sol
        end do
        if (nside1 /= nside2) then
            if (ordering == 1) then
                call udgrade_ring(sample_T_low,nside2, sample_T,nside1)
                call convert_nest2ring(nside2, sample_T_low)
            else
                call udgrade_nest(sample_T_low,nside2, sample_T,nside1)
            end if
        else
             sample_T =  sample_T_low
        end if

        deallocate(band_low)
        deallocate(fg_amp_low)
        deallocate(T_low)
        deallocate(rms_low)
        deallocate(sample_T_low)

    end function sample_T
    ! ------------------------------------------------------------

    subroutine sample_joint_amp(npix, map_n, method)
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
        integer(i4b),              intent(in)     :: npix, map_n
        character(len=*),          intent(in)     :: method
        real(dp), allocatable, dimension(:,:,:)   :: T_nu, T_nu_T, covar, A_1, A_2
        real(dp), allocatable, dimension(:,:)     :: A, c_1, c_2, dats, norm
        real(dp), allocatable, dimension(:,:)     :: mat_l, mat_u, A_inv
        real(dp), allocatable, dimension(:)       :: b, c, d, samp, rand
        integer(i4b)                              :: x, y, z, nfit1, nfit2, w

        ! Load which components to jointly fit for arary allocation
        ! vv These will not change based off of components
        x = npix
        y = 0
        z = nbands

        do i = 1, size(joint_comps)
            if (joint_comps(i) == 'synch') then
                y = y + npix
            else if (joint_comps(i) == 'dust') then
                y = y + npix
            else if (joint_comps(i) == 'template01') then
                ! Must have at most nbands -1 templates to avoid degeneracy
                ! Count how many bands are not being fit
                nfit1 = 0
                do j = 1, nbands
                    if (j_corr01(j) .eqv. .true.) then
                        nfit1 = nfit1 + 1
                    end if
                end do
                y = y + nfit1
            else if (joint_comps(i) == 'template02') then
                ! Must have at most nbands -1 templates to avoid degeneracy
                ! Count how many bands are not being fit
                nfit2 = 0
                do j = 1, nbands
                    if (j_corr02(j) .eqv. .true.) then
                        nfit2 = nfit2 + 1
                    end if
                end do
                y = y + nfit2
            end if
        end do

        allocate(T_nu(x,y,z),T_nu_T(y,x,z),dats(x,z))
        allocate(A_1(y,x,z),A_2(y,y,z))
        allocate(A(y,y),b(y),c(y),d(y),A_inv(y,y))
        allocate(mat_l(y,y),mat_u(y,y))
        allocate(samp(y),rand(y), norm(y,y))
        allocate(covar(x,x,z),c_1(x,z),c_2(y,z))

        ! Initialize arrays
        covar(:,:,:)      = 0.d0
        T_nu(:,:,:)       = 0.d0
        A_1(:,:,:)        = 0.d0
        A_2(:,:,:)        = 0.d0
        A(:,:)            = 0.d0
        A_inv(:,:)        = 0.d0
        b(:)              = 0.d0
        c(:)              = 0.d0
        d(:)              = 0.d0
        rand(:)           = 0.d0
        samp(:)           = 0.d0
        norm(:,:)         = 0.d0
        dats(:,:)         = 0.d0
        mat_l(:,:)        = 0.d0
        mat_u(:,:)        = 0.d0
        c_1(:,:)          = 0.d0
        c_2(:,:)          = 0.d0

        ! Fill data and covariance arrays
        do i=1, x
            do j=1,z
                covar(i,i,j) = rmss(i-1,map_n,j)**2
                dats(i,j)    = maps(i-1,map_n,j)
            end do
        end do

        ! Fill template matrix
        w = 0 
        do m = 1, size(joint_comps)
            if (joint_comps(m) == 'synch') then
                do i = 1, x
                    frgrnd(1)%p(1) = beta_s(i-1,map_n)
                    do j = 1, z
                        T_nu(i,w+i,j) = compute_spectrum(frgrnd(1),bands(j)%nu)
                        if (i == 1) then
                            
                        end if
                    end do
                end do
                w = w + x
            else if (joint_comps(m) == 'dust') then
                do i = 1, x
                    frgrnd(3)%p(1) = beta_d(i-1,map_n)
                    frgrnd(3)%p(2) = T_d(i-1,map_n)
                    do j = 1, z
                        T_nu(i,w+i,j) = compute_spectrum(frgrnd(3),bands(j)%nu)
                    end do
                end do
                w = w + x
            end if
        end do

        do m = 1, size(joint_comps)
            if (joint_comps(m) == 'template01') then
                do i = 1, x
                    l = 1
                    do j = 1, z
                        if (j_corr01(j) .eqv. .true.) then
                            T_nu(i,w+l,j) = template_01(i-1,map_n)
                            l = l + 1
                        end if
                    end do
                end do
                w = w + l  
            else if (joint_comps(m) == 'template02') then
                do i = 1, x
                    l = 1
                    do j = 1, z
                        if (j_corr02(j) .eqv. .true.) then
                            T_nu(i,w+l,j) = template_01(i-1,map_n)
                            l = l + 1
                        end if
                    end do
                end do
                w = w + l
            end if
        end do

        ! Computing the LHS and RHS of the linear equation
        do j=1, nbands
            T_nu_T(:,:,j) = transpose(T_nu(:,:,j))
            c_1(:,j)      = matmul(inv(covar(:,:,j)),dats(:,j))
            c_2(:,j)      = matmul(T_nu_T(:,:,j),c_1(:,j))
            A_1(:,:,j)    = matmul(T_nu_T(:,:,j),inv(covar(:,:,j)))
            A_2(:,:,j)    = matmul(A_1(:,:,j),T_nu(:,:,j)) 
            A(:,:)        = A(:,:) + A_2(:,:,j)
            c(:)          = c(:) + c_2(:,j)
        end do

        ! Computation
        if (trim(method) == 'cholesky') then
            if (mod(iter,output_iter) .EQ. 0) then
                write(*,*) 'Joint sampling using Cholesky Decomp'
            end if
            call cholesky_decomp(A,mat_l,y)
            mat_u  = transpose(mat_l)
            call forward_sub(mat_l,d,c)
            call backward_sub(mat_u,b,d)
        else if (trim(method) == 'cg') then
            if (mod(iter,output_iter) .EQ. 0) then
                write(*,*) 'Joint sampling using CG'
            end if
            call compute_cg(A,b,c,y)
        else if (trim(method) == 'lu') then
            if (mod(iter,output_iter) .EQ. 0) then
                write(*,*) 'Joint sampling using LU Decomp'
            end if
            call LUDecomp(A,mat_l,mat_u,y)
            call forward_sub(mat_l,d,c)
            call backward_sub(mat_u,b,d)
        end if

        ! Draw a sample by cholesky decompsing A^-1, and multiplying 
        ! the subsequent lower triangular by a vector of random numbers
        A_inv = inv(A)
        call cholesky_decomp(A_inv,norm,y)
        do i = 1, y
            rand(i) = rand_normal(0.d0,1.d0)
        end do
        samp  = matmul(norm,rand)
        b     = b!  + samp

        ! Output amplitudes to the appropriate variables
        do m = 1, size(joint_comps)
            if (joint_comps(m) == 'synch') then
                do i = 1, x
                    fg_amp(i-1,map_n,frgrnd(1)%loc,1) = b(i)
                end do
                b = b(x+1:)
            else if (joint_comps(m) == 'dust') then
                do i = 1, x
                    fg_amp(i-1,map_n,frgrnd(3)%loc,3) = b(i)
                end do
                b = b(x+1:)
            end if
        end do
        do m = 1, size(joint_comps)
            if (joint_comps(m) == 'template01') then
                l = 1
                do while (l .lt. (nfit1))
                    do j= 1, z
                        if (j_corr01(j) .eqv. .true.) then
                            fg_amp(:,map_n,j,2) = b(l)
                            l = l + 1
                        end if
                    end do
                end do
                b = b(l+nfit1:)
            else if (joint_comps(m) == 'template02') then
                l = 1
                do while (l .lt. (nfit2))
                    do j= 1, z
                        if (j_corr02(j) .eqv. .true.) then
                            fg_amp(:,map_n,j,2) = b(l)
                            l = l + 1
                        end if
                    end do
                end do
                b = b(l+nfit2:)
            end if
        end do

        ! Sure to deallocate all arrays here to free up memory
        deallocate(A_1)
        deallocate(A_2)
        deallocate(A)
        deallocate(A_inv)
        ! deallocate(b)
        deallocate(c)
        deallocate(c_1)
        deallocate(c_2)
        deallocate(covar)
        deallocate(dats)
        deallocate(d)
        deallocate(norm)
        deallocate(mat_l)
        deallocate(mat_u)
        deallocate(rand)
        deallocate(samp)
        deallocate(T_nu)
        deallocate(T_nu_T)

    end subroutine sample_joint_amp

    subroutine write_maps(nm)
        implicit none

        integer(i4b), intent(in)    :: nm
        real(dp), dimension(npix,1) :: map

        write(iter_str, '(i0.5)') iter
        if (output_fg .eqv. .true.) then
            do j = 1, nbands
                title = trim(direct) // trim(bands(j)%label) // '_dust_fit_'// trim(tqu(nm)) & 
                        // '_' // trim(iter_str) // '.fits'
                map(:,1)   = dust_map(:,nm,j)
                call write_bintab(map,npix,1, header, nlheader, trim(title))
                title = trim(direct) // trim(bands(j)%label) // '_synch_amplitude_' //  trim(tqu(nm)) &
                        // '_' // trim(iter_str) // '.fits'
                map(:,1)   = fg_amp(:,nm,j,1)
                call write_bintab(map,npix,1, header, nlheader, trim(title))
            end do
        else 
            title = trim(direct) // trim(bands(frgrnd(1)%loc)%label) // '_synch_amplitude_' //  trim(tqu(nm)) &
                    // '_' // trim(iter_str) // '.fits'
            map(:,1)   = fg_amp(:,nm,frgrnd(1)%loc,1)
            call write_bintab(map,npix,1, header, nlheader, trim(title))
        end if
        do j = 1, nbands
            title = trim(direct) // trim(bands(j)%label) // '_residual_' // trim(tqu(nm)) & 
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
                chi_map(i,nm) = chi_map(i,nm) + (maps(i,nm,j) - synch_map(i,nm,j) - dust_map(i,nm,j))**2.d0/rmss(i,nm,j)**2.d0
            end do
        end do
        chi_map(:,nm) = chi_map(:,nm)/(nbands+3)
        title = trim(direct) // 'chisq_' // trim(tqu(nm)) // '_' // trim(iter_str) // '.fits'
        map(:,1)   = chi_map(:,nm)
        call write_bintab(map,npix,1, header, nlheader, trim(title))

    end subroutine write_maps


    subroutine write_data
        implicit none

        title = trim(direct) // 'pixel_100_A_d.dat'
        inquire(file=title,exist=exist)
        if (exist) then
            open(30,file=title, status="old",position="append", action="write")
        else
            open(30,file=title, status="new", action="write")
        endif
        write(30,*) dust_map(100,k,frgrnd(1)%loc)
        close(30)

        title = trim(direct) // 'pixel_100_A_s.dat'
        inquire(file=title,exist=exist)
        if (exist) then
            open(31,file=title, status="old",position="append", action="write")
        else
            open(31,file=title, status="new", action="write")
        endif
        write(31,*) fg_amp(100,k,frgrnd(1)%loc,1)
        close(31)

        title = trim(direct) // 'pixel_100_beta_s.dat'
        inquire(file=title,exist=exist)
        if (exist) then
            open(32,file=title, status="old",position="append", action="write")
        else
            open(32,file=title, status="new", action="write")
        endif
        write(32,*) beta_s(100,k)
        close(32)

        title = trim(direct) // 'total_chisq.dat'
        inquire(file=title,exist=exist)
        if (exist) then
            open(33,file=title, status="old",position="append", action="write")
        else
            open(33,file=title, status="new", action="write")
        endif
        call compute_chisq(fg_amp,k)
        write(33,*) chisq
        close(33)

        inquire(file=trim(direct) // 'dust_' // trim(tqu(k)) // '_amplitudes.dat',exist=exist)
        if (exist) then
            open(34,file = trim(direct) // 'dust_' // trim(tqu(k)) // '_amplitudes.dat', status="old", &
                        position="append", action="write")
        else
            open(34,file = trim(direct) // 'dust_' // trim(tqu(k)) // '_amplitudes.dat', status="new", action="write")
        endif
        write(34,'(6(E17.8))') dust_amps
        close(34)

    end subroutine write_data
  
! Still need to rewrite vv

    subroutine compute_chisq(amp,map_n)
      use healpix_types
      implicit none
      real(dp), dimension(0:npix-1,nmaps,nbands,3), intent(in)   :: amp
      integer(i4b), intent(in)                                   :: map_n
      real(dp)                                                   :: s, signal
      integer(i4b)                                               :: m,n
  
      chisq = 0.d0
      do i = 0, npix-1
        do j = 1, nbands
            s = 0.d0
            signal = amp(i,map_n,j,1)*mask(i,1) + amp(i,map_n,j,2)*template_01(i,map_n)*mask(i,1)
            s = s + signal
            chisq = chisq + (((maps(i,map_n,j) - s)**2))/(rmss(i,map_n,j)**2)*mask(i,1)
        end do
      end do 
      chisq = chisq/(sum(mask(:,1))+nbands+3) ! n-1 dof, npix + nbands + A_s + A_dust + A_cmb + \beta_s
    end subroutine compute_chisq

  end program dust_fit
