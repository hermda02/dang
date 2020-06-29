program dust_fit
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
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
    
    ! Constants
    real(dp)           :: k_B     = 1.3806503d-23
    real(dp)           :: h       = 1.0545726691251021d-34 * 2.d0 * pi
    real(dp)           :: c       = 2.99792458d8
    real(dp)           :: T_CMB   = 2.7255d0
  
    integer(i4b)       :: i, j, k, l, iter, npix, nside, nmaps, ordering, loc
    integer(i4b)       :: beta_samp_nside, nlheader, niter, nbands, nfgs, iterations
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
    real(dp)                                     :: chisq, temp_norm_01
    character(len=80), dimension(180)            :: header
    character(len=80), dimension(3)              :: tqu
    character(len=80), allocatable, dimension(:) :: bands
    character(len=5)                             :: iter_str
    logical(lgt), allocatable, dimension(:)      :: j_corr

    real(dp), allocatable, dimension(:,:)        :: mat_test, mat_l, mat_u
    real(dp), allocatable, dimension(:)          :: x, b, d

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

    niter             = 1000       ! # of MC-MC iterations
    iterations        = 100        ! # of iterations in the samplers
    output_iter       = 100        ! Output maps every <- # of iterations
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
    allocate(template_01(0:npix-1,nmaps), template_02(0:npix-1,nmaps), dust_amps(nbands),j_corr(nbands))
    allocate(maps(0:npix-1,nmaps,nbands), rmss(0:npix-1,nmaps,nbands), nodust(0:npix-1,nmaps,nbands))
    allocate(mask(0:npix-1,1), res(0:npix-1,nmaps,nbands), chi_map(0:npix-1,nmaps))
    allocate(fg_amp(0:npix-1,nmaps,nbands,nfgs), beta_s(0:npix-1,nmaps))
    allocate(T_d(0:npix-1,nmaps), beta_d(0:npix-1,nmaps), HI(0:npix-1,nmaps))
    allocate(cmb_map(0:npix-1,nmaps,nbands), dust_map(0:npix-1,nmaps,nbands), synch_map(0:npix-1,nmaps,nbands))
    allocate(nuz(nbands), bands(nbands), par(2))
    allocate(map(0:npix-1,nmaps))
    allocate(rms(0:npix-1,nmaps))
    !----------------------------------------------------------------------------------------------------------
    beta_s     = -3.10d0    ! Synchrotron beta initial guess
    beta_d     = 1.60d0     ! Dust beta initial guess

    bands(1)   = ('norm_pol_020_')
    bands(2)   = ('norm_pol_045_')
    bands(3)   = ('norm_pol_070_')
    bands(4)   = ('norm_pol_100_')
    bands(5)   = ('norm_pol_200_')
    bands(6)   = ('norm_pol_353_')

    nuz(1)     = 20.0d0
    nuz(2)     = 45.0d0
    nuz(3)     = 70.0d0
    nuz(4)     = 100.0d0
    nuz(5)     = 200.0d0
    nuz(6)     = 353.0d0

    loc        = minloc(abs(nuz-nu_ref_s),1)
    !----------------------------------------------------------------------------------------------------------
    ! Read maps

    do j = 1, nbands
        call read_bintab('data/test_data/norm_pol/' // trim(bands(j)) // 'rms_n0004.fits',&
        rms,npix,nmaps,nullval,anynull,header=header)
        rmss(:,:,j) = rms
        call read_bintab('data/test_data/norm_pol/' // trim(bands(j)) // 'noised_n0004.fits', &
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
    j_corr(1) = .false.
    j_corr(2) = .true.
    j_corr(3) = .true.
    j_corr(4) = .true.
    j_corr(5) = .true.
    j_corr(6) = .true.      
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Metropolis-Hastings testing things
    allocate(chi(iterations*niter+1), tump(iterations*niter+1), accept(iterations*niter+1), prob(iterations*niter+1))
    !----------------------------------------------------------------------------------------------------------
    temp_norm_01 = maxval(template_01)
    template_01 = template_01/temp_norm_01
    !----------------------------------------------------------------------------------------------------------
    ! Calculation portion
    !----------------------------------------------------------------------------------------------------------

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

        call compute_chisq(fg_amp,k)

        write(*,*) 'Initial Chisq = ', chisq

        do iter = 1, niter
        
            ! write(*,*) 'Iteration', iter
            ! write(*,*) '-----------------------'
            ! write(*,*) ''


            ! write(*,*) 'Jointly Sampling Amplitudes' 
            call sample_joint_amp(npix,k,'cholesky')  ! Method possibilities are 'cg', 'LU', and 'cholesky'
            dust_amps = fg_amp(0,k,:,2)

            call compute_chisq(fg_amp,k)

            write(*,fmt='(i6, a, f10.3)') iter, ' - chisq: ' , chisq

            ! -------------------------------------------------------------------------------------------------------------------
            ! Extrapolating A_synch to bands
            do i = 0, npix-1
                par(1) = beta_s(i,k)
                do j = 1, nbands
                    fg_amp(i,k,j,1) = fg_amp(i,k,loc,1)*compute_spectrum('synch',nuz(j),par)
                end do
            end do
            synch_map(:,k,:)        = fg_amp(:,k,:,1)

            ! -------------------------------------------------------------------------------------------------------------------
            ! Applying A_d to make dust maps
            do j = 1, nbands
                dust_map(:,k,j) = fg_amp(:,k,j,2)*template_01(:,k)
            end do
            ! -------------------------------------------------------------------------------------------------------------------

            ! nodust = maps-dust_map
            ! -------------------------------------------------------------------------------------------------------------------
            ! call sample_index(nodust,'synch',beta_samp_nside,k)
            ! do i = 0, npix-1
            !    par(1) = beta_s(i,k)
            !    do j = 1, nbands
            !        fg_amp(i,k,j,1) = fg_amp(i,k,loc,1)*compute_spectrum('synch',nuz(j),par)
            !    end do
            ! end do
            ! synch_map(:,k,:)        = fg_amp(:,k,:,1)
            ! -------------------------------------------------------------------------------------------------------------------

            res       = maps - synch_map - dust_map

            call compute_chisq(fg_amp,k)

            if (mod(iter, 1) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, f10.3, a, f7.3, a, f8.4, a, 6e10.3)')&
                 iter, " - chisq: " , chisq, " - A_s: ",&
                 fg_amp(100,k,loc,1),  " - beta_s: ",&
                 sum(beta_s(:,k))/npix, ' - A_d: ', dust_amps*temp_norm_01
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

    function compute_spectrum(fg_type,freq, params)
        implicit none

        character(len=*), intent(in)       :: fg_type
        real(dp),         intent(in)       :: freq
        real(dp), dimension(:), intent(in) :: params
        real(dp)                           :: y, rj_cmb ! Conversion factor from K_{RJ} -> K_{CMB}
        real(dp)                           :: z, compute_spectrum

        ! Simple function for computing foreground spectra

        y = h*(freq*1.0d9) / (k_B*T_CMB)
        rj_cmb = ((exp(y)-1)**2.d0)/(y**2.d0*exp(y))

        if (trim(fg_type) == 'synch') then
            compute_spectrum = (freq/nu_ref_s)**params(1)!*rj_cmb
        else if (trim(fg_type) == 'dust') then
            z = h / (k_B*params(2))
            compute_spectrum = (exp(z*nu_ref_d*1d9)-1.d0) / (exp(z*freq*1d9)-1.d0) * (freq/nu_ref_d)**(params(1)+1.d0)!*rj_cmb
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
            ! Following code is for sampling and only accepting the best fit solution
            ! -----------------------------------------------------------------------
            ! chi_00 = 0.d0
            ! do i = 0, npix-1
            !     chi_00 = chi_00 + (data(i)-amp*template(i))**2.d0/cov(i)*t_mask(i)
            ! end do
    
            ! do l = 1, iterations
            !     chi_0 = 0.d0
            !     sam = amp + rand_normal(0.d0,1.d0)/sqrt(norm)
            !     do i = 0, npix-1
            !         chi_0 = chi_0 + (data(i)-amp*template(i))**2.d0/cov(i)*t_mask(i)
            !     end do
            !     if (chi_0 < chi_0) then
            !         amp = sam
            !         chi = chi_0
            !     end if
            ! end do
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
                spec          = compute_spectrum(type,nuz(j),pars)
                sum1          = sum1 + (data(i,j)*spec)/cov(i,j)
                sum2          = sum2 + (spec)**2.d0/cov(i,j)
                norm(i)       = norm(i) + ((spec)**2.d0)/cov(i,j)
            end do
            norm(i)           = norm(i)/nbands
            amp               = sum1/sum2
    
            sample_spec_amp(i) = amp + rand_normal(0.d0,1.d0)/sqrt(norm(i))

            ! Following code is for sampling and only accepting the best fit solution
            ! -----------------------------------------------------------------------
            ! chi_00 = 0.d0
            ! do j = 1, nbands
            !     chi_00 = chi_00 + (band(i,map_n,j)-(A(i,map_n))*compute_spectrum(type,nuz(j),i))**2.d0/cov(i,map_n,j)
            ! end do
            ! chi = 0.d0
            ! do j = 1, nbands
            !     chi = chi + (band(i,map_n,j)-(sum1/sum2)*compute_spectrum(type,nuz(j),i))**2.d0/cov(i,map_n,j)
            ! end do
            ! ! write(*,*) 
            ! do l = 1, iterations
            !     chi_0 = 0.d0
            !     sam = amp + rand_normal(0.d0,1.d0)/sqrt(norm(i,map_n))
            !     do j = 1, nbands
            !         chi_0 = chi_0 + (band(i,map_n,j)-(sam*compute_spectrum(type,nuz(j),i)))**2.d0/cov(i,map_n,j)
            !     end do
            !     if (chi_0 < chi_0) then
            !         amp = sam
            !         chi = chi_0
            !     end if
            ! end do
            ! if (chi_0 < chi_00) then
            !     sample_s_amp(i) = amp
            ! else
            !     sample_s_amp(i) = A(i,map_n)
            ! end if
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
                a = a + (((fg_amp_low(i,map_n,loc) * compute_spectrum(type,nuz(j),pars)) &
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
                    tmp(j) = fg_amp_low(i,map_n,loc)*compute_spectrum(type,nuz(j),pars)
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
        real(dp), allocatable, dimension(:,:)     :: A, c_1, c_2, dats, norm, A_inv
        real(dp), allocatable, dimension(:,:)     :: mat_l, mat_u, samp_test
        real(dp), allocatable, dimension(:)       :: b, c, d, samp, rand
        real(dp)                                  :: chi0, chi_prop
        integer(i4b)                              :: x, y, z, nskip

        nskip = 0

        do j =1, nbands
           if (j_corr(j) .eqv. .false.) then
              nskip = nskip + 1
           end if
        end do

        chi0 = chisq

        x = npix
        y = npix+nbands-nskip
        z = nbands

        allocate(T_nu(x,y,z),T_nu_T(y,x,z),dats(x,z))
        allocate(A_1(y,x,z),A_2(y,y,z))
        allocate(A(y,y),b(y),c(y),d(y))
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
        do i=1,x
            l = 1
            do j=1, z
                T_nu(i,i,j) = (nuz(j)/nu_ref_s)**beta_s(i-1,map_n)
                if (j_corr(j) .eqv. .true.) then
                    T_nu(i,x+l,j) = template_01(i-1,map_n)
                    l = l + 1
                end if
            end do
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
        else if (trim(method) == 'LU') then
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

        samp  = matmul(norm,rand) !matmul(norm,rand)

        b = b + samp

        ! Output amplitudes to the appropriate variables
        do i =1, x
            fg_amp(i-1,map_n,loc,1) = b(i)
        end do
        l = 1
        do while (l .lt. (nbands-nskip))
            do j= 1, z
                if (j_corr(j) .eqv. .true.) then
                    fg_amp(:,map_n,j,2) = b(x+l)
                    l = l + 1
                end if
            end do
        end do

        ! Sure to deallocate all arrays here to free up memory
        deallocate(A_1)
        deallocate(A_2)
        deallocate(A)
        deallocate(A_inv)
        deallocate(b)
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

    subroutine cholesky_decomp(mat,low,n)
        implicit none
        real(dp), dimension(:,:), intent(in)  :: mat
        real(dp), dimension(:,:), intent(out) :: low
        integer(i4b),             intent(in)  :: n
        integer(i4b)                          :: ip
        real(dp)                              :: s

        low(:,:)   = 0.d0

        do i = 1, n
           low(i,i) = sqrt(mat(i,i) - dot_product(low(i,1:i-1),low(i,1:i-1)) )
           do j = i+1, n
              low(j,i) = (mat(j,i) - dot_product(low(j,1:i-1),low(i,1:i-1)))/low(i,i)
           end do
        end do

        ! This code would be used for the LDU decomp:
        ! -------------------------------------------
        ! do i = 1, n
        !     mat_s(i,i) = low(i,i)
        ! end do

        ! low  = matmul(low,inv(mat_s))
        ! diag = mat_s**2

        ! deallocate(mat_s)
        ! -------------------------------------------

    end subroutine cholesky_decomp

    function inv(A) result(Ainv)
        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(size(A,1),size(A,2)) :: Ainv
      
        real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info
      
        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI
      
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
      
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)
      
        if (info /= 0) then
           stop 'Matrix is numerically singular!'
        end if
      
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)
      
        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if
    end function inv

    subroutine LUDecomp(A,L,U,n)
        
        ! Using Doolittle's method for LU decomposition of matrix A
        ! L and U can then be used to solve the matrix equation Ax = b by solving
        ! Lv = b (using forward substitution), Ux = v (backwards substitution).
        
        implicit none        
        real(dp), dimension(:,:), intent(in)  :: A
        real(dp), dimension(:,:), intent(out) :: L
        real(dp), dimension(:,:), intent(out) :: U
        integer(i4b), intent(in)              :: n
        real(dp)                              :: sum

        L = 0.d0
        U = 0.d0

        do i = 1, n
            do m = i, n
                sum = 0
                do j = 1, i
                    sum = sum + (L(i,j) * U(j,m))
                U(i,m) = A(i,m) - sum
                end do
            end do
            do m = i, n
                if (i == m) then
                    L(i,i) = 1
                else
                    sum = 0
                    do j = 1, i
                        sum = sum + (L(m,j)*U(j,i))
                    end do
                    L(m,i) = (A(m,i) - sum)/U(i,i)
                end if
            end do
        end do

    end subroutine LUDecomp

    subroutine forward_sub(L,xf,bf)
        ! Forward substitution to solve the matrix equation Lx=b
        implicit none
        real(dp), dimension(:,:), intent(in)  :: L
        real(dp), dimension(:),   intent(in)  :: bf
        real(dp), dimension(:),   intent(out) :: xf
        integer(i4b)                          :: n,m

        n = size(bf)

        do i = 1, n
            xf(i) = bf(i)
            do m = 1, i-1
                xf(i) = xf(i)-l(i,m)*xf(m) 
            end do
            xf(i) = xf(i)/L(i,i)
        end do

    end subroutine forward_sub

    subroutine backward_sub(U,xb,bb)
        ! Backward substitution to solve the matrix equation Ux=b
        implicit none
        real(dp), dimension(:,:), intent(in)  :: U
        real(dp), dimension(:),   intent(in)  :: bb
        real(dp), dimension(:),   intent(out) :: xb
        integer(i4b)                          :: n, m

        n = size(bb)
        xb(n) = bb(n)/U(n,n)

        do i = n-1, 1, -1
            xb(i) = bb(i)
            do m = n, i+1, -1
                xb(i) = xb(i) - U(i,m)*xb(m)
            end do
            xb(i) = xb(i)/U(i,i)
        end do

    end subroutine backward_sub

    subroutine compute_cg(A,x,b,n)
        
        ! Implementation of the canned algorithm (B2) outlined in Jonathan Richard Shewuck (1994)
        ! "An introduction to the Conjugate Gradient Method Without the Agonizing Pain"

        implicit none

        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(:), intent(in)   :: b
        real(dp), dimension(:), intent(out)  :: x
        integer(i4b), intent(in)             :: n
        real(dp), allocatable, dimension(:)  :: r, q, d
        real(dp)                             :: epsil, alpha, beta, delta_0
        real(dp)                             :: delta_old, delta_new
        integer(i4b)                         :: i_max

        allocate(r(n),q(n),d(n))

        x(:) = 0.0d0
        i_max = 10

        i = 0
        epsil = 1.0d-16

        r = b - matmul(A,x)
        d = r
        delta_new = sum(r*r)
        delta_0   = delta_new

        do while( (i .lt. i_max) .and. (delta_new .gt. (epsil**2)*delta_0))
            q = matmul(A,d)
            alpha = delta_new/(sum(d*q))
            x = x + alpha*d
            if (mod(i,50) == 0) then
                r = b - matmul(A,x)
            else
                r = r - alpha*q
            end if
            delta_old = delta_new
            delta_new = sum(r*r)
            beta = delta_new/delta_old
            d = r + beta*d
            i = i + 1
        end do

        deallocate(r)
        deallocate(q)
        deallocate(d)

    end subroutine compute_cg

    subroutine write_maps(nm)
        implicit none

        integer(i4b), intent(in)    :: nm
        real(dp), dimension(npix,1) :: map

        write(iter_str, '(i0.5)') iter
        if (output_fg .eqv. .true.) then
            do j = 1, nbands
                title = trim(direct) // trim(bands(j)) // 'dust_fit_'// trim(tqu(nm)) & 
                        // '_' // trim(iter_str) // '.fits'
                map(:,1)   = dust_map(:,nm,j)
                call write_bintab(map,npix,1, header, nlheader, trim(title))
                title = trim(direct) // trim(bands(j)) // 'synch_amplitude_' //  trim(tqu(nm)) &
                        // '_' // trim(iter_str) // '.fits'
                map(:,1)   = fg_amp(:,nm,j,1)
                call write_bintab(map,npix,1, header, nlheader, trim(title))
            end do
        else 
            title = trim(direct) // trim(bands(loc)) // 'synch_amplitude_' //  trim(tqu(nm)) &
                    // '_' // trim(iter_str) // '.fits'
            map(:,1)   = fg_amp(:,nm,loc,1)
            call write_bintab(map,npix,1, header, nlheader, trim(title))
        end if
        do j = 1, nbands
            title = trim(direct) // trim(bands(j)) // 'residual_' // trim(tqu(nm)) & 
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
        write(30,*) dust_map(100,k,loc)
        close(30)

        title = trim(direct) // 'pixel_100_A_s.dat'
        inquire(file=title,exist=exist)
        if (exist) then
            open(31,file=title, status="old",position="append", action="write")
        else
            open(31,file=title, status="new", action="write")
        endif
        write(31,*) fg_amp(100,k,loc,1)
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

    subroutine likelihood(s_amp,d_amp,bet,map_n)
        implicit none

        real(dp),                    intent(in) :: s_amp, bet
        real(dp), dimension(nbands), intent(in) :: d_amp
        integer(i4b),                intent(in) :: map_n
        real(dp), dimension(nbands)             :: actual
        real(dp), dimension(200)                :: A_s, A_s_ln, A_s_like
        real(dp), dimension(200)                :: A_d, A_d_ln, A_d_like
        real(dp), dimension(2)                  :: pars
        real(dp)                                :: del_a_s, del_a_d, signal
        integer(i4b)                            :: x, y, z

        write(iter_str, '(i0.5)') iter
        del_a_s = s_amp/100.d0

        do x = 1, 200
            A_s(x) = (x*del_a_s)! + (49.d0/100)*s_amp
        end do
        write(*,*) s_amp

        pars(1) = bet

        do x = 1, 200
            signal = 0.d0
            do z = 1, nbands
                signal    = A_s(x)*compute_spectrum('synch',nuz(j),pars) + fg_amp(100,map_n,z,2)*template_01(100,map_n)
                A_s_ln(x) = A_s_ln(x) + (((maps(100,map_n,z) - signal)**2.d0)/(rmss(100,map_n,z)**2))
            end do
            A_s_like(x) = exp((-0.5d0)*A_s_ln(x))
        end do 

        open(50,file=trim(direct) // 'likelihood_A_s_' // trim(iter_str) //'.dat')
        do x = 1, 200
            write(50,*) A_s_like(x)
        end do 
        close(50)

        open(51,file=trim(direct) // 'A_s_' // trim(iter_str) //'.dat')
        do x = 1, 200
            write(51,*) A_s(x)
        end do
        close(51)

        open(52,file=trim(direct) // 'ln_A_s_' // trim(iter_str) //'.dat')
        do x = 1, 200
            write(52,*) A_s_ln(x)
        end do 
        close(52)

        actual(1) = 1.644429d-2
        actual(2) = 8.39923d-3
        actual(3) = 1.19689d-3
        actual(4) = 5.890824d-2
        actual(5) = 2.9665593d-1

        do z = 1, nbands
            
            ! write(*,*) d_amp(z)
            ! del_a_d = d_amp(z)/100.d0

            ! do x = 1, 200
            !     A_d(x) = (x*del_a_d) !+ (49.d0/100)*d_amp(z)
            ! end do

            write(*,*) d_amp(z)
            del_a_d = actual(z)/100.d0

            do x = 1, 200
                A_d(x) = (x*del_a_d)! + (49.d0/100)*actual(z)
            end do

            do x = 1, 200
                signal = 0.d0
                do y = 0, npix-1
                    pars(1) = beta_s(y,map_n)
                    signal    = fg_amp(y,map_n,loc,map_n)*compute_spectrum('synch',nuz(j),pars) + A_d(x)*template_01(y,map_n)
                    A_d_ln(x) = A_d_ln(x) + (((maps(y,map_n,z) - signal)**2.d0)/(rmss(y,map_n,z)**2))
                end do
                A_d_like(x) = exp((-0.5d0)*A_d_ln(x))
            end do 

            open(50,file=trim(direct) //  trim(bands(z)) // 'likelihood_A_d_' // trim(iter_str) //'.dat')
            do x = 1, 200
                write(50,*) A_d_like(x)
            end do 
            close(50)

            open(51,file=trim(direct) // trim(bands(z)) //'A_d_' // trim(iter_str) //'.dat')
            do x = 1, 200
                write(51,*) A_d(x)
            end do
            close(51)

            open(52,file=trim(direct) // trim(bands(z)) // 'ln_A_d_' // trim(iter_str) //'.dat')
            do x = 1, 200
                write(52,*) A_d_ln(x)
            end do 
            close(52)
        end do

    end subroutine likelihood
  end program dust_fit
