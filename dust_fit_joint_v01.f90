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
  
    integer(i4b)       :: i, j, k, l, iter, npix, nside, nmaps, ordering, loc, j_corr
    integer(i4b)       :: beta_samp_nside, nlheader, niter, nbands, iterations
    integer(i4b)       :: output_iter, like_iter
    real(dp)           :: nullval
    real(dp)           :: missval = -1.6375d30
    logical(lgt)       :: anynull, double_precision, test, exist, output_fg
  
    character(len=128) :: dust_file, mask_file, arg1, synch_file
    character(len=128) :: mapfile, title, direct
  
    real(dp), allocatable, dimension(:,:,:,:)    :: fg_amp
    real(dp), allocatable, dimension(:,:,:)      :: maps, rmss, model, res
    real(dp), allocatable, dimension(:,:,:)      :: dust_free, synch_free
    real(dp), allocatable, dimension(:,:,:)      :: synch_map, dust_map, cmb_map
    real(dp), allocatable, dimension(:,:)        :: dust_temp, synch_temp, map, rms
    real(dp), allocatable, dimension(:,:)        :: beta_s, T_d, beta_d, chi_map, mask, HI
    real(dp), allocatable, dimension(:)          :: nuz, chi, tump, accept, prob, par
    real(dp)                                     :: nu_ref_s, beta_s_std, beta_s_mu
    real(dp)                                     :: nu_ref_d, beta_d_std, beta_d_mu, T_d_mu
    character(len=80),  dimension(180)           :: header
    character(len=80),  dimension(3)             :: tqu
    character(len=80), allocatable, dimension(:) :: bands
    character(len=5)                             :: iter_str
  
    dust_file  = 'test_data/npipe6v20_353_map_Q_n0008.fits'
    mask_file  = '~/data/masks/mask_fullsky_n0008.fits'
    synch_file = 'test_data/sim_synch_030_n0008.fits'
    tqu(1)     = 'T'
    tqu(2)     = 'Q'
    tqu(3)     = 'U'
    i          = getsize_fits(dust_file, nside=nside, ordering=ordering, nmaps=nmaps)
    npix       = nside2npix(nside) 
    nbands     = 5
    nlheader   = size(header)
    nmaps      = 1

    call getarg(1,arg1)
    direct = arg1

    !----------------------------------------------------------------------------------------------------------
    allocate(dust_temp(0:npix-1,nmaps), synch_temp(0:npix-1,nmaps))
    allocate(maps(0:npix-1,nmaps,nbands), rmss(0:npix-1,nmaps,nbands), model(0:npix-1,nmaps,nbands))
    allocate(res(0:npix-1,nmaps,nbands),chi_map(0:npix-1,nmaps),mask(0:npix-1,1))
    allocate(fg_amp(0:npix-1,nmaps,nbands,3), T_d(0:npix-1,nmaps), beta_d(0:npix-1,nmaps))
    allocate(dust_free(0:npix-1,nmaps,nbands), synch_free(0:npix-1,nmaps,nbands))
    allocate(cmb_map(0:npix-1,nmaps,nbands), dust_map(0:npix-1,nmaps,nbands), synch_map(0:npix-1,nmaps,nbands))
    allocate(nuz(nbands),beta_s(0:npix-1,nmaps),bands(nbands), HI(0:npix-1,nmaps), par(2))

    allocate(map(0:npix-1,nmaps))
    allocate(rms(0:npix-1,nmaps))
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    niter             = 10         ! # of MC-MC iterations
    iterations        = 100        ! # of iterations in the samplers
    output_iter       = 1          ! Output maps every <- # of iterations
    like_iter         = 1000       ! Output likelihood test every <- # of iterations
    nu_ref_s          = 30.0d0     ! Synchrotron reference frequency
    nu_ref_d          = 353.d0     ! Dust reference frequency
    beta_s            = -3.10d0    ! Synchrotron beta initial guess
    beta_s_mu         = -3.10d0    ! \beta_synch Gaussian prior mean
    beta_s_std        = 0.01d0     ! \beta_synch Gaussian prior std
    beta_samp_nside   = 8          ! \beta_synch nside sampling
    j_corr            = 4          ! which band number has been dust corrected
    output_fg         = .true.     ! Option for outputting foregrounds for all bands
    test              = .true.     ! Testing Metropolis-Hasting apparatus
    !----------------------------------------------------------------------------------------------------------

    bands(1)   = ('sim_data_020_')
    bands(2)   = ('sim_data_030_')
    bands(3)   = ('sim_data_045_')
    bands(4)   = ('sim_data_070_')
    bands(5)   = ('sim_data_100_')

    nuz(1)     = 20.0d0
    nuz(2)     = 30.0d0
    nuz(3)     = 45.0d0
    nuz(4)     = 70.0d0
    nuz(5)     = 100.0d0

    loc        = minloc(abs(nuz-nu_ref_s),1)

    !----------------------------------------------------------------------------------------------------------
    ! Read maps

    do j = 1, nbands
        call read_bintab('test_data/' // trim(bands(j)) // 'rms_n0008.fits',rms,npix,nmaps,nullval,anynull,header=header)
        rmss(:,:,j) = rms
        call read_bintab('test_data/' // trim(bands(j)) // 'n0008.fits', &
        map,npix,nmaps,nullval,anynull,header=header)
        maps(:,:,j) = map
    end do

    call read_bintab(dust_file,dust_temp,npix,nmaps,nullval,anynull,header=header)
    call read_bintab(mask_file,mask,npix,1,nullval,anynull,header=header)
    call read_bintab(synch_file,synch_temp,npix,nmaps,nullval,anynull,header=header)

    !----------------------------------------------------------------------------------------------------------
    ! Metropolis-Hastings testing things

    allocate(chi(iterations*niter+1), tump(iterations*niter+1), accept(iterations*niter+1), prob(iterations*niter+1))

    !----------------------------------------------------------------------------------------------------------

    ! fg_amp(:,:,loc,1) = synch_temp

    ! fg_amp(:,1,1,2) = 1.644429d-2
    ! fg_amp(:,1,2,2) = 8.39923d-3
    ! fg_amp(:,1,3,2) = 1.19689d-3
    ! fg_amp(:,1,4,2) = 5.890824d-2
    ! fg_amp(:,1,5,2) = 2.9665593d-1

    ! do j = 1, nbands
    !     dust_map(:,1,j) = fg_amp(:,1,j,2)*dust_temp(:,1)
    ! end do

    !----------------------------------------------------------------------------------------------------------
    ! Calculation portion
    !----------------------------------------------------------------------------------------------------------

    do k = 1, nmaps
        
        ! do i = 0, npix-1
        !     par(1) = beta_s(i,k)
        !     do j = 1, nbands
        !         fg_amp(i,k,j,1) = fg_amp(i,k,loc,1)*compute_spectrum('synch',nuz(j),par)
        !     end do
        ! end do
        ! synch_map(:,k,:)        = fg_amp(:,k,:,1)

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
        
            ! write(*,*) 'Iteration', iter
            ! write(*,*) '-----------------------'
            ! write(*,*) ''

            call sample_joint_amp(npix,k)

            ! lower = 0.0
            ! call cholesky_decomp(mat,lower,3)


            ! write(*,*) 'Sampling A_synch'
            ! ! -------------------------------------------------------------------------------------------------------------------
            ! fg_amp(:,k,loc,1)       = sample_spec_amp((maps(:,k,:)-dust_map(:,k,:)),'synch', rmss(:,k,:))
        
            do i = 0, npix-1
                par(1) = beta_s(i,k)
                do j = 1, nbands
                    fg_amp(i,k,j,1) = fg_amp(i,k,loc,1)*compute_spectrum('synch',nuz(j),par)
                end do
            end do
            synch_map(:,k,:)        = fg_amp(:,k,:,1)

            ! ! write(*,*) 'Sampling dust amplitudes'
            ! ! -------------------------------------------------------------------------------------------------------------------
            do j = 1, nbands
            !     fg_amp(:,k,j,2) = temp_fit((maps(:,k,j)-synch_map(:,k,j)),dust_temp(:,k),rmss(:,k,j),mask(:,1))
                dust_map(:,k,j) = fg_amp(:,k,j,2)*dust_temp(:,k)
            end do

            ! ! -------------------------------------------------------------------------------------------------------------------

            ! ! write(*,*) 'Sampling Beta at nside ', beta_samp_nside
            ! ! -------------------------------------------------------------------------------------------------------------------
            ! call sample_index((maps-dust_map),'synch',beta_samp_nside,k)
            ! ! beta_s(:,k)             = sample_beta((maps-dust_map), npix, k, beta_std, beta_s, beta_samp_nside)
            ! do i = 0, npix-1
            !     par(1) = beta_s(i,k)
            !     do j = 1, nbands
            !         fg_amp(i,k,j,1) = fg_amp(i,k,loc,1)*compute_spectrum('synch',nuz(j),par)
            !     end do
            ! end do
            ! synch_map(:,k,:)        = fg_amp(:,k,:,1)
            ! -------------------------------------------------------------------------------------------------------------------

            res                     = maps - synch_map - dust_map

            if (mod(iter, 10) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, f6.3, a, f7.3, a, f8.4, a, 5e10.3)')&
                 iter, " - chisq: " , compute_chisq(fg_amp,k), " - A_s: ",&
                 fg_amp(350,k,2,1),  " - beta_s: ", sum(beta_s(:,k))/npix, ' - A_d: ', fg_amp(0,k,:,2)
            end if

            call write_data
            if (mod(iter,output_iter) .EQ. 0) then
                call write_maps(k)
            end if

            if (mod(iter,like_iter) .EQ. 0) then
                call likelihood(fg_amp(350,k,2,1),fg_amp(350,k,:,2),beta_s(350,k),k)
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
        norm  = sum(cov(:)*template(:)**2.d0)/sum(t_mask)
        
        sum1 = 0.d0
        sum2 = 0.d0
    
        do i=0,npix-1
            sum1   = sum1 + (data(i)*template(i))*cov(i)*t_mask(i)
            sum2   = sum2 + template(i)**2.d0*cov(i)*t_mask(i)
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
        ! Samples spectral amplitude (per pixel), following the spectrum of foreground 'type'. Returns a full map of amplitudes.
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
                sum1          = sum1 + (data(i,j)*spec)*cov(i,j)
                sum2          = sum2 + (spec)**2.d0*cov(i,j)
                norm(i)       = norm(i) + cov(i,j)*(spec)**2.d0
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
        real(dp)                                               :: mu, sigma

        real(dp)                                               :: naccept   
        logical                                                :: exist

        map2fit = data
        cov     = rmss*rmss

        if (trim(type) == 'synch') then 
            prior(1) = beta_s_mu
            prior(2) = beta_s_std
            indx     = beta_s
        else if (trim(type) == 'dust') then 
            prior(1) = beta_d_mu
            prior(2) = beta_d_std
            indx     = beta_d
        end if

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
        do i = 0, npix2-1
            a       = 0.d0
            sol     = indx_low(i,map_n)
            pars(1) = indx_low(i,map_n)

            ! Chi-square from the most recent Gibbs chain update
            do j = 1, nbands
                a = a + (((fg_amp_low(i,map_n,loc) * compute_spectrum(type,nuz(j),pars)) &
                      - data_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
            end do
            c   = a
            if (test .eqv. .true.) then
                if (i == 100) then
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
                        if (i == 100) then
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
                                if (i == 100) then
                                    naccept = naccept + 1
                                end if
                            end if
                        end if
                    end if
                    ! Metropolis testing apparatus
                    !-----------------------------
                    if (test .eqv. .true.) then
                        if (i == 100) then
                            prob(l+1)   = p
                        end if
                    end if
                end if
                if (test .eqv. .true.) then
                    if (i == 100) then 
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
                if (i == 100) then
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
                    if (i == 100) then

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
                            if (i == 100) then
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
                                    if (i == 100) then
                                        naccept_t_d = naccept_t_d + 1
                                    end if
                                end if
                            end if
                        end if
                        ! Metropolis testing apparatus
                        !-----------------------------
                        if (test .eqv. .true.) then
                            if (i == 100) then
                                prob((iter-1)*100+l+1)   = p
                            end if
                        end if
                    end if
                    if (test .eqv. .true.) then
                        if (i == 100) then 
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


    subroutine sample_joint_amp(npix, map_n)
        !------------------------------------------------------------------------
        ! Solving the matrix equation Ab = c                                    |
        !                                                                       |
        ! For this purpose, we solve the equation:                              |
        !         sum_nu ((T_nu)^T N^-1 T_nu)amp = sum_nu ((T_nu)^T d_nu)       |
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
        real(dp), allocatable, dimension(:,:,:)   :: T_nu, T_nu_T, covar, A_1, A_2
        real(dp), allocatable, dimension(:,:)     :: A, lower, upper, c_1, dats
        real(dp), allocatable, dimension(:)       :: b, c, d
        integer(i4b)                              :: x, y, z

        x = npix
        y = npix+nbands-1
        z = nbands

        allocate(T_nu(x,y,z),T_nu_T(y,x,z),dats(x,z))
        allocate(A_1(y,x,z),A_2(y,y,z)) 
        allocate(A(y,y),b(y),c(y),d(y))
        allocate(lower(y,y),upper(y,y))
        allocate(covar(x,x,z),c_1(y,z))

        do i=1, x
            do j=1,z
                covar(i,i,j) = rmss(i-1,map_n,j)**2
                dats(i,j)    = maps(i-1,map_n,j)
            end do
        end do

        ! Matrix equation initialization
        do i=1,x
            l = 1
            do j=1, z
                T_nu(i,i,j)      = (nuz(j)/nu_ref_s)**beta_s(i-1,map_n)
                if (j /= j_corr) then
                    T_nu(i,x+l,j) = dust_temp(i-1,map_n)
                    l = l + 1
                end if
            end do
        end do


        do j=1, nbands
            T_nu_T(:,:,j) = transpose(T_nu(:,:,j))
            c_1(:,j)      = matmul(T_nu_T(:,:,j),dats(:,j))
            A_1(:,:,j)    = matmul(T_nu_T(:,:,j),covar(:,:,j))
            A_2(:,:,j)    = matmul(A_1(:,:,j),T_nu(:,:,j)) 
            A(:,:)        = A(:,:) + A_2(:,:,j)
            c(:)          = c(:) + c_1(:,j)
        end do

        call cholesky_decomp(A,lower,y)
        upper = transpose(lower)
        call forward_sub(lower,d,c)
        call backward_sub(upper,b,d)

        do i =1, x
            fg_amp(i-1,map_n,loc,1) = b(i)
        end do
        l = 1
        do while (l .lt. (nbands-1))
            do j= 1, z
                if (j /= j_corr) then
                    ! write(*,*) l
                    fg_amp(:,map_n,j,2) = b(x+l)
                    write(*,*) b(x+l)
                    l = l + 1
                end if
            end do
        end do
        ! stop

        deallocate(T_nu,T_nu_T)
        ! deallocate(amp_vec)
        deallocate(A_1)
        deallocate(A_2)
        deallocate(A)
        deallocate(b)
        deallocate(c)
        deallocate(d)
        deallocate(lower,upper)
        deallocate(covar,c_1)

    end subroutine sample_joint_amp
  
    subroutine cholesky_decomp(A,lower,n)
        implicit none
        real(dp), dimension(:,:), intent(in)  :: A
        real(dp), dimension(:,:), intent(out) :: lower
        integer(i4b),             intent(in)  :: n
        integer(i4b)                          :: ip
        real(dp)                              :: s

        do i = 1, n
            do j = 1, i
                s = 0
                if (j == i) then
                    do ip = 1, j-1
                        s = s + (lower(ip,j))**2
                    end do
                    lower(j,j) = sqrt(A(j,j)-s)
                else
                    do ip = 1, j-1
                        s = lower(ip,i)*lower(ip,j)
                    end do
                    lower(j,i) = (A(j,i)-s)/lower(j,j)
                end if
            end do
        end do

    end subroutine cholesky_decomp

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

        title = trim(direct) // 'pixel_350_A_d.dat'
        inquire(file=title,exist=exist)
        if (exist) then
            open(30,file=title, status="old",position="append", action="write")
        else
            open(30,file=title, status="new", action="write")
        endif
        write(30,*) dust_map(350,k,loc)
        close(30)

        title = trim(direct) // 'pixel_350_A_s.dat'
        inquire(file=title,exist=exist)
        if (exist) then
            open(31,file=title, status="old",position="append", action="write")
        else
            open(31,file=title, status="new", action="write")
        endif
        write(31,*) fg_amp(350,k,loc,1)
        close(31)

        title = trim(direct) // 'pixel_350_beta_s.dat'
        inquire(file=title,exist=exist)
        if (exist) then
            open(32,file=title, status="old",position="append", action="write")
        else
            open(32,file=title, status="new", action="write")
        endif
        write(32,*) beta_s(350,k)
        close(32)

        title = trim(direct) // 'total_chisq.dat'
        inquire(file=title,exist=exist)
        if (exist) then
            open(33,file=title, status="old",position="append", action="write")
        else
            open(33,file=title, status="new", action="write")
        endif
        write(33,*) compute_chisq(fg_amp,k)
        close(33)

        inquire(file=trim(direct) // 'dust_' // trim(tqu(k)) // '_amplitudes.dat',exist=exist)
        if (exist) then
            open(34,file = trim(direct) // 'dust_' // trim(tqu(k)) // '_amplitudes.dat', status="old", &
                        position="append", action="write")
        else
            open(34,file = trim(direct) // 'dust_' // trim(tqu(k)) // '_amplitudes.dat', status="new", action="write")
        endif
        write(34,'(12(E17.8))') fg_amp(1,k,:,2)
        close(34)

    end subroutine write_data
  
! Still need to rewrite vv

    function compute_chisq(amp,map_n)
      use healpix_types
      implicit none
      real(dp), dimension(0:npix-1,nmaps,nbands,3), intent(in)   :: amp
      integer(i4b), intent(in)                                   :: map_n
      real(dp)                                                   :: compute_chisq, s, signal, chisq
      integer(i4b)                                               :: m,n
  
      chisq = 0.d0
      do i = 0, npix-1
        do j = 1, nbands
            s = 0.d0
            signal = amp(i,map_n,j,1)*mask(i,1) + amp(i,map_n,j,2)*dust_temp(i,map_n)*mask(i,1)
            s = s + signal
            chisq = chisq + (((maps(i,map_n,j) - s)**2))/(rmss(i,map_n,j)**2)*mask(i,1)
        end do
      end do 
      compute_chisq = chisq/(sum(mask(:,1))+nbands+3) ! n-1 dof, npix + nbands + A_s + A_dust + A_cmb + \beta_s
    end function compute_chisq

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
                signal    = A_s(x)*compute_spectrum('synch',nuz(j),pars) + fg_amp(350,map_n,z,2)*dust_temp(350,map_n)
                A_s_ln(x) = A_s_ln(x) + (((maps(350,map_n,z) - signal)**2.d0)/(rmss(350,map_n,z)**2))
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
                    signal    = fg_amp(y,map_n,loc,map_n)*compute_spectrum('synch',nuz(j),pars) + A_d(x)*dust_temp(y,map_n)
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