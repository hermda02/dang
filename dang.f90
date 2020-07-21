program dang
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
    use init_mod
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
      
    integer(i4b)       :: i, j, k, l, iter, npix, nside, nmaps, ordering, ln
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
    real(dp), allocatable, dimension(:)          :: dust_amps
    real(dp)                                     :: chisq, temp_norm_01, temp_norm_02
    character(len=80), dimension(180)            :: header
    character(len=80), dimension(3)              :: tqu
    character(len=80), allocatable, dimension(:) :: joint_comps
    character(len=10)                            :: solver
    character(len=5)                             :: iter_str
    logical(lgt), allocatable, dimension(:)      :: j_corr01, j_corr02

    real(dp), allocatable, dimension(:,:)        :: mat_test, mat_l, mat_u
    real(dp), allocatable, dimension(:)          :: x, b, d

    ! Object Orient
    type(params)       :: par

    call read_param_file(par)

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
    nbands            = par%numband
    nfgs              = par%ncomp
    nlheader          = size(header)
    nmaps             = 1

    niter             = par%ngibbs       ! # of MC-MC iterations
    iterations        = par%nsample      ! # of iterations in the samplers
    output_iter       = par%iter_out     ! Output maps every <- # of iterations
    output_fg         = par%output_fg    ! Option for outputting foregrounds for all bands
    direct            = par%outdir       ! Output directory name
    beta_samp_nside   = 4          ! \beta_synch nside sampling
    !----------------------------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------------------------
    ! Array allocation
    allocate(template_01(0:npix-1,nmaps), template_02(0:npix-1,nmaps), j_corr01(nbands), j_corr02(nbands))
    allocate(maps(0:npix-1,nmaps,nbands), rmss(0:npix-1,nmaps,nbands), nodust(0:npix-1,nmaps,nbands))
    allocate(mask(0:npix-1,1), res(0:npix-1,nmaps,nbands), chi_map(0:npix-1,nmaps))
    allocate(fg_amp(0:npix-1,nmaps,nbands,nfgs), beta_s(0:npix-1,nmaps))
    allocate(T_d(0:npix-1,nmaps), beta_d(0:npix-1,nmaps))
    allocate(cmb_map(0:npix-1,nmaps,nbands), dust_map(0:npix-1,nmaps,nbands), synch_map(0:npix-1,nmaps,nbands))
    allocate(dust_amps(nbands))
    allocate(map(0:npix-1,nmaps))
    allocate(rms(0:npix-1,nmaps))
    !----------------------------------------------------------------------------------------------------------
    beta_s     = -3.10d0    ! Synchrotron beta initial guess
    beta_d     = 1.60d0     ! Dust beta initial guess

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

    deallocate(map,rms)

    call read_bintab(template_file_01,template_01,npix,nmaps,nullval,anynull,header=header)
    call read_bintab(mask_file,mask,npix,1,nullval,anynull,header=header)
    call read_bintab(template_file_02,template_02,npix,nmaps,nullval,anynull,header=header)

    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Choose which bands for fitting templates
    j_corr01(1) = .true.
    j_corr01(2) = .true.
    j_corr01(3) = .true.
    j_corr01(4) = .true.
    j_corr01(5) = .false.

    j_corr02(1) = .false.
    j_corr02(2) = .false.
    j_corr02(3) = .false.
    j_corr02(4) = .false.
    j_corr02(5) = .false.
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    ! Normalize template to avoid large values in the matrix equation
    temp_norm_01 = maxval(template_01)
    template_01  = template_01/temp_norm_01
    temp_norm_02 = maxval(template_02)
    template_02  = template_02/temp_norm_02
    !----------------------------------------------------------------------------------------------------------
    ! Joint Sampler Info

    allocate(joint_comps(2))

    joint_comps(1) = 'synch'
    joint_comps(2) = 'template01'

    solver = 'lu' ! Method possibilities are 'cg', 'lu', and 'cholesky'
    !----------------------------------------------------------------------------------------------------------
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

        do iter = 1, niter
            call sample_joint_amp(npix,k,trim(solver))  
            dust_amps = fg_amp(0,k,:,2)

            ! -------------------------------------------------------------------------------------------------------------------
            ! Extrapolating A_synch to bands
            if (ANY(joint_comps=='synch')) then
                do i = 0, npix-1
                    do j = 1, nbands
                        fg_amp(i,k,j,1) = fg_amp(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,1,par%dat_nu(j),i,k)
                    end do
                end do
                synch_map(:,k,:)        = fg_amp(:,k,:,1)
            end if

            ! Extrapolating A_dust to bands
            if (ANY(joint_comps=='dust')) then
                do i = 0, npix-1
                    do j = 1, nbands
                        fg_amp(i,k,j,3) = fg_amp(i,k,par%fg_ref_loc(2),3)*compute_spectrum(par,2,par%dat_nu(j),i,k)
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
            call sample_index(par,nodust,beta_samp_nside,1,k)
            do i = 0, npix-1
                do j = 1, nbands
                    fg_amp(i,k,j,1) = fg_amp(i,k,par%fg_ref_loc(1),1)*compute_spectrum(par,1,par%dat_nu(j),i,k)
                end do
            end do
            synch_map(:,k,:)        = fg_amp(:,k,:,1)
            ! -------------------------------------------------------------------------------------------------------------------

            res       = maps - synch_map - dust_map

            if (mod(iter, 1) == 0 .or. iter == 1) then
                write(*,fmt='(i6, a, f10.3, a, f7.3, a, f8.4, a, 6e10.3)')&
                 iter, " - chisq: " , chisq, " - A_s: ",&
                 fg_amp(100,k,par%fg_ref_loc(1),1),  " - beta_s: ",&
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

    function compute_spectrum(self, comp, freq, pix, mapn, param)
        implicit none
        class(params)                  :: self
        real(dp),           intent(in) :: freq
        integer(i4b),       intent(in) :: comp
        integer(i4b),       intent(in) :: pix
        integer(i4b),       intent(in) :: mapn
        real(dp), optional             :: param
        real(dp)                       :: y, rj_cmb ! Conversion factor from K_{RJ} -> K_{CMB}
        real(dp)                       :: z, compute_spectrum

        y = h*(freq*1.0d9) / (k_B*T_CMB)
        rj_cmb = ((exp(y)-1)**2.d0)/(y**2.d0*exp(y))

        if (trim(self%fg_label(comp)) == 'synch') then
           if (present(param)) then
              compute_spectrum = (freq/self%fg_nu_ref(comp))**param
           else 
              compute_spectrum = (freq/self%fg_nu_ref(comp))**beta_s(pix,mapn)
           end if
        else if (trim(self%fg_label(comp)) == 'mbb') then
           z = h / (k_B*T_d(pix,mapn))
           compute_spectrum = (exp(z*self%fg_nu_ref(comp)*1d9)-1.d0) / &
                (exp(z*freq*1d9)-1.d0) * (freq/self%fg_nu_ref(comp))**(beta_d(pix,mapn)+1.d0)!*rj_cmb
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
  
    function sample_spec_amp(self, data, noise, comp, mapn)
        !------------------------------------------------------------------------
        ! Samples spectral amplitude (per pixel), following the spectrum of foreground 'type'. Returns a full map of amplitudes.
        !------------------------------------------------------------------------
        implicit none
  
        class(params)                                    :: self
        real(dp), dimension(0:npix-1,nbands), intent(in) :: data, noise
        integer(i4b),                         intent(in) :: comp
        integer(i4b),                         intent(in) :: mapn
        real(dp)                                         :: sum1, sum2, spec
        real(dp)                                         :: chi, chi_0, chi_00, p
        real(dp)                                         :: amp, num, t, sam
        real(dp), dimension(2)                           :: pars
        real(dp), dimension(nbands)                      :: tmp
        real(dp), dimension(0:npix-1)                    :: norm
        real(dp), dimension(0:npix-1,nbands)             :: cov
        real(dp), dimension(0:npix-1)                    :: sample_spec_amp

        cov = noise*2

        do i = 0, npix-1
            sum1 = 0.0d0
            sum2 = 0.0d0
            do j = 1, nbands
                spec          = compute_spectrum(self,comp,par%dat_nu(j),i,mapn)
                sum1           = sum1 + (data(i,j)*spec)/cov(i,j)
                sum2           = sum2 + (spec)**2.d0/cov(i,j)
                norm(i)        = norm(i) + ((spec)**2.d0)/cov(i,j)
            end do
            norm(i)            = norm(i)!/nbands
            amp                = sum1/sum2
            sample_spec_amp(i) = amp + rand_normal(0.d0,1.d0)/sqrt(norm(i))
        end do

    end function sample_spec_amp
  
    subroutine sample_index(self, data, nside2, comp, map_n)
        implicit none
  
        class(params)                                          :: self
        integer(i4b),                               intent(in) :: map_n, nside2, comp
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
        if (trim(self%fg_label(comp)) == 'synch') then 
            indx     = beta_s
        else if (trim(self%fg_label(comp)) == 'dust') then 
            indx     = beta_d
        end if


        if (mod(iter,output_iter) .EQ. 0) then
           write(*,*) 'Sampling ' // trim(self%fg_label(comp)) // ' beta at nside', beta_samp_nside
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
            a         = 0.d0
            sol       = indx_low(i,map_n)

            ! Chi-square from the most recent Gibbs chain update
            do j = 1, nbands
                a = a + (((fg_amp_low(i,map_n,par%fg_ref_loc(1)) * compute_spectrum(self,comp,par%dat_nu(j),i,map_n,sol)) &
                      - data_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
            end do
            c = a

            do l = 1, iterations

                ! Sampling from the prior
                t         = rand_normal(self%fg_gauss(comp,1,1), self%fg_gauss(comp,1,2))
                b         = 0.d0

                do j = 1, nbands
                    tmp(j) = fg_amp_low(i,map_n,par%fg_ref_loc(1))*compute_spectrum(self,comp,par%dat_nu(j),i,map_n,t)
                    b      = b + ((tmp(j)-data_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
                end do
                b = b

                if (b < c .and. t .lt. self%fg_uni(comp,1,2) .and. t .gt. self%fg_uni(comp,1,1)) then
                    sam = t
                    c   = b
                else
                    x(2) = exp(0.5d0*(c-b))
                    p = minval(x)
                    call RANDOM_NUMBER(num)
                    if (num < p) then
                        if (t .lt. self%fg_uni(comp,1,2) .and. t .gt. self%fg_uni(comp,1,1)) then
                            sam = t
                            c   = b
                        end if
                    end if
                end if
            end do
            sol = sam
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

        if (trim(self%fg_label(comp)) == 'synch') then 
            beta_s(:,k) = indx_sample
        else if (trim(self%fg_label(comp)) == 'dust') then 
            beta_d(:,k) = indx_sample
        end if
        deallocate(data_low)
        deallocate(fg_amp_low)
        deallocate(indx_low)
        deallocate(rms_low)
        deallocate(indx_sample_low)

    end subroutine sample_index


    ! This architecture of this function has not been modified yet
    function sample_HI_T(self, band, npix, map_n, sigma, T_map, nside2)
        implicit none
  
        class(params)                                          :: self
        integer(i4b), intent(in)                               :: npix, map_n, nside2
        integer(i4b)                                           :: nside1, npix2
        real(dp), intent(in)                                   :: sigma
        real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: band
        real(dp), dimension(0:npix-1,nmaps), intent(in)        :: T_map
        real(dp), dimension(0:npix-1)                          :: sample_hi_T
        real(dp), dimension(0:npix-1,nmaps,nbands)             :: map2fit, cov
        real(dp), dimension(0:npix-1,nmaps)                    :: te
        real(dp), allocatable, dimension(:,:,:)                :: band_low, fg_amp_low, rms_low
        real(dp), allocatable, dimension(:,:)                  :: T_low
        real(dp), allocatable, dimension(:)                    :: sample_T_low
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
                    a = a + (((fg_amp_low(i,map_n,j) * HI(i,1)*planck(par%dat_nu(j)*1.d9,sol)) &
                        - band_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2
                end do
                c   = a

                do l = 1, iterations
                    ! Begin sampling from the prior
                    t = rand_normal(self%fg_gauss(2,2,1), self%fg_gauss(2,2,2))
                    b = 0.d0
                    do j = 1, nbands
                        tmp(j) = fg_amp_low(i,map_n,j)*HI(i,1)*planck(par%dat_nu(j)*1.d9,t)
                        b      = b + ((tmp(j)-band_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
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
                call udgrade_ring(sample_T_low,nside2, sample_hi_T,nside1)
                call convert_nest2ring(nside2, sample_T_low)
            else
                call udgrade_nest(sample_T_low,nside2, sample_hi_T,nside1)
            end if
        else
             sample_HI_T =  sample_T_low
        end if

        deallocate(band_low)
        deallocate(fg_amp_low)
        deallocate(T_low)
        deallocate(rms_low)
        deallocate(sample_T_low)

    end function sample_HI_T
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
                    do j = 1, z
                        T_nu(i,w+i,j) = compute_spectrum(par,1,par%dat_nu(j),i-1,map_n)
                    end do
                end do
                w = w + x
            else if (joint_comps(m) == 'dust') then
                do i = 1, x
                    do j = 1, z
                        T_nu(i,w+i,j) = compute_spectrum(par,2,par%dat_nu(j),i-1,map_n)
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
                            T_nu(i,w+l,j) = template_02(i-1,map_n)
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
        b     = b  + samp

        ! Output amplitudes to the appropriate variables
        do m = 1, size(joint_comps)
            if (joint_comps(m) == 'synch') then
                do i = 1, x
                    fg_amp(i-1,map_n,par%fg_ref_loc(1),1) = b(i)
                end do
                b = b(x+1:)
            else if (joint_comps(m) == 'dust') then
                do i = 1, x
                    fg_amp(i-1,map_n,par%fg_ref_loc(2),3) = b(i)
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
                        else
                           fg_amp(:,map_n,j,2) = 0.d0
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

    subroutine write_maps(nm)
        implicit none

        integer(i4b), intent(in)    :: nm
        real(dp), dimension(npix,1) :: map

        write(iter_str, '(i0.5)') iter
        if (output_fg .eqv. .true.) then
            do j = 1, nbands
                title = trim(direct) // trim(par%dat_label(j)) // '_dust_fit_'// trim(tqu(nm)) & 
                        // '_' // trim(iter_str) // '.fits'
                map(:,1)   = dust_map(:,nm,j)
                call write_bintab(map,npix,1, header, nlheader, trim(title))
                title = trim(direct) // trim(par%dat_label(j)) // '_synch_amplitude_' //  trim(tqu(nm)) &
                        // '_' // trim(iter_str) // '.fits'
                map(:,1)   = fg_amp(:,nm,j,1)
                call write_bintab(map,npix,1, header, nlheader, trim(title))
            end do
        else 
            title = trim(direct) // trim(par%dat_label(par%fg_ref_loc(1))) // '_synch_amplitude_' //  trim(tqu(nm)) &
                    // '_' // trim(iter_str) // '.fits'
            map(:,1)   = fg_amp(:,nm,par%fg_ref_loc(1),1)
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
        write(30,*) dust_map(100,k,par%fg_ref_loc(1))
        close(30)

        title = trim(direct) // 'pixel_100_A_s.dat'
        inquire(file=title,exist=exist)
        if (exist) then
            open(31,file=title, status="old",position="append", action="write")
        else
            open(31,file=title, status="new", action="write")
        endif
        write(31,*) fg_amp(100,k,par%fg_ref_loc(1),1)
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

  end program dang
