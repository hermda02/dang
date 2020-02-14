program ame_fit
    use healpix_types
    use pix_tools
    use fitstools
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
    ! side the WMAP data, we iteratively fit CMB (global), synchrotron (per pixel), and dust (global).    |
    ! Once this has been done, estimates to the level of AME polarization will be made using the lowest   |
    ! frequency bands used here.                                                                          |
    !-----------------------------------------------------------------------------------------------------|  
    
    ! Constants
    real(dp)           :: k_b     = 1.3806503d-23
    real(dp)           :: h       = 1.0545726691251021d-34 * 2.d0 * pi
    real(dp)           :: c       = 2.99792458d8
    real(dp)           :: T_CMB   = 2.7255d0
  
    integer(i4b)       :: i, j, k, l, iter, npix, nside, nmaps, ordering, nlheader, niter, nbands, iterations
    real(dp)           :: nullval
    real(dp)           :: missval = -1.6375d30
    logical(lgt)       :: anynull, double_precision, test
  
    character(len=128) :: dust_map, synch_map, cmb_map
    character(len=128) :: mapfile, title
  
    real(dp), allocatable, dimension(:,:,:,:)    :: fg_amp, temp_amp
    real(dp), allocatable, dimension(:,:,:)      :: maps, rmss, model, res, ame_map
    real(dp), allocatable, dimension(:,:,:)      :: ame_free, synch, synch_free
    real(dp), allocatable, dimension(:,:)        :: dust_temp, synch_temp, cmb_temp, map, rms
    real(dp), allocatable, dimension(:,:)        :: A_ame, A_cmb, A_synch, beta_s, chi_map
    real(dp), allocatable, dimension(:)          :: frequencies
    real(dp)                                     :: nu_ref, beta_std, ki, ki_0
    character(len=80),  dimension(180)           :: header
    character(len=80),  dimension(3)             :: tqu
    character(len=80), allocatable, dimension(:) :: bands
    character(len=4)                             :: iter_str
  
    dust_map  = 'test_data/npipe6v20_353_map_QUADCOR_ZODICOR_DIPCOR_n0064_uK.fits'
    synch_map = 'test_data/synch_amp_init_bp_v2.00_n0064.fits'
    tqu(1)    = 'T'
    tqu(2)    = 'Q'
    tqu(3)    = 'U'
  
    i         = getsize_fits(dust_map, nside=nside, ordering=ordering, nmaps=nmaps)
    npix      = nside2npix(nside) 
    nbands    = 7
  
    !----------------------------------------------------------------------------------------------------------
    allocate(dust_temp(0:npix-1,nmaps), synch_temp(0:npix-1,nmaps), cmb_temp(0:npix-1,nmaps))
    allocate(maps(0:npix-1,nmaps,nbands), rmss(0:npix-1,nmaps,nbands), model(0:npix-1,nmaps,nbands))
    allocate(res(0:npix-1,nmaps,nbands), ame_map(0:npix-1,nmaps,nbands),chi_map(0:npix-1,nmaps))
    allocate(fg_amp(0:npix-1,nmaps,nbands,2),temp_amp(0:npix-1,nmaps,nbands,2))
    allocate(ame_free(0:npix-1,nmaps,nbands), synch_free(0:npix-1,nmaps,nbands), synch(0:npix-1,nmaps,nbands))
    allocate(A_ame(nbands,nmaps), A_cmb(nbands,nmaps), A_synch(0:npix-1,nmaps))
    allocate(frequencies(nbands),beta_s(0:npix-1,nmaps),bands(nbands))
  
    allocate(map(0:npix-1,nmaps))
    allocate(rms(0:npix-1,nmaps))
      !----------------------------------------------------------------------------------------------------------

    niter          = 5       ! # of MC-MC iterations
    iterations     = 100     ! # of iterations in the \beta sampler
    nu_ref         = 44.1d0  ! Synchrotron reference frequency
    beta_s         = -3.1d0  ! \beta_synch Gaussian prior mean
    beta_std       = 1.0d0   ! \beta_synch Gaussian prior std
    test           = .false. ! beta sampler testing

    bands(1)       = ('bp_030_')
    bands(2)       = ('bp_044_')
    bands(3)       = ('bp_070_')
    bands(4)       = ('wmap9_K_')
    bands(5)       = ('wmap9_Ka_')
    bands(6)       = ('wmap9_Q1_')
    bands(7)       = ('wmap9_V1_')

    frequencies(1) = 28.4d0
    frequencies(2) = 44.1d0
    frequencies(3) = 70.3d0
    frequencies(4) = 22.8d0
    frequencies(5) = 33.d0
    frequencies(6) = 40.6d0
    frequencies(7) = 60.8d0
  
    !----------------------------------------------------------------------------------------------------------
    ! Read maps

    do j = 1, nbands
        call read_bintab('test_data/' // trim(bands(j)) // 'rms_n0064.fits',rms,npix,nmaps,nullval,anynull,header=header)
        rmss(:,:,j) = rms
    end do


    do j = 1, nbands
        call read_bintab('test_data/' // trim(bands(j)) // 'map_n0064_60arcmin_nocmbff.fits', &
        map,npix,nmaps,nullval,anynull,header=header)
        maps(:,:,j) = map
    end do

    call read_bintab(synch_map,synch_temp,npix,nmaps,nullval,anynull,header=header)
    call read_bintab(dust_map,dust_temp,npix,nmaps,nullval,anynull,header=header)

    maps(:,1,4) = maps(:,1,4) + 15.48542
    maps(:,1,5) = maps(:,1,5) + 1.222034
    maps(:,1,6) = maps(:,1,6) - 0.8288112
    maps(:,1,7) = maps(:,1,7) - 2.336047

    dust_temp(:,1)  = dust_temp(:,1) + 1175.8888823383916
    synch_temp(:,1) = synch_temp(:,1)*1.4434668d-6
    synch_temp(:,2) = synch_temp(:,2)*0.2995193844767811
    synch_temp(:,3) = synch_temp(:,3)*0.2995193844767811

    !----------------------------------------------------------------------------------------------------------
  
    nlheader = size(header)
  
    ! fg_amp(:,:,2,1) = synch_temp
    fg_amp(:,:,:,2) = 0.d0

    do k = 1, nmaps    
        do i = 0, npix-1
            do j = 1, nbands
            fg_amp(i,k,j,1) = fg_amp(i,k,2,1)*compute_fg(1,frequencies(j),beta_s(i,k))
            end do
        end do
    end do
    do k = 1, nmaps   
        write(*,*) 'Mean A_synch: ', sum(fg_amp(:,k,2,1))/npix
    end do
    
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------

    do k = 1, nmaps
        if (k == 1) then
            write(*,*) 'Sampling Temperature'
            write(*,*) '-----------------------'
        else if (k == 2) then
            write(*,*) 'Stokes Q'
            write(*,*) '-----------------------'
        else if (k == 3) then
            write(*,*) 'Stokes U'
            write(*,*) '-----------------------'
        end if    

        do iter = 1, niter
        
            write(*,*) 'Iteration', iter
            write(*,*) '-----------------------'
            write(*,*) ''
            
            write(*,*) 'Sampling AME amplitudes'

            synch(:,k,:)      = fg_amp(:,k,:,1)
            synch_free(:,k,:) = maps(:,k,:) - synch(:,k,:)
        
            do j = 1, nbands
                fg_amp(:,k,j,2)       = temp_fit(synch_free(:,:,j),dust_temp(:,:),rmss(:,:,j),k)
                ame_map(:,k,j)        = fg_amp(:,k,j,2)*dust_temp(:,k)
                ame_free(:,k,j)       = maps(:,k,j) - ame_map(:,k,j)
            end do
        
            write(*,*) 'Chi^2 ', compute_chisq(fg_amp,k)
            write(*,*) ''

            write(*,*) 'Sampling A_synch'
        
            temp_amp(:,k,2,1) = sample_s_amp(ame_free,fg_amp(:,:,3,1), rmss(:,:,:), k)
        
            do i = 0, npix-1
                do j = 1, nbands
                    temp_amp(i,k,j,1) = temp_amp(i,k,2,1)*compute_fg(1,frequencies(j),beta_s(i,k))
                end do
            end do

            fg_amp(:,k,:,1) = temp_amp(:,k,:,1)


            write(*,*) 'Chi^2 ', compute_chisq(fg_amp,k)
            write(*,*) ''

            write(*,*) 'Sampling Beta'

            beta_s(:,k)  = sample_beta(ame_free, npix, k, beta_std, beta_s)
        
            ! write(*,*) 'Mean synch beta: '
            ! write(*,*) sum(beta_s(:,k))/npix
            write(*,*) 'Chi^2 ', compute_chisq(fg_amp,k)

            do i = 0, npix-1
                do j = 1, nbands
                    fg_amp(i,k,j,1) = fg_amp(i,k,2,1)*compute_fg(1,frequencies(j),beta_s(i,k))
                end do
            end do

            synch(:,k,:)      = fg_amp(:,k,:,1)

            do j = 1, nbands
                ame_map(:,k,j)    = fg_amp(:,k,j,2)*dust_temp(:,k)
            end do
            model = synch + ame_map
            
            res = maps - model
        
            write(*,*) ''

            if (mod(iter,1) .EQ. 0) then

                write(iter_str, '(i0.4)') iter
                
                do j = 1, nbands
                    title =  trim(bands(j)) // '_ame_fit_'// trim(tqu(k)) // '_' // trim(iter_str) // '.fits'
                    call write_maps(ame_map(:,k,j),npix,1,title)
                    title = trim(bands(j)) // '_synch_amplitude_' //  trim(tqu(k)) // '_' // trim(iter_str) // '.fits'
                    call write_maps(fg_amp(:,k,j,1),npix,1,title)
                    title = trim(bands(j)) // 'residual_' // trim(tqu(k)) // '_' // trim(iter_str) // '.fits'
                    call write_maps(res(:,k,j),npix,1,title)
                end do
                title = 'synch_beta_' // trim(tqu(k)) // '_' // trim(iter_str) // '.fits'
                call write_maps(beta_s(:,k),npix,1,title)
                do i = 0, npix-1
                    do j = 1, nbands
                        chi_map(i,k) = chi_map(i,k) + (maps(i,k,j) - (synch(i,k,j)+ame_map(i,k,j)))**2.d0/rmss(i,k,j)**2.d0
                    end do
                end do
                chi_map(:,k) = chi_map(:,k)/nbands
                title = 'chisq_' // trim(tqu(k)) // '_' // trim(iter_str) // '.fits'
                call write_maps(chi_map(:,k),npix,1,title)

            end if
        end do    
    end do
  
  contains
  
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
  
    function temp_fit(band,temp,rms,map_n)
        implicit none
    
        integer(i4b), intent(in)            :: map_n
        real(dp), dimension(0:npix-1,nmaps) :: band, temp, rms, cov
        real(dp)                            :: temp_fit, norm, chiz
        real(dp)                            :: amplitude, amp, new
        real(dp)                            :: sum1, sum2, chi_0
    
        cov = rms*2
    
        ! Uncertainty:
        norm  = sum(cov(:,map_n)*temp(:,map_n)**2.d0)/npix
        
        sum1 = 0.d0
        sum2 = 0.d0
    
        do i=0,npix-1
            sum1   = sum1 + (band(i,map_n)*temp(i,map_n))*cov(i,map_n)
            sum2   = sum2 + temp(i,map_n)**2.d0*cov(i,map_n)
        end do

        amp   = sum1/sum2

        chi_0 = 0.d0
        do i = 0, npix-1
            chi_0 = chi_0 + (band(i,map_n)-amp*temp(i,map_n))**2.d0/cov(i,map_n)
        end do

        do l = 1, iterations
            chiz = 0.d0
            new = amp + rand_normal(0.d0,1.d0)/sqrt(norm)
            do i = 0, npix-1
                chiz = chiz + (band(i,map_n)-new*temp(i,map_n))**2.d0/cov(i,map_n)
            end do
            if (chiz < chi_0) then
                amp = new
                chi_0 = chiz
            end if
        end do
        temp_fit  = amp
  
    end function temp_fit
  
    function sample_s_amp(band,A,rms,map_n)
      implicit none
  
      integer(i4b),                               intent(in) :: map_n
      real(dp), dimension(0:npix-1,nmaps),        intent(in) :: A
      real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: band, rms
      real(dp)                                               :: sum1, sum2, powerl, chi_0, chi_00, ch
      real(dp)                                               :: amp, new, chiz, damp
      real(dp), dimension(0:npix-1,nmaps)                    :: norm
      real(dp), dimension(0:npix-1,nmaps,nbands)             :: cov
      real(dp), dimension(0:npix-1,nmaps,nbands,2)           :: trmp
      real(dp), dimension(0:npix-1)                          :: sample_s_amp
      integer(i4b)                                           :: fail

      cov = rms*2

      do i = 0, npix-1
        sum1 = 0.0d0
        sum2 = 0.0d0
        do j = 1, nbands
            powerl = compute_fg(1,frequencies(j),beta_s(i,map_n))
            sum1   = sum1 + (band(i,map_n,j)*powerl)*cov(i,map_n,j)
            sum2   = sum2 + (powerl)**2.d0*cov(i,map_n,j)
            norm(i,map_n) = norm(i,map_n) + cov(i,map_n,j)*(powerl)**2.d0
        end do
        norm(i,map_n) = norm(i,map_n)/nbands
        amp   = sum1/sum2

        chi_00 = 0.d0
        do j = 1, nbands
            chi_00 = chi_00 + (band(i,map_n,j)-(A(i,map_n))*compute_fg(1,frequencies(j),beta_s(i,map_n)))**2.d0/cov(i,map_n,j)
        end do
        chi_0 = 0.d0
        do j = 1, nbands
            chi_0 = chi_0 + (band(i,map_n,j)-(sum1/sum2)*compute_fg(1,frequencies(j),beta_s(i,map_n)))**2.d0/cov(i,map_n,j)
        end do
        do l = 1, iterations
            chiz = 0.d0
            new = amp + rand_normal(0.d0,1.d0)/sqrt(norm(i,map_n))
            do j = 1, nbands
                chiz = chiz + (band(i,map_n,j)-(new*compute_fg(1,frequencies(j),beta_s(i,map_n))))**2.d0/cov(i,map_n,j)
            end do
            if (chiz < chi_0) then
                amp = new
                chi_0 = chiz
            end if
        end do
        if (isnan(amp)) then
            write(*,*) i
            stop
        end if
        if (chi_0 < chi_00) then
            sample_s_amp(i) = amp
        else
            sample_s_amp(i) = A(i,map_n)
        end if
      end do
  
    end function sample_s_amp
  
    function sample_beta(map2fit, npix, map_n, sigma, beta)
        implicit none
  
        integer(i4b), intent(in)                               :: npix, map_n
        real(dp), intent(in)                                   :: sigma
        real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: map2fit
        real(dp), dimension(0:npix-1,nmaps), intent(in)        :: beta
        real(dp), dimension(0:npix-1)                          :: sample_beta
        real(dp), dimension(nbands)                            :: signal, tmp
        real(dp), dimension(2)                                 :: x
        real(dp)                                               :: a, b, c, num, temp, r, p
    
        real(dp), dimension(iterations+1)                    :: chi, tump, accept, prob
        real(dp)                                             :: naccept    
        logical                                              :: exist
    
        x(1) = 1.d0
        do i = 0, npix-1
            a = 0.d0
            temp = beta(i,map_n)
            do j = 1, nbands
              a = a + (((fg_amp(i,map_n,j,1) * compute_fg(1,frequencies(j),temp)) - map2fit(i,map_n,j))**2.d0)/rmss(i,map_n,j)**2
            end do
            c    = a/5
  
            if (test .eqv. .true.) then
                prob(1)   = 1.d0
                chi(1)    = a
                tump(1)   = temp
                accept(1) = 0.d0
                naccept   = 0
            end if
  
            do l = 1, iterations
                r = rand_normal(temp, sigma)
                b = 0.d0
                do j = 1, nbands
                    tmp(j) = fg_amp(i,map_n,j,1)*compute_fg(1,frequencies(j),r)
                    b      = b + ((tmp(j)-map2fit(i,map_n,j))**2.d0)/rmss(i,map_n,j)**2.d0
                end do
                b = b/5
                x(2) = exp(b-c)
                p = minval(x)
                call RANDOM_NUMBER(num)
                if (num > p) then
                    if (r .lt. -2.5 .and. r .gt. -3.5) then
                        temp = r
                        c    = b
                
                        if (test .eqv. .true.) then
                            naccept = naccept + 1
                
                        end if
                    end if
                end if

                ! Metropolis testing apparatus
                !-----------------------------
                if (test .eqv. .true.) then
                    if (i == 1234) then
                        prob(l+1)   = p 
                        chi(l+1)    = c
                        tump(l+1)   = temp
                        accept(l+1) = naccept/l
                    end if
                end if
                !-----------------------------
                
            end do

            ! Metropolis testing apparatus
            !-----------------------------
            if (test .eqv. .true.) then
                if (i == 1234) then
                    inquire(file='prob.dat',exist=exist)
                    if (exist) then
                        open(40,file = 'prob.dat', status="old", position="append", action="write")
                    else
                        open(40,file = 'prob.dat', status="new", action="write")
                    endif
    
                    inquire(file='chi.dat',exist=exist)
                    if (exist) then
                        open(41,file = 'chi.dat', status="old", position="append", action="write")
                    else
                        open(41,file = 'chi.dat', status="new", action="write")
                    endif
    
                    inquire(file='temps.dat',exist=exist)
                    if (exist) then
                        open(42,file = 'temps.dat', status="old", position="append", action="write")
                    else
                        open(42,file = 'temps.dat', status="new", action="write")
                    endif
    
                    inquire(file='accept.dat',exist=exist)
                    if (exist) then
                        open(43,file = 'accept.dat', status="old", position="append", action="write")
                    else
                        open(43,file = 'accept.dat', status="new", action="write")
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
            sample_beta(i) = temp
        end do
    end function sample_beta
  
    subroutine write_maps(map,npix,nmap,title)
      implicit none
  
      integer(i4b), intent(in)       :: npix, nmap
      character(len=80), intent(in) :: title
      real(dp), dimension(npix,nmap) :: map
      call write_bintab(map, npix, nmap, header, nlheader, trim(title))

    end subroutine write_maps
  
    function compute_chisq(amp,map_n)
      use healpix_types
      implicit none
      real(dp), dimension(0:npix-1,nmaps,nbands,2), intent(in)   :: amp
      integer(i4b), intent(in)                                   :: map_n
      real(dp)                                                   :: compute_chisq, s, signal, chisq
      integer(i4b)                                               :: m,n
  
      chisq = 0.d0
      do m = 0, npix-1
        do n = 1, nbands
            s = 0.d0
            signal = amp(m,map_n,n,1)*(compute_fg(1,frequencies(n),beta_s(m,map_n))) &
                   + amp(m,map_n,n,2)*dust_temp(m,map_n)
            s = s + signal
            chisq = chisq + ((maps(m,map_n,n) - s)**2)/(rmss(m,map_n,n)**2)
        end do
      end do 
      compute_chisq = chisq/(npix+nbands-1)
  
    end function compute_chisq
  
    function compute_fg(comp,freq,param)
      implicit none
      integer(i4b), intent(in)           :: comp
      real(dp), intent(in)               :: freq, param
      real(dp)                           :: compute_fg
      
      if (comp == 1) then
          compute_fg  = (freq/nu_ref)**param
      else if (comp == 2) then
          compute_fg = param
      end if
  
    end function compute_fg
    
  end program ame_fit