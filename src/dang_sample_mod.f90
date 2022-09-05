module dang_sample_mod
  use healpix_types
  use pix_tools
  use fitstools
  use udgrade_nr
  use dang_util_mod
  use dang_param_mod
  use dang_linalg_mod
  use dang_component_mod
  use dang_data_mod
  use dang_bp_mod
  use dang_cg_mod
  implicit none
  
  private :: i, j, k, l
  integer(i4b) :: i, j, k, l
  
contains

  subroutine sample_spectral_parameters(ddata)
    ! ============================================ |
    ! This function simply loops over each of the  |
    ! parameters, checking to see if their indices |
    ! ought to be sampled.                         |
    ! ============================================ |

    implicit none

    type(dang_data),             intent(in) :: ddata
    type(dang_comps),            pointer    :: c
    integer(i4b)                            :: i, j, k

    ! Loop over all foregrounds and sample for indices
    do i = 1, ncomp
       c => component_list(i)%p
       if (c%nindices == 0) cycle
       do j = 1, c%nindices
          if (c%sample_index(j)) then
             do k = 1, c%nflag(j)
                if (iand(c%pol_flag(j,k),1) .ne. 0) then
                   write(*,*) 'Sampling spectral index for ', trim(c%label), ', poltype = I.'
                   call sample_index_mh(ddata,c,j,1)
                else if (iand(c%pol_flag(j,k),2) .ne. 0) then
                   write(*,*) 'Sampling spectral index for ', trim(c%label), ', poltype = Q.'
                   call sample_index_mh(ddata,c,j,2)
                else if (iand(c%pol_flag(j,k),4) .ne. 0) then
                   write(*,*) 'Sampling spectral index for ', trim(c%label), ', poltype = U.'
                   call sample_index_mh(ddata,c,j,3)
                else if (iand(c%pol_flag(j,k),8) .ne. 0) then
                   write(*,*) 'Sampling spectral index for ', trim(c%label), ', poltype = Q+U.'
                   call sample_index_mh(ddata,c,j,-1)
                else if (iand(c%pol_flag(j,k),15) .ne. 0) then
                   write(*,*) 'Sampling spectral index for ', trim(c%label), ', poltype = I+Q+U.'
                   call sample_index_mh(ddata,c,j,-2)
                else
                   write(*,*) "There is something wrong with the poltype flag"
                   write(*,*) "for component ", trim(c%label)
                end if
             end do
          end if
       end do
    end do

  end subroutine sample_spectral_parameters
  
  subroutine sample_index_mh(ddata,c,nind,map_n)
    ! ============================================== |
    ! Implementation of the Metropolis-Hastings      |
    ! sampling algorithm. The user is responsible    |
    ! for determining the step size for the sampler. |
    !                                                |
    ! This routine is set up to sample fullsky, or   |
    ! per-pixel, individually for each poltype, or a |
    ! combination of poltypes (Q+U,T+Q+U).           |
    !------------------------------------------------!
    !                                                |
    ! Inputs:                                        |
    !   ddata: type(dang_data)                       |
    !   c: type(dang_comps) - pointer to component   |
    !   nind: integer - index to sample              |
    !   map_n: integer - map to sample, or flag for  |
    !                    poltype combinations        |
    !                                                |
    ! ============================================== |
    implicit none
    
    type(dang_data),             intent(in) :: ddata
    type(dang_comps),   pointer, intent(in) :: c
    integer(i4b),                intent(in) :: nind
    integer(i4b),                intent(in) :: map_n
    type(dang_comps),   pointer             :: c2
    
    real(dp), allocatable, dimension(:,:,:) :: data, model
    real(dp), allocatable, dimension(:,:)   :: index
    real(dp), allocatable, dimension(:)     :: sample, theta

    integer(i4b)                            :: i, j, k
    integer(i4b)                            :: l, m, n, mh_mode
    real(dp)                                :: lnl, lnl_old, lnl_new
    real(dp)                                :: diff, ratio, num

    ! Here is an object we'll call data, which we will correct to be the
    ! portion of the sky signal we wish to fit to
    allocate(data(0:npix-1,nmaps,nbands))
    allocate(model(0:npix-1,nmaps,nbands))
    model(:,:,:) = 0.d0
    data(:,:,:)  = ddata%sig_map

    ! Loop over components, remove all others
    do l = 1, ncomp
       c2 => component_list(l)%p
       if (c2%label /= c%label) then
          if (c2%type /='template') then
             do i = 0, npix-1
                do k = 1, nmaps
                   do j = 1, nbands
                      data(i,k,j) = data(i,k,j) - c2%amplitude(i,k)*c2%eval_sed(j,i,k)
                   end do
                end do
             end do
          else if (c2%type == 'template') then
             do i = 0, npix-1
                do k = 1, nmaps
                   do j = 1, nbands
                      data(i,k,j) = data(i,k,j) - c2%template(i,k)*c2%template_amplitudes(j,k)
                   end do
                end do
             end do
          end if
       end if
    end do

    ! Initialize index map from previous Gibbs iteration
    allocate(sample(c%nindices),theta(c%nindices))
    
    ! Little extra section here for mode determination
    !=================================================
    ! mh_mode :
    !    1 --- corresponds to sampling for the poltype which is input to the routine
    !    2 --- sampling for Q+U jointly
    !    3 --- sampling for I+Q+U jointly
    !=================================================
    if (map_n == -1) then
       mh_mode = 2
    else if (map_n == -2) then
       mh_mode = 3
    else
       mh_mode = 1
    end if
    
    ! Now time to begin sampling
    ! Index mode 1 corresponds to full sky value for the spectral parameter
    if (c%index_mode(nind) == 1) then
       write(*,*) 'Sampling fullsky'
       lnl = 0.d0
       ! Ensure proper handling of poltypes
       if (mh_mode == 1) then
          do l = 1, c%nindices
             sample(l)   = c%indices(0,map_n,l)
          end do
       else if (mh_mode == 2) then
          do l = 1, c%nindices
             sample(l)   = c%indices(0,2,l)
          end do
       else if (mh_mode == 3) then
          do l = 1, c%nindices
             sample(l)   = c%indices(0,1,l)
          end do
       end if
       
       ! Define the model to toss into the likelihood evaluation
       ! Ensure proper handling of poltypes
       if (mh_mode == 1) then
          do i = 0, npix-1
             do j = 1, nbands
                model(i,map_n,j) = c%eval_signal(j,i,map_n,sample)
             end do
          end do
       else if (mh_mode == 2) then
          do i = 0, npix-1
             do k = 2, 3
                do j = 1, nbands
                   model(i,k,j) = c%eval_signal(j,i,k,sample)
                end do
             end do
          end do
       else if (mh_mode == 3) then
          do i = 0, npix-1
             do k = 1, 3
                do j = 1, nbands
                   model(i,k,j) = c%eval_signal(j,i,k,sample)
                end do
             end do
          end do
       end if
       
       ! Evaluate the lnL (already includes the -0.5 out front)
       lnl = evaluate_lnL(data,ddata%rms_map,model,map_n,-1,ddata%masks(:,1))
       
       if (trim(c%prior_type(nind)) == 'gaussian') then
          lnl_old = lnl + log(eval_normal_prior(sample(nind),c%gauss_prior(nind,1),c%gauss_prior(nind,2)))
       else if (trim(c%prior_type(nind)) == 'uniform') then
          lnl_old = lnl
       end if

       ! Now we do the real sampling
       do l = 1, nsample
          
          ! Update theta with the new sample
          theta(nind) = sample(nind) + rand_normal(0.d0,c%gauss_prior(nind,2))

          ! If the prior is outside the allow prior range, move on
          if (theta(nind) > c%uni_prior(nind,2) .or. theta(nind) < c%uni_prior(nind,1)) cycle
          
          ! Evaluate model for likelihood evaluation
          ! Ensure proper handling of poltypes
          if (mh_mode == 1) then
             do i = 0, npix-1
                do j = 1, nbands
                   model(i,map_n,j) = c%eval_signal(j,i,map_n,theta)
                end do
             end do
          else if (mh_mode == 2) then
             do i = 0, npix-1
                do k = 2, 3
                   do j = 1, nbands
                      model(i,k,j) = c%eval_signal(j,i,k,theta)
                   end do
                end do
             end do
          else if (mh_mode == 3) then
             do i = 0, npix-1
                do k = 1, 3
                   do j = 1, nbands
                      model(i,k,j) = c%eval_signal(j,i,k,theta)
                   end do
                end do
             end do
          end if
          
          ! Evaluate the lnL (already includes the -0.5 out front)
          lnl = evaluate_lnL(data,ddata%rms_map,model,map_n,-1,ddata%masks(:,1))
          
          ! Accept/reject
          if (trim(c%prior_type(nind)) == 'gaussian') then
             lnl_new = lnl + log(eval_normal_prior(theta(nind),c%gauss_prior(nind,1),c%gauss_prior(nind,2)))
          else if (trim(c%prior_type(nind)) == 'uniform') then
             lnl_new = lnl
          end if
          
          ! Compute the difference, and get the likelihood ratios
          diff  = lnl_new - lnl_old
          ratio = exp(diff)
          
          ! Accept/reject step
          if (trim(ml_mode) == 'optimize') then
             if (ratio > 1.d0) then
                sample(nind) = theta(nind)
                lnl_old      = lnl_new
             end if
          else if (trim(ml_mode) == 'sample') then
             call RANDOM_NUMBER(num)
             if (ratio > num) then
                sample(nind) = theta(nind)
                lnl_old      = lnl_new
             end if
          end if
       end do
       ! Ensure proper handling of poltypes
       if (mh_mode == 1) then
          c%indices(:,map_n,nind) = sample(nind)
       else if (mh_mode == 2) then
          c%indices(:,2:3,nind) = sample(nind)
       else if (mh_mode == 3) then
          c%indices(:,1:3,nind) = sample(nind)
       end if

    else if (c%index_mode(nind) == 2) then
       write(*,*) 'Sampling per-pixel'
       ! Pixel-by-pixel

       !$OMP PARALLEL PRIVATE(i,j,k,l,lnl,sample,theta,lnl_old,lnl_new,diff,ratio)
       !$OMP DO SCHEDULE(static)
       do i = 0, npix-1
          if (ddata%masks(i,1) == 0.d0 .or. ddata%masks(i,1) == 0.d0) cycle

          ! Initialize the MH chain
          lnl      = 0.d0
          ! Ensure proper handling of poltypes
          if (mh_mode == 1) then
             do l = 1, c%nindices
                sample(l)   = c%indices(i,map_n,l)
             end do
          else
             do l = 1, c%nindices
                sample(l)   = c%indices(i,2,l)
             end do
          end if

          ! Define the model to toss into the likelihood evaluation
          ! Ensure proper handling of poltypes
          if (mh_mode == 1) then
             do j = 1, nbands
                model(i,map_n,j) = c%eval_signal(j,i,map_n,sample)
             end do
          else if (mh_mode == 2) then
             do k = 2, 3
                do j = 1, nbands
                   model(i,k,j) = c%eval_signal(j,i,k,sample)
                end do
             end do
          else if (mh_mode == 3) then
             do k = 1, 3
                do j = 1, nbands
                   model(i,k,j) = c%eval_signal(j,i,k,sample)
                end do
             end do
          end if
          
          ! Evaluate the lnL (already includes the -0.5 out front)
          lnl = evaluate_lnL(data,ddata%rms_map,model,map_n,i,ddata%masks(:,1))
          
          if (trim(c%prior_type(nind)) == 'gaussian') then
             lnl_old = lnl + log(eval_normal_prior(sample(nind),c%gauss_prior(nind,1),c%gauss_prior(nind,2)))
          else if (trim(c%prior_type(nind)) == 'uniform') then
             lnl_old = lnl
          end if

          ! Now we do the real sampling
          do l = 1, nsample

             ! Update theta with the new sample
             theta(nind) = sample(nind) + rand_normal(0.d0,c%step_size(nind))

             ! If the prior is outside the allow prior range, move on
             if (theta(nind) > c%uni_prior(nind,2) .or. theta(nind) < c%uni_prior(nind,1)) cycle

             ! Evaluate model for likelihood evaluation
             ! Ensure proper handling of poltypes
             if (mh_mode == 1) then
                do j = 1, nbands
                   model(i,map_n,j) = c%eval_signal(j,i,map_n,theta)
                end do
             else if (mh_mode == 2) then
                do k = 2, 3
                   do j = 1, nbands
                      model(i,k,j) = c%eval_signal(j,i,k,theta)
                   end do
                end do
             else if (mh_mode == 3) then
                do k = 1, 3
                   do j = 1, nbands
                      model(i,k,j) = c%eval_signal(j,i,k,theta)
                   end do
                end do
             end if
             
             ! Evaluate likelihood of sample
             lnl = evaluate_lnL(data,ddata%rms_map,model,map_n,i,ddata%masks(:,1))
             
             ! With the prior of course
             if (trim(c%prior_type(nind)) == 'gaussian') then
                lnl_new = lnl + & 
                     & log(eval_normal_prior(theta(nind),c%gauss_prior(nind,1),c%gauss_prior(nind,2)))
             else if (trim(c%prior_type(nind)) == 'uniform') then
                lnl_new = lnl
             end if
             
             ! Accept/reject
             diff  = lnl_new - lnl_old
             ratio = exp(diff)
             if (trim(ml_mode) == 'optimize') then
                if (ratio > 1.d0) then
                   sample(nind) = theta(nind)
                   lnl_old      = lnl_new
                end if
             else if (trim(ml_mode) == 'sample') then
                call RANDOM_NUMBER(num)
                if (ratio > num) then
                   sample(nind) = theta(nind)
                   lnl_old      = lnl_new
                end if
             end if
          end do
          ! Ensure proper handling of poltypes
          if (mh_mode == 1) then
             c%indices(i,map_n,nind) = sample(nind)
          else if (mh_mode == 2) then
             c%indices(i,2:3,nind) = sample(nind)
          else if (mh_mode == 3) then
             c%indices(i,1:3,nind) = sample(nind)
          end if
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    end if

    deallocate(data,model)

  end subroutine sample_index_mh

  ! Note that the other evaluation routines do NOT include the -0.5, working to make them obsolete
  function evaluate_lnL(data,rms,model,map_n,pixel,mask) result(lnL)
    !==========================================================================
    ! Inputs:
    !         data:  array(real(dp)) - data with which we compare the model
    !         rms:   array(real(dp)) - noise associated with the data
    !         model: array(real(dp)) - the model to compare to the data
    !         map_n: integer         - poltype for likelihood evaluation
    !         pixel: integer         - pixel number for likelihood evaluation
    !
    ! Output:
    !         lnL: real(dp)
    !
    ! if pixel < 0,  evaluate full sky, otherwise sample that pixel
    ! if map_n > 0,  evaluate that poltype
    ! if map_n = -1, evaluate Q+U jointly
    ! if map_n = -2, evaluate T+Q+U jointly
    !==========================================================================
    implicit none

    real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: data, rms, model
    real(dp), dimension(0:npix-1),              intent(in) :: mask
    integer(i4b),               intent(in) :: map_n, pixel
    integer(i4b)                           :: i,j,k
    real(dp)                               :: lnL

    ! Initialize the result to null
    lnL = 0.d0

    ! For your simplest, pol type by pol type case
    if (map_n > 0) then
       if (pixel > -1) then
          do j = 1, nbands
             lnL = lnL - 0.5d0*((data(pixel,map_n,j)-model(pixel,map_n,j))/rms(pixel,map_n,j))**2
          end do
       else
          do i = 0, npix-1
             if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
             do j = 1, nbands
                lnL = lnL - 0.5d0*((data(i,map_n,j)-model(i,map_n,j))/rms(i,map_n,j))**2
             end do
          end do
       end if
    ! For Q+U joint sampling
    else if (map_n == -1) then
       if (pixel > -1) then
          do j = 1, nbands
             do k = 2, 3
                lnL = lnL - 0.5d0*((data(pixel,k,j)-model(pixel,k,j))/rms(pixel,k,j))**2
             end do
          end do
       else
          do i = 0, npix-1
             if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
             do j = 1, nbands
                do k = 2, 3
                   lnL = lnL - 0.5d0*((data(i,k,j)-model(i,k,j))/rms(i,k,j))**2
                end do
             end do
          end do
       end if
    ! For T+Q+U joint sampling
    else if (map_n == -2) then
       if (pixel > -1) then
          do j = 1, nbands
             do k = 1, 3
                lnL = lnL - 0.5d0*((data(pixel,k,j)-model(pixel,k,j))/rms(pixel,k,j))**2
             end do
          end do
       else
          do i = 0, npix-1
             if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
             do j = 1, nbands
                do k = 1, 3
                   lnL = lnL - 0.5d0*((data(i,k,j)-model(i,k,j))/rms(i,k,j))**2
                end do
             end do
          end do
       end if
    end if
  end function evaluate_lnL

  subroutine sample_band_calibrators(dpar,ddata)
    !==========================================================================
    ! Loop over each band, and check to see which bands need 
    ! gains and offsets fit
    !
    ! Gains and offsets will only ever be fit in total intensity!
    !==========================================================================
    implicit none
    type(dang_data),  intent(in) :: ddata
    type(dang_params)            :: dpar
    integer(i4b)                 :: j

    if (iter > 1) then
       do j = 1, nbands
          if (dpar%fit_gain(j)) then
             call sample_band_gain(ddata, dpar, j)!, 1)
          end if
          if (dpar%fit_offs(j)) then
             call sample_band_offset(ddata, dpar, j)!, 1)
          end if
       end do
    end if

  end subroutine sample_band_calibrators

  subroutine sample_band_gain(ddata, dpar, band, sample)
    !==========================================================================
    ! Solve a simple linear maximum likelihood equation to find the 
    ! gain for band `band`.
    !==========================================================================
    implicit none
    type(dang_data)                     :: ddata
    class(dang_params)                  :: dpar
    integer(i4b),            intent(in) :: band
    integer(i4b), optional,  intent(in) :: sample

    type(dang_comps),   pointer         :: c
    real(dp), allocatable, dimension(:) :: calibrator, mask
    real(dp), allocatable, dimension(:) :: calibratee
    real(dp), allocatable, dimension(:) :: N_inv

    real(dp)                            :: sum1, sum2, sample_vec, gain

    integer(i4b)                        :: i, j, k

    ! Allocate all of our necessary arrays
    allocate(calibrator(0:npix-1))
    allocate(calibratee(0:npix-1))
    allocate(N_inv(0:npix-1))
    allocate(mask(0:npix-1))

    mask = ddata%masks(:,1)

    ! Make sure our mask doesn't have any missvals, as we'll be multiplying by it
    do i = 0, ddata%npix-1
       if (mask(i) == missval) then
          mask(i) = 0.d0
       end if
    end do

    ! Set the map to be calibrated as a dummy array
    calibratee = ddata%sig_map(:,1,band)-ddata%offset(band)

    ! Define N_inv
    do i = 0, npix-1
       N_inv(i) = 1.d0/(ddata%rms_map(i,1,band)**2)
    end do

    ! Check to see which component we calibrate against
    if (trim(dpar%band_calib(band)) == 'all') then
       ! If all, calibrate against the sky model
       do i = 0, npix-1
          calibrator(i) = ddata%sky_model(i,1,band)
       end do
    else
       ! Find the correct component
       do i = 1, ncomp
          if (trim(component_list(i)%p%label) == trim(dpar%band_calib(band))) then
             c => component_list(i)%p
          end if
          if (.not. associated(c)) then
             write(*,*) "Error: cannot fit band gain to ", trim(dpar%band_calib(band))
             write(*,*) "Component label not found!"
             stop
          end if
       end do
       ! Set calibrator as the component map
       do i = 0, npix-1
          calibrator(i) = c%eval_signal(band,i,1)
       end do
    end if

    ! This is just a shorthand version for solving:
    !
    ! (calibrator^t N_{inv,band} calibrator_band)*gain = (calibrator*N_{inv,band}*data)
    sum1 = sum(mask*calibrator*N_inv*calibratee)
    sum2 = sum(mask*calibrator*N_inv*calibrator)


    ! Add this term to the RHS if we're sampling:
    ! (calibrator*N_{inv,band}^1/2)*eta
    if (present(sample)) then
       sample_vec = sum(mask*calibrator*sqrt(N_inv))*rand_normal(0.d0,1.d0)
       sum2 = sum2 + sample_vec
    end if

    ! And finally solve
    gain = sum1/sum2
    
    ! And save
    ddata%gain(band) = gain

    deallocate(mask,calibrator,calibratee,N_inv)

  end subroutine sample_band_gain

  subroutine sample_band_offset(ddata, dpar, band, sample)
    ! Use a linear regression to solve for the slope between two maps
    !
    implicit none
    type(dang_data)                     :: ddata
    class(dang_params)                  :: dpar
    integer(i4b),            intent(in) :: band
    integer(i4b), optional,  intent(in) :: sample

    type(dang_comps),   pointer         :: c
    real(dp), allocatable, dimension(:) :: calibrator, mask
    real(dp), allocatable, dimension(:) :: calibratee
    real(dp), allocatable, dimension(:) :: N_inv

    real(dp)                            :: sum1, sum2, sample_vec, offset

    integer(i4b)                        :: i, j, k

    ! Allocate all of our necessary arrays
    allocate(calibrator(0:npix-1))
    allocate(calibratee(0:npix-1))
    allocate(N_inv(0:npix-1))
    allocate(mask(0:npix-1))

    mask = ddata%masks(:,1)

    ! Find the correct component
    do i = 1, ncomp
       if (trim(component_list(i)%p%label) == trim(dpar%band_calib(band))) then
          c => component_list(i)%p
       end if
       if (.not. associated(c)) then
          write(*,*) "Error: cannot fit band gain to ", trim(dpar%band_calib(band))
          write(*,*) "Component label not found!"
          stop
       end if
    end do

    ! Set the map to be calibrated as a dummy array
    calibratee = ddata%gain(band)*ddata%sig_map(:,1,band)

    ! Define N_inv
    ! Make sure our mask doesn't have any missvals, as we'll be multiplying by it
    do i = 0, npix-1
       if (mask(i) == missval) then
          mask(i)       = 0.d0
          N_inv(i)      = 0.d0
          calibrator(i) = 0.d0
          calibratee(i) = 0.d0
       else 
          N_inv(i)      = 1.d0/(ddata%rms_map(i,1,band)**2)
          calibrator(i) = c%eval_signal(band,i,1)
       end if
    end do

    ! Weighted least squares
    offset = sqrt(sum(mask*((calibrator-calibratee)*N_inv)**2))/sum(mask)

    ddata%offset(band) = offset

  end subroutine sample_band_offset
  
  subroutine calc_hi_gain_offset(dpar, dat, comp, map_n, fg, band)
    class(dang_params)                  :: dpar
    type(dang_data)                     :: dat
    type(dang_comps)                    :: comp
    integer(i4b),            intent(in) :: map_n
    integer(i4b),            intent(in) :: band
    integer(i4b),            intent(in) :: fg
    real(dp), allocatable, dimension(:) :: map1, map2, mask
    real(dp), allocatable, dimension(:) :: gain, offset
    
    
    allocate(map1(0:dat%npix-1))
    allocate(map2(0:dat%npix-1))
    allocate(mask(0:dat%npix-1))
    
    allocate(gain(0:dpar%nsample))
    allocate(offset(0:dpar%nsample))
    
    gain(0)   = 1.d0
    offset(0) = 0.d0
    
    map1 = 0.d0
    map2 = 0.d0
    mask = 0.d0
    
    mask = dat%masks(:,1)
    
    ! map1 is the map we calibrate against here, being a component map.
    ! Must ensure the map is calculated prior to gain fitting
    ! for the HI fit, we use the foreground model to fit
    
    ! For the HI fit, we know that as the dust emission goes to 0, so does HI column density.
    ! This relationship is approximately linear below HI ~ 4e20 cm-2. In this fit we fit the gain
    ! to the full HI model presented, but the offset will be determined using the dust-HI relationship
    
    if (trim(dpar%mode) == 'hi_fit') then
       do i = 0, dat%npix-1
          if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
          map1(i) = comp%HI_amps(band)*comp%HI(i,1)*planck(bp(band),comp%T_d(i,1))
       end do
    else
       map1 = dat%fg_map(:,map_n,band,fg)
    end if
    
    map2 = dat%sig_map(:,map_n,band)
    
    do i = 1, dpar%nsample-1
       
       ! Fit gain to the SED first (HI in the HI fit case)
       
       ! if (trim(dpar%mode) == 'hi_fit') then
       !    offset(i) = sum(mask(:)*(map2(:) - gain(i-1)*comp%HI(:,1)))/sum(mask(:))
       ! else
       offset(i) = sum(mask(:)*(map2(:) - gain(i-1)*map1(:)))/sum(mask(:))
       ! end if                 
       ! Calculate the offset using the HI map from the calibrated band
       
       gain(i)   = sum(mask(:)*(map1(:)*(map2(:)-offset(i))))/sum(mask(:)*(map1(:)**2))
       
       if (i > dpar%nsample) then
          exit
       end if
       
       if (i > 10) then
          if (abs(gain(i)-gain(i-10))/abs(gain(i)) < 1d-6 .and. abs(offset(i) - offset(i-10)) < 1d-6) exit
       end if
    end do
    
    ! write(*,*) gain(i-1)
    ! write(*,*) offset(i-1)
    
    dat%gain(band)   = gain(i-1)
    dat%offset(band) = offset(i-1)
    
  end subroutine calc_hi_gain_offset

  function eval_jeffreys_prior(dpar, dat, comp, map_n, ind, val, pixel) result(prob)
    implicit none

    class(dang_params)                    :: dpar
    type(dang_comps),       intent(inout) :: comp
    type(dang_data)                       :: dat
    real(dp),               intent(in)    :: val
    integer(i4b),           intent(in)    :: map_n, ind
    integer(i4b), optional, intent(in)    :: pixel
    real(dp)                              :: prob, sum, ss
    integer(i4b)                          :: i, j, k
    
    prob = 0.d0
    sum = 0.d0

    if (trim(dpar%fg_label(ind)) == 'synch') then
       ! Is this evaluated for a single pixel?
       if (present(pixel)) then
          ! If map_n = -1, then sum over poltypes
          if (map_n == -1) then
             do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
                do j = 1, nbands
                   ss  = dat%fg_map(pixel,k,0,ind)*(dpar%band_nu(j)/dpar%fg_nu_ref(ind))**val
                   sum = sum + (((1.0/dat%rms_map(pixel,k,j))**2)*(ss/dat%fg_map(pixel,k,0,ind))*&
                        & log(dpar%band_nu(j)/dpar%fg_nu_ref(ind)))**2.0
                end do
             end do
          else
             do j = 1, nbands
                ss  = dat%fg_map(pixel,map_n,0,ind)*(dpar%band_nu(j)/dpar%fg_nu_ref(ind))**val
                sum = sum + (((1.0/dat%rms_map(pixel,map_n,j))**2)*(ss/dat%fg_map(pixel,map_n,0,ind))*&
                     & log(dpar%band_nu(j)/dpar%fg_nu_ref(ind)))**2.0
             end do
          end if
       else
          ! If map_n = -1, then sum over poltypes
          if (map_n == -1) then
             do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
                do i = 0, npix-1
                   if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) cycle
                   do j = 1, nbands
                      ss  = dat%fg_map(i,k,0,ind)*(dpar%band_nu(j)/dpar%fg_nu_ref(ind))**val
                      sum = sum + (((1.0/dat%rms_map(i,k,j))**2)*(ss/dat%fg_map(i,k,0,ind))*&
                           & log(dpar%band_nu(j)/dpar%fg_nu_ref(ind)))**2.0
                   end do
                end do
             end do
          else
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) cycle
                do j = 1, nbands
                   ss  = dat%fg_map(i,k,0,ind)*(dpar%band_nu(j)/dpar%fg_nu_ref(ind))**val
                   sum = sum + (((1.0/dat%rms_map(i,k,j))**2)*(ss/dat%fg_map(i,k,0,ind))*& 
                        & log(dpar%band_nu(j)/dpar%fg_nu_ref(ind)))**2.0
                end do
             end do
          end if
       end if
    end if
    prob = sqrt(sum)

  end function eval_jeffreys_prior
  
end module dang_sample_mod
