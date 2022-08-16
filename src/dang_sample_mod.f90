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
             ! Should probably include some modes for varying
             ! poltype combinations
             ! if (c%joint_index(l)) then
             ! 
             write(*,*) 'Sampling spectral index for ', trim(c%label)
             call sample_index_mh(ddata,c,j,-1)
          end if
       end do
    end do

  end subroutine sample_spectral_parameters
  
  subroutine sample_index_mh(ddata,c,nind,map_n)
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
          else
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
    !=================================================
    if (map_n == -1) then
       mh_mode = 2
    else
       mh_mode = 1
    end if
    
    ! Now time to begin sampling
    if (c%index_mode(nind) == 1) then
       lnl = 0.d0
       ! Ensure proper handling of poltypes
       if (mh_mode == 1) then
          do l = 1, c%nindices
             sample(l)   = c%indices(0,map_n,l)
          end do
       else
          do l = 1, c%nindices
             sample(l)   = c%indices(0,2,l)
          end do
       end if
       
       ! Define the model to toss into the likelihood evaluation
       ! Ensure proper handling of poltypes
       if (mh_mode == 1) then
          do i = 0, npix-1
             do j = 1, nbands
                model(i,map_n,j) = c%amplitude(i,map_n)*c%eval_sed(j,i,map_n,sample)
             end do
          end do
       else
          do i = 0, npix-1
             do k = 1, nmaps
                do j = 1, nbands
                   model(i,k,j) = c%amplitude(i,k)*c%eval_sed(j,i,k,sample)
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
          
          ! Evaluate model for likelihood evaluation
          ! Ensure proper handling of poltypes
          if (mh_mode == 1) then
             do i = 0, npix-1
                do j = 1, nbands
                   model(i,map_n,j) = c%amplitude(i,map_n)*c%eval_sed(j,i,map_n,theta)
                end do
             end do
          else
             do i = 0, npix-1
                do k = 1, nmaps
                   do j = 1, nbands
                      model(i,k,j) = c%amplitude(i,k)*c%eval_sed(j,i,k,theta)
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
       stop
       ! Ensure proper handling of poltypes
       if (mh_mode == 1) then
          c%indices(:,map_n,nind) = sample(nind)
       else
          c%indices(:,2:3,nind)   = sample(nind)
       end if

    else if (c%index_mode(nind) == 2) then
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
                model(i,map_n,j) = c%amplitude(i,map_n)*c%eval_sed(j,i,map_n,sample)
             end do
          else
             do k = 1, nmaps
                do j = 1, nbands
                   model(i,k,j) = c%amplitude(i,k)*c%eval_sed(j,i,k,sample)
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
             theta(nind) = sample(nind) + rand_normal(0.d0,c%gauss_prior(nind,2))

             ! Evaluate model for likelihood evaluation
             ! Ensure proper handling of poltypes
             if (mh_mode == 1) then
                do j = 1, nbands
                   model(i,map_n,j) = c%amplitude(i,map_n)*c%eval_sed(j,i,map_n,theta)
                end do
             else
                do k = 1, nmaps
                   do j = 1, nbands
                      model(i,k,j) = c%amplitude(i,k)*c%eval_sed(j,i,k,theta)
                   end do
                end do
             end if
             
             ! Evaluate
             lnl = evaluate_lnL(data,ddata%rms_map,model,map_n,i,ddata%masks(:,1))
             
             ! Accept/reject
             if (trim(c%prior_type(nind)) == 'gaussian') then
                lnl_new = lnl + & 
                     & log(eval_normal_prior(theta(nind),c%gauss_prior(nind,1),c%gauss_prior(nind,2)))
             else if (trim(c%prior_type(nind)) == 'uniform') then
                lnl_new = lnl
             end if
             
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
          else
             c%indices(i,2:3,nind) = sample(nind)
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

  ! function sample_fg_amp(dpar, dat, comp, ind, map_n)
  !   !------------------------------------------------------------------------
  !   ! Samples spectral amplitude (per pixel), following the spectrum of foreground 'type'. Returns a full map of amplitudes.
  !   !------------------------------------------------------------------------
  !   implicit none
    
  !   class(dang_params)                   :: dpar
  !   type(dang_comps)                     :: comp
  !   type(dang_data)                      :: dat
  !   integer(i4b),             intent(in) :: ind
  !   integer(i4b),             intent(in) :: map_n
  !   integer(i4b)                         :: f
  !   real(dp)                             :: sum1, sum2, spec
  !   real(dp)                             :: amp, num, t, sam
  !   real(dp), dimension(0:npix-1,nbands) :: map2fit
  !   real(dp), dimension(0:npix-1)        :: sample_fg_amp, norm
    
  !   map2fit = dat%sig_map(:,map_n,:)

  !   norm = 0.d0

  !   ! remove all other fg signals
  !   do f = 1, nfgs
  !      if (f /= ind) then
  !         map2fit(:,:) = map2fit(:,:) - dat%fg_map(:,map_n,1:,f)
  !      end if
  !   end do
    
  !   ! sum_nu ((T_nu)^T N_nu^-1 T_nu)amp = sum_nu ((T_nu)^T N_nu^-1 d_nu)  |
  !   ! sum_nu ((T_nu)^T N_nu^-1 T_nu)amp = sum_nu ((T_nu)^T N_nu^-1 d_nu)  + (T_nu)^T N_nu^{-1/2} eta|
    
  !   do i = 0, npix-1
  !      sum1    = 0.0d0
  !      sum2    = 0.0d0
  !      do j = 1, nbands
  !         spec    = compute_spectrum(dpar,comp,bp(j),ind,i,map_n)
  !         sum1    = sum1 + (map2fit(i,j)*spec)/dat%rms_map(i,map_n,j)**2.d0
  !         sum2    = sum2 + (spec)**2.d0/dat%rms_map(i,map_n,j)**2.d0
  !         norm(i) = norm(i) + spec/dat%rms_map(i,map_n,j)
  !      end do
  !      if (trim(dpar%ml_mode) == 'sample') then
  !         amp        = sum1/sum2 + rand_normal(0.d0,1.d0)*norm(i)/sum2
  !      else if (trim(dpar%ml_mode) == 'optimize') then
  !         amp        = sum1/sum2
  !      end if
  !      sample_fg_amp(i) = amp
  !   end do
  ! end function sample_fg_amp

  subroutine template_fit(dpar, dat, comp, map_n, temp_num)
    !------------------------------------------------------------------------
    ! Simple linear fit of a template to data map with a sampling term
    !------------------------------------------------------------------------
    implicit none
    
    type(dang_params)                       :: dpar
    type(dang_comps)                        :: comp
    type(dang_data)                         :: dat
    real(dp), dimension(0:npix-1,nbands)    :: cov, nos, map
    real(dp), allocatable, dimension(:,:,:) :: map2fit
    integer(i4b),                intent(in) :: map_n
    integer(i4b), optional,      intent(in) :: temp_num
    real(dp)                                :: temp, sum1, sum2, norm
    integer(i4b)                            :: i, j, k, n

    nos = dat%rms_map(:,map_n,:)
    cov = nos**2.d0
    allocate(map2fit(0:npix-1,nmaps,nbands))
    
    map2fit = dat%sig_map
    
    if (trim(dpar%mode) == 'comp_sep') then
       write(*,*) "Sampling for template "//trim(dpar%temp_label(temp_num))//", pol = "//trim(tqu(map_n))//"."
       do j = 1, dpar%numinc
          sum1 = 0.d0
          sum2 = 0.d0
          norm = 0.d0
          temp = 0.d0
          if (dpar%temp_corr(temp_num,j)) then
             ! Remove the other foregrounds from the input data before fitting
             do n = 1, dpar%ncomp
                map2fit(:,:,j) = map2fit(:,:,j) - dat%fg_map(:,:,j,n)
             end do
             do n = 1, dpar%ntemp
                if (n == temp_num) cycle
                map2fit(:,:,j) = map2fit(:,:,j) - dat%fg_map(:,:,j,dpar%ncomp+n)
             end do
             ! Calculate template amplitude
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) cycle
                sum1 = sum1 + ((map2fit(i,map_n,j)*dat%temps(i,map_n,temp_num))/cov(i,j))
                sum2 = sum2 + (dat%temps(i,map_n,temp_num)**2.d0)/cov(i,j)
                norm = norm + dat%temps(i,map_n,temp_num)/dat%rms_map(i,map_n,j)
             end do
             if (trim(dpar%ml_mode) == 'sample') then
                dat%temp_amps(j,map_n,temp_num) = sum1/sum2 + rand_normal(0.d0,1.d0)/sqrt(norm)
             else if (trim(dpar%ml_mode) == 'optimize') then
                dat%temp_amps(j,map_n,temp_num) = sum1/sum2
             end if
          end if
       end do
    else if (trim(dpar%mode) == 'hi_fit') then
       do j = 1, dpar%numinc
          sum1 = 0.d0
          sum2 = 0.d0
          norm = 0.d0
          temp = 0.d0
          if (dpar%temp_corr(1,j)) then
             do i = 0, npix-1
                if (comp%HI(i,1) > dpar%thresh) cycle
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) cycle
                temp = comp%HI(i,1)*planck(bp(j),comp%T_d(i,1))
                sum1 = sum1 + (((dat%sig_map(i,map_n,j)-dat%offset(j))/dat%gain(j))*temp)/dat%rms_map(i,map_n,j)**2.d0
                sum2 = sum2 + (temp)**2.d0/dat%rms_map(i,map_n,j)**2.d0
                norm = norm + temp/dat%rms_map(i,map_n,j)
             end do
          end if

          if (trim(dpar%ml_mode) == 'sample') then
             sum1 = sum1 + norm*rand_normal(0.d0,1.d0)
          end if
          comp%HI_amps(j) = sum1/sum2
      end do
    end if
    
  end subroutine template_fit

  subroutine sample_index_metropolis(data,rms,model,mask,theta_0,p_theta,step,map_n,pixel,lnL,sample)
    implicit none

    real(dp), dimension(0:npix-1,nmaps,nbands), intent(in)    :: data, rms
    real(dp), dimension(2),                     intent(in)    :: p_theta
    real(dp), dimension(0:npix-1),              intent(in)    :: mask
    real(dp),                                   intent(inout) :: theta_0
    real(dp),                                   intent(inout) :: lnL
    real(dp),                                   intent(in)    :: step
    integer(i4b),                               intent(in)    :: pixel, map_n
    logical(lgt),                               intent(in)    :: sample

    real(dp), dimension(0:npix-1,nmaps,nbands)                :: model_map
    real(dp)                                                  :: t, a, b, lnL_sample
    real(dp)                                                  :: diff, ratio, num, lnl_new
    integer(i4b)                                              :: j, k

    interface 
       function model(pixel,map_n,band,theta) 
         use healpix_types
         implicit none
         real(dp),     intent(in) :: theta
         integer(i4b), intent(in) :: pixel, map_n, band
         real(dp)                 :: model
       end function model
    end interface

    lnl_new = 0.d0
    ! if map_n == -1, sample Q and U jointly
    if (map_n == -1) then
       do k = 2, 3
          ! Draw a sample
          t = theta_0 + rand_normal(0.d0,step)
          ! Evaluate model
          do j = 1, nbands
             model_map(pixel,k,j) = model(pixel,k,j,t)
          end do
          ! Assess fit and consider the prior
          lnl_new = lnl_new + evaluate_lnL(data,rms,model_map,k,pixel,mask) + &
               & log(eval_normal_prior(t,p_theta(1), p_theta(2)))
       end do
    else
       ! Draw a sample
       t = theta_0 + rand_normal(0.d0,step)
       ! Evaluate model
       do j = 1, nbands
          model_map(pixel,map_n,j) = model(pixel,map_n,j,t)
       end do
       ! Assess fit and consider the prior
       lnl_new = evaluate_lnL(data,rms,model_map,map_n,pixel,mask) + &
            & log(eval_normal_prior(t,p_theta(1), p_theta(2)))
    end if
    ! Compare to previous fit
    diff = lnl_new - lnL
    ratio = exp(diff)
    
    if (sample) then
       call RANDOM_NUMBER(num)
       if (ratio > num) then
          theta_0  = t
          lnL      = lnl_new
       end if
    else 
       if (ratio > 1.d0) then
          theta_0  = t
          lnl      = lnl_new
       end if
    end if
    
  end subroutine sample_index_metropolis
    
  ! This architecture of this function has not been verified yet
  subroutine sample_HI_T(dpar, ddata, dcomps, map_n, output)
    implicit none
    
    class(dang_params)                         :: dpar
    type(dang_data)                            :: ddata
    type(dang_comps)                           :: dcomps
    integer(i4b), intent(in)                   :: map_n
    integer(i4b)                               :: nside1, npix2, nside2
    real(dp), dimension(0:npix-1,nmaps,nbands) :: cov, model
    real(dp), dimension(0:npix-1,nmaps)        :: te
    real(dp), dimension(0:npix-1)              :: te_sample
    real(dp), allocatable, dimension(:,:,:)    :: maps_low, cov_low
    real(dp), allocatable, dimension(:,:)      :: T_low
    real(dp), allocatable, dimension(:)        :: sample_T_low
    real(dp), dimension(2)                     :: p_Td
    real(dp)                                   :: Td, lnl, s

    integer(i4b), optional, intent(in)         :: output

    logical(lgt)                               :: test

    test = .false.
    if (present(output)) then
       test = .true.
    end if
    
    te      = dcomps%T_d
    cov     = ddata%rms_map*ddata%rms_map

    p_Td(1) = dpar%HI_Td_mean
    p_Td(2) = dpar%HI_Td_std
    
    nside2 = 16
    
    nside1 = npix2nside(npix)
    if (nside1 == nside2) then
       npix2 = npix
    else
       npix2 = nside2npix(nside2)
    end if
    allocate(maps_low(0:npix2-1,nmaps,nbands))
    allocate(T_low(0:npix2-1,nmaps),cov_low(0:npix2-1,nmaps,nbands))
    allocate(sample_T_low(0:npix2-1))
    
    do j = 1, nbands
       maps_low(:,:,j)   = ddata%sig_map(:,:,j)
       cov_low(:,:,j)    = cov(:,:,j)
    end do
    T_low = te
    
    ! Metropolis algorithm
    
    !!$OMP PARALLEL PRIVATE(i,j,k,l,Td,lnl,t)
    !!$OMP DO SCHEDULE(STATIC)
    do i = 0, npix2-1
       if (ddata%masks(i,1) == 0.d0  .or. ddata%masks(i,1) == missval) then
          sample_T_low(i) = missval
          cycle
       else
          Td = T_low(i,map_n)
          ! Chi-square from the most recent Gibbs chain update
          do j = 1, nbands
             model(i,map_n,j) = model_HI_Td(i,map_n,j,Td)
          end do
          
          lnl = evaluate_lnL(maps_low,ddata%rms_map,model,1,i,ddata%masks(:,1)) + &
               & log(eval_normal_prior(Td,dpar%HI_Td_mean, dpar%HI_Td_std))

          s = dpar%HI_Td_step
          do l = 1, dpar%nsample
             call sample_index_metropolis(maps_low,ddata%rms_map,model_HI_Td,ddata%masks(:,1),Td,p_Td,s,map_n,i,lnl,.true.)
          end do
          sample_T_low(i) = Td
       end if
       ! stop
    end do
    !!$OMP END DO
    !!$OMP END PARALLEL

    ! if (nside1 /= nside2) then
    !    if (ordering == 1) then
    !       call udgrade_ring(sample_T_low, nside2, te_sample, nside1)
    !       !call convert_nest2ring(nside2, sample_T_low)
    !    else
    !       call udgrade_nest(sample_T_low, nside2, te_sample, nside1)
    !    end if
    ! else
    te_sample =  sample_T_low
    ! end if
    dcomps%T_d(:,1) = te_sample
    
    deallocate(maps_low)
    deallocate(T_low)
    deallocate(cov_low)
    deallocate(sample_T_low)

  contains

    function model_HI_Td(pixel,map_n,band,Td)
      use healpix_types
      implicit none
      integer(i4b), intent(in) :: pixel, map_n, band
      real(dp),     intent(in) :: Td

      real(dp) :: model_HI_Td
      
      model_HI_Td = ddata%gain(band)*dcomps%HI_amps(band)*dcomps%HI(pixel,1)*&
           & planck(bp(band),Td)+ddata%offset(band)

    end function model_HI_Td

  end subroutine sample_HI_T
  ! ------------------------------------------------------------
  
  subroutine sample_band_gain(dpar, dat, comp, map_n, band, fg, sample)
    class(dang_params)                  :: dpar
    type(dang_data)                     :: dat
    type(dang_comps)                    :: comp
    integer(i4b),            intent(in) :: map_n
    integer(i4b),            intent(in) :: band
    integer(i4b),            intent(in) :: fg
    integer(i4b), optional,  intent(in) :: sample
    real(dp), allocatable, dimension(:) :: map1, map2, mask, noise, N_inv
    real(dp)                            :: norm, gain

    real(dp), allocatable, dimension(:,:) :: map3
    
    allocate(map1(0:dat%npix-1))
    allocate(map2(0:dat%npix-1))
    allocate(map3(0:dat%npix-1,1))
    allocate(mask(0:dat%npix-1))
    allocate(noise(0:dat%npix-1))
    allocate(N_inv(0:dat%npix-1))
    
    map1 = 0.d0
    map2 = 0.d0
    mask = 0.d0
    
    mask = dat%masks(:,1)

    ! Make sure our mask doesn't have any missvals, as we'll be multiplying by it
    do i = 0, dat%npix-1
       if (mask(i) == missval) then
          mask(i) = 0.d0
       end if
    end do
    
    ! map1 is the map we calibrate against here, being a component map.
    ! Must ensure the map is calculated prior to gain fitting
    ! for the HI fit, we use the foreground model to fit
    
    if (trim(dpar%mode) == 'hi_fit') then
       do i = 0, dat%npix-1
          if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
          map1(i) = comp%HI_amps(band)*comp%HI(i,1)*planck(bp(band),comp%T_d(i,1))
       end do
    else
       map1 = dat%fg_map(:,map_n,band,fg)
    end if
    
    map2  = dat%sig_map(:,map_n,band)-dat%offset(band)
    noise = dat%rms_map(:,map_n,band)
    N_inv = 1.d0/(noise**2)
    
    ! Super simple - find the multiplicative factor by finding the maximum likelihood
    ! solution through a linear fit to the foreground map.
    
    gain = sum(mask*map1*N_inv*map2)/sum(mask*map1*N_inv*map1)
    norm = sqrt(sum(mask*map1*map1*N_inv))
    
    ! Sample variable can be any number. If it's present, we sample!
    if (present(sample)) then
       gain = gain + rand_normal(0.d0,1.d0)/norm
    end if

    ! Save to data type variable corresponding to the band.
    dat%gain(band) = gain
    
  end subroutine sample_band_gain
  
  subroutine sample_band_offset(dpar, dat, comp, map_n, band, fg, sample)
    class(dang_params)                  :: dpar
    type(dang_data)                     :: dat
    type(dang_comps)                    :: comp
    integer(i4b),            intent(in) :: map_n
    integer(i4b),            intent(in) :: band
    integer(i4b),            intent(in) :: fg
    integer(i4b), optional,  intent(in) :: sample
    real(dp), allocatable, dimension(:) :: map1, map2, mask
    real(dp)                            :: norm, offset
    real(dp)                            :: n, x, x2, y, y2, xy
    
    n  = 0.d0
    x  = 0.d0
    x2 = 0.d0
    y  = 0.d0
    y2 = 0.d0
    xy = 0.d0
    
    
    allocate(map1(0:dat%npix-1))
    allocate(map2(0:dat%npix-1))
    allocate(mask(0:dat%npix-1))
    
    mask = dat%masks(:,1)

    do i = 0, dat%npix-1
       if (mask(i) == missval) then
          mask(i) = 0.d0
       end if
    end do
    
    ! map1 is the map we calibrate against here, being a component map.
    ! Must ensure the map is calculated prior to gain fitting
    ! for the HI fit, we use the foreground model to fit (amplitude*HI*B_nu(T))
    
    if (trim(dpar%mode) == 'hi_fit') then
       do i = 0, dat%npix-1
          if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
          map1(i) = comp%HI(i,1)
       end do
    else
       map1 = dat%fg_map(:,map_n,band,fg)
    end if
    
    map2 = dat%sig_map(:,map_n,band)/dat%gain(band)
    
    ! offset = sum(mask(:)*(map2(:) - dat%gain(band)*map1(:)))/sum(mask(:))
    
    do i = 0, dat%npix-1
       if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
       n  = n + 1.d0
       x  = x + map1(i)
       x2 = x2 + map1(i)*map1(i)
       y  = y + map2(i)
       y2 = y2 + map2(i)*map2(i)
       xy = xy + map1(i)*map2(i)
    end do
    
    offset = (y*x2 - x*xy)/(n*x2 - x**2)
    
    
    norm   = sum(mask(:))*sum(mask(:)*(map1(:)**2))-(sum(mask(:)*map1(:))**2)
    
    if (present(sample)) then
       offset = offset + rand_normal(0.d0,1.d0)/norm
    end if
    
    dat%offset(band) = offset
    
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
                ! write(*,*) 'k = ', k
                do j = 1, nbands
                   ! write(*,*) 'j = ', j
                   ss  = dat%fg_map(pixel,k,0,ind)*(dpar%band_nu(j)/dpar%fg_nu_ref(ind))**val
                   sum = sum + (((1.0/dat%rms_map(pixel,k,j))**2)*(ss/dat%fg_map(pixel,k,0,ind))*&
                        & log(dpar%band_nu(j)/dpar%fg_nu_ref(ind)))**2.0
                end do
                ! write(*,*) ''
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
