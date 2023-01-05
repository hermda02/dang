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

  subroutine sample_spectral_parameters(dpar,ddata)
    ! ============================================ |
    ! This function simply loops over each of the  |
    ! parameters, checking to see if their indices |
    ! ought to be sampled.                         |
    ! ============================================ |
    implicit none

    type(dang_data)              :: ddata
    type(dang_params)            :: dpar
    type(dang_comps), pointer    :: c
    integer(i4b)                 :: i, j, k

    logical(lgt)                 :: sampled

    sampled = .false.

    ! Loop over all foregrounds and sample for indices
    do i = 1, ncomp
       c => component_list(i)%p
       ! if there are no indices, don't even look
       if (c%nindices == 0) cycle
       ! if none get sampled, don't even look
       if  (any(c%sample_index)) then
          sampled = .true.
       else
          cycle
       end if
       do j = 1, c%nindices
          if (c%sample_index(j)) then
             do k = 1, c%nflag(j)
                if (iand(c%pol_flag(j,k),1) .ne. 0) then
                   write(*,*) 'Sampling spectral index ', trim(c%ind_label(j)), ' for ', trim(c%label), ', poltype = I.'
                   call sample_index_mh(ddata,c,j,1)
                else if (iand(c%pol_flag(j,k),2) .ne. 0) then
                   write(*,*) 'Sampling spectral index ', trim(c%ind_label(j)), ' for ', trim(c%label), ', poltype = Q.'
                   call sample_index_mh(ddata,c,j,2)
                else if (iand(c%pol_flag(j,k),4) .ne. 0) then
                   write(*,*) 'Sampling spectral index ', trim(c%ind_label(j)), ' for ', trim(c%label), ', poltype = U.'
                   call sample_index_mh(ddata,c,j,3)
                else if (iand(c%pol_flag(j,k),8) .ne. 0) then
                   write(*,*) 'Sampling spectral index ', trim(c%ind_label(j)), ' for ', trim(c%label), ', poltype = Q+U.'
                   call sample_index_mh(ddata,c,j,-1)
                else if (iand(c%pol_flag(j,k),0) .ne. 0) then
                   write(*,*) 'Sampling spectral index ', trim(c%ind_label(j)), ' for ', trim(c%label), ', poltype = I+Q+U.'
                   call sample_index_mh(ddata,c,j,-2)
                else
                   write(*,*) "There is something wrong with the poltype flag"
                   write(*,*) "for component ", trim(c%label)
                end if
             end do
          end if
       end do
    end do

    if (sampled) then
       call ddata%update_sky_model
       call write_stats_to_term(ddata,dpar,iter)
    end if

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
    !                                                |
    ! Modes defined as follows:                      |
    ! if map_n > 0,  evaluate that poltype           |
    ! if map_n = -1, evaluate Q+U jointly            |
    ! if map_n = -2, evaluate T+Q+U jointly          |
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
    logical(lgt)                            :: sample_it
    
    real(dp), allocatable, dimension(:)     :: sample, theta
    integer(i4b),          dimension(2)     :: map_inds

    ! Arrays for the native resolution data/rms
    real(dp), allocatable, dimension(:,:,:) :: data_raw, rms_raw
    real(dp), allocatable, dimension(:,:)   :: mask_raw
    real(dp), allocatable, dimension(:,:)   :: index_full_res

    ! And what's used in the routines:
    real(dp), allocatable, dimension(:,:,:) :: data, rms
    real(dp), allocatable, dimension(:,:,:) :: model
    real(dp), allocatable, dimension(:,:)   :: mask, index_map

    integer(i4b)                            :: i, j, k, sample_npix
    integer(i4b)                            :: l, m, n

    ! Metropolis parameters
    real(dp)                                :: lnl, lnl_old, lnl_new
    real(dp)                                :: diff, ratio, num

    ! real(dp), dimension(1000)               :: theta_grid, lnl_grid

    real(dp)                                :: t1, t2, t3, t4, t5, t6 ! Timing variables

    ! This parameter is set in case we want to sample from the prior
    ! If we do, we turn this to false so we don't get into the sampling loop
    sample_it = .true.
    
    ! Little extra section here for mode determination
    !=================================================
    ! if map_n > 0,  evaluate that poltype
    ! if map_n = -1, evaluate Q+U jointly
    ! if map_n = -2, evaluate T+Q+U jointly
    !=================================================
    if (map_n == -1) then
       map_inds(1) = 2; map_inds(2) = 3
    else if (map_n == -2) then
       map_inds(1) = 1; map_inds(2) = 3
    else
       map_inds(1) = map_n; map_inds(2) = map_n
    end if
    !=================================================

    ! Here is an object we'll call data, which we will correct to be the
    ! portion of the sky signal we wish to fit to
    allocate(data_raw(0:npix-1,nmaps,nbands))
    allocate(rms_raw(0:npix-1,nmaps,nbands))
    allocate(mask_raw(0:npix-1,nmaps))
    
    ! Initialize data and model
    do j = 1, nbands
       data_raw(0:,1,j)   = ddata%sig_map(0:,1,j)/ddata%gain(j)
       data_raw(0:,2:3,j) = ddata%sig_map(0:,2:3,j)
       rms_raw(0:,:,j)    = ddata%rms_map(0:,:,j)
    end do
    mask_raw = ddata%masks

    ! Loop over components, remove all others
    do l = 1, ncomp
       c2 => component_list(l)%p
       if (c2%label /= c%label) then
          do i = 0, npix-1
             do k = 1, nmaps
                do j = 1, nbands
                   data_raw(i,k,j) = data_raw(i,k,j) - c2%eval_signal(j,i,k)
                end do
             end do
          end do
       end if
    end do

    ! Set up data structures at the appropriate nside
    sample_npix = 12*c%sample_nside(nind)**2
    allocate(data(0:sample_npix-1,nmaps,nbands))   
    allocate(rms(0:sample_npix-1,nmaps,nbands))   
    allocate(mask(0:sample_npix-1,nmaps))   

    if (nside == c%sample_nside(nind)) then
       data = data_raw
       rms  = rms_raw 
       mask = mask_raw
    else
       ! udgrade the mask
       call udgrade_mask(mask_raw,nside,mask,c%sample_nside(nind),0.5d0)
       
       ! and data and rms maps
       do j = 1, nbands
          call udgrade_ring(data_raw(:,:,j),nside,data(:,:,j),c%sample_nside(nind))
          call udgrade_rms(rms_raw(:,:,j),nside,rms(:,:,j),c%sample_nside(nind))
       end do
    end if
    deallocate(data_raw,rms_raw,mask_raw)

    ! Allocate and initialize empty arrays
    allocate(model(0:sample_npix-1,nmaps,nbands))   
    allocate(index_map(0:sample_npix-1,nmaps))   
    allocate(index_full_res(0:npix-1,nmaps))
    model(:,:,:)        = 0.d0
    index_map(:,:)      = 0.d0
    index_full_res(:,:) = 0.d0

    allocate(sample(c%nindices),theta(c%nindices))

    ! Now time to begin sampling
    !======================================================================
    ! Index mode 1 corresponds to full sky value for the spectral parameter
    if (c%index_mode(nind) == 1) then
       write(*,*) 'Sampling fullsky'

       lnl = 0.d0

       ! Initialize index map from previous Gibbs iteration
       ! Ensure proper handling of poltypes
       do l = 1, c%nindices
          sample(l) = c%indices(0,map_inds(1),l)
       end do

       ! Init theta
       theta = sample
       
       ! Define the model to toss into the likelihood evaluation
       call update_sample_model(model,c,map_inds,sample)
       
       ! Evaluate the lnL (already includes the -0.5 out front)
       if (c%lnl_type(nind) == 'chisq') then
          lnl = evaluate_lnL(data,rms,model,map_inds,-1,mask(:,1))
       else if (c%lnl_type(nind) == 'marginal') then
          lnl = evaluate_marginal_lnL(data,rms,model,map_inds,-1,mask(:,1))
       else if (c%lnl_type(nind) == 'prior') then
          sample_it = .false.
          sample(nind) = rand_normal(c%gauss_prior(nind,1),c%gauss_prior(nind,2))
       end if
       
       if (trim(c%prior_type(nind)) == 'gaussian') then
          lnl_old = lnl + log(eval_normal_prior(sample(nind),c%gauss_prior(nind,1),c%gauss_prior(nind,2)))
       else if (trim(c%prior_type(nind)) == 'uniform') then
          lnl_old = lnl
       end if

       ! Now we do the real sampling
       if (sample_it) then
          do l = 1, nsample
             
             ! Update theta with the new sample
             ! Evaluate model for likelihood evaluation
             theta(nind) = sample(nind) + rand_normal(0.d0,c%step_size(nind))
             if (theta(nind) .lt. c%uni_prior(nind,1) .or. theta(nind) .gt. c%uni_prior(nind,2)) cycle
             
             call update_sample_model(model,c,map_inds,theta)
             
             ! Evaluate the lnL (already includes the -0.5 out front)
             if (c%lnl_type(nind) == 'chisq') then
                lnl = evaluate_lnL(data,rms,model,map_inds,-1,mask(:,1))
             else if (c%lnl_type(nind) == 'marginal') then
                lnl = evaluate_marginal_lnL(data,rms,model,map_inds,-1,mask(:,1))
             end if
             
             ! Ensure prior contributes
             if (trim(c%prior_type(nind)) == 'gaussian') then
                lnl_new = lnl + log(eval_normal_prior(theta(nind),c%gauss_prior(nind,1),c%gauss_prior(nind,2)))
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
       end if

       ! Cast the final sample back to the dummy index map
       index_map(:,map_inds(1):map_inds(2)) = sample(nind)

    ! Index mode 2 corresponds to per-pixel values for the spectral parameter
    else if (c%index_mode(nind) == 2) then
       write(*,*) 'Sampling per-pixel at nside ', c%sample_nside(nind)
       ! Pixel-by-pixel

       !!$OMP PARALLEL PRIVATE(i,j,k,l,lnl,sample,theta,lnl_old,lnl_new,diff,ratio,model,num,sample_it) SHARED(index_map)
       !!$OMP DO SCHEDULE(static)
       do i = 0, sample_npix-1
          if (mask(i,1) == 0.d0 .or. mask(i,1) == 0.d0) cycle

          ! Initialize the MH chain
          lnl      = 0.d0

          ! Ensure proper handling of poltypes
          do l = 1, c%nindices
             sample(l) = c%indices(i,map_inds(1),l)
          end do

          ! Init theta
          theta = sample

          ! Define the model to toss into the likelihood evaluation
          call update_sample_model(model,c,map_inds,sample,i)

          ! Evaluate the lnL (already includes the -0.5 out front)
          if (c%lnl_type(nind) == 'chisq') then
             lnl = evaluate_lnL(data,rms,model,map_inds,i,mask(:,1))
          else if (c%lnl_type(nind) == 'marginal') then
             lnl = evaluate_marginal_lnL(data,rms,model,map_inds,i,mask(:,1))
          else if (c%lnl_type(nind) == 'prior') then
             sample_it = .false.
             sample(nind) = rand_normal(c%gauss_prior(nind,1),c%gauss_prior(nind,2))
          end if

          if (trim(c%prior_type(nind)) == 'gaussian') then
             lnl_old = lnl + log(eval_normal_prior(sample(nind),c%gauss_prior(nind,1),c%gauss_prior(nind,2)))
          else if (trim(c%prior_type(nind)) == 'uniform') then
             lnl_old = lnl
          end if

          ! Now we do the real sampling
          if (sample_it) then
             do l = 1, nsample
                
                ! Update theta with the new sample
                ! Evaluate model for likelihood evaluation
                theta(nind) = sample(nind) + rand_normal(0.d0,c%step_size(nind))
                if (theta(nind) .lt. c%uni_prior(nind,1) .or. theta(nind) .gt. c%uni_prior(nind,2)) cycle
                call update_sample_model(model,c,map_inds,theta,i)
                
                ! Evaluate likelihood of sample
                if (c%lnl_type(nind) == 'chisq') then
                   lnl = evaluate_lnL(data,rms,model,map_inds,i,mask(:,1))
                else if (c%lnl_type(nind) == 'marginal') then
                   lnl = evaluate_marginal_lnL(data,rms,model,map_inds,i,mask(:,1))
                end if
                
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
          end if
          ! Ensure proper handling of poltypes
          ! Cast the final sample back to the component index map
          index_map(i,map_inds(1):map_inds(2)) = sample(nind)
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL
       !!$OMP BARRIER
       call udgrade_ring(index_map,c%sample_nside(nind),index_full_res,nside)
       c%indices(:,map_inds(1):map_inds(2),nind) = index_full_res(:,map_inds(1):map_inds(2))
    end if

  end subroutine sample_index_mh

  subroutine sample_calibrators(ddata)
    !=======================================================================|
    !                                                                       |  
    ! Loop through bands determining the calibration for each band.         |
    ! Determines both gain and offset relative to the full sky model.       |
    !                                                                       |
    !=======================================================================|
    implicit none
    
    type(dang_data),             intent(in) :: ddata
    integer(i4b)                            :: j

    if (any(ddata%fit_gain(:)) .or. any(ddata%fit_offset(:))) write(*,*) "Sampling band calibrators"
    do j = 1, nbands
       if (ddata%fit_gain(j)) then
          call fit_band_gain(ddata, 1, j)
       end if
       if (ddata%fit_offset(j)) then
          call fit_band_offset(ddata, 1, j)
       end if
    end do

  end subroutine sample_calibrators

  subroutine update_sample_model(model,c,map_inds,sample,pixel)
    !=======================================================================|
    !                                                                       |  
    ! This routine simply takes in an array, the component                  |
    ! pointer, and an indices array, and creates a model.                   |
    !                                                                       |
    !=======================================================================|
    ! Inputs:                                                               |
    !    model: array(real(dp)) - the model to compare to the data          |
    !    c: type(dang_comps) - pointer to component                         |
    !    map_inds: array(integer) - poltype for likelihood evaluation       |
    !    sample: array(real(dp))  - array of spectral indices for component |
    !    pixel: integer           - pixel number for likelihood evaluation  |
    !                                                                       |
    !=======================================================================|

    implicit none

    real(dp),  dimension(0:,:,:), intent(inout) :: model 
    type(dang_comps),    pointer, intent(in)    :: c
    integer(i4b),   dimension(2), intent(in)    :: map_inds
    real(dp),       dimension(:), intent(in)    :: sample
    integer(i4b),   optional,     intent(in)    :: pixel

    integer(i4b)                                :: i, j, k, sample_npix

    sample_npix = size(model,DIM=1)

    if (present(pixel)) then
       do k = map_inds(1), map_inds(2)
          do j = 1, nbands
             model(pixel,k,j) = c%eval_signal(j,pixel,k,sample)
          end do
       end do
    else
       !$OMP PARALLEL PRIVATE(i,j,k)
       !$OMP DO SCHEDULE(static)
       do i = 0, sample_npix-1
          do k = map_inds(1), map_inds(2)
             do j = 1, nbands
                model(i,k,j) = c%eval_signal(j,i,k,sample)
             end do
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    end if

  end subroutine update_sample_model

  function evaluate_marginal_lnL(data,rms,model,map_inds,pixel,mask) result(lnL)
    !==========================================================================++|
    !                                                                            |
    ! For this implementation, we want to marginalize over the foreground        |
    ! amplitude values in order to draw from a more broad distribution of        |
    ! spectral parameters.                                                       |
    !                                                                            |
    ! The equation we wish to solve for in the likelihood evaluation here is:    |
    !                                                                            |
    ! lnL(beta) = sum_nu -0.5*(T_nu^t N_nu^-1 d_nu)^t ( T_nu^t N_nu^-1 T_nu)^-1  |
    !                 (T_nu^t N_nu^-1 d_nu) + 0.5 ln |(T_nu^t N_nu^-1 T_nu)^-1|  |
    !                                                                            |
    ! This technique was used in the BeyondPlanck analysis here:                 |
    !                  https://arxiv/org/abs/2201.08188                          |
    !                                                                            |
    !============================================================================|
    ! Inputs:                                                                    |
    !         data:  array(real(dp))   - data with which we compare the model    |
    !         rms:   array(real(dp))   - noise associated with the data          |
    !         model: array(real(dp))   - the model to compare to the data        |
    !         map_inds: array(integer) - poltype for likelihood evaluation       |
    !         pixel: integer           - pixel number for likelihood evaluation  |
    !                                                                            |
    ! Output:                                                                    |
    !         lnL: real(dp)                                                      |
    !                                                                            |
    ! if pixel < 0,  evaluate full sky, otherwise sample that pixel              |
    ! if map_n > 0,  evaluate that poltype                                       |
    ! if map_n = -1, evaluate Q+U jointly                                        |
    ! if map_n = -2, evaluate T+Q+U jointly                                      |
    !============================================================================|
    implicit none

    real(dp), dimension(0:,:,:),  intent(in) :: data, rms, model
    real(dp), dimension(0:),      intent(in) :: mask
    integer(i4b),   dimension(2), intent(in) :: map_inds
    integer(i4b),                 intent(in) :: pixel
    integer(i4b)                             :: i,j,k
    real(dp)                                 :: lnL
    real(dp)                                 :: TNd, TNT, invTNT

    real(dp), allocatable, dimension(:)      :: TN

    integer(i4b), dimension(2,2)             :: inds


    ! Initialize the result to null
    lnL = 0.d0

    ! Initialize the inds array to condense the following lines:
    inds(1,:) = map_inds

    if (pixel > -1) then
       ! For a specific pixel
       inds(2,:) = pixel
    else
       ! For all pixels
       inds(2,1) = 0; inds(2,2) = size(data,DIM=1)-1
    end if

    ! So for the pixel-by-pixel case, this whole thing will be easy because all of 
    ! our crafted matrices will be fully diagonal. Looking at this on paper it 
    ! seems that unless there are cross-band template terms.
    ! TN = model*noise^{-1}
    ! TNd = TN*data

    do j = 1, nbands
       do k = inds(1,1), inds(1,2)
          TN     = model(inds(2,1):inds(2,2),k,j)/rms(inds(2,1):inds(2,2),k,j)**2
          TNd    = sum(TN*data(inds(2,1):inds(2,2),k,j))
          TNT    = sum(TN*model(inds(2,1):inds(2,2),k,j))
          invTNT = 1.d0/TNT
          ! write(*,*) j, k, invTNT, TNd
          ! write(*,*) " ", -0.5d0*TNd*invTNT*TNd
          ! write(*,*) " ", -0.5d0*sum(data(inds(2,1):inds(2,2),k,j)**2/rms(inds(2,1):inds(2,2),k,j)**2)
          ! write(*,*) " ", 0.5*log(invTNT)
          ! write(*,*) ""

          lnL    = lnL - 0.5d0*TNd*invTNT*TNd

          ! If we use the determinant
          ! lnL = lnL + 0.5*log(invTNT)
       end do
    end do
    ! stop

  end function evaluate_marginal_lnL

  function evaluate_lnL(data,rms,model,map_inds,pixel,mask) result(lnL)
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

    real(dp), dimension(0:,:,:),  intent(in) :: data, rms, model
    real(dp), dimension(0:),      intent(in) :: mask
    integer(i4b),   dimension(2), intent(in) :: map_inds
    integer(i4b),                 intent(in) :: pixel
    integer(i4b)                             :: i,j,k
    integer(i4b),   dimension(2,2)           :: inds
    real(dp)                                 :: lnL

    ! Initialize the result to null
    lnL = 0.d0

    ! Initialize the inds array to condense the following lines:
    inds(1,:) = map_inds

    if (pixel > -1) then
       ! For a specific pixel
       inds(2,:) = pixel
    else
       ! For all pixels
       inds(2,1) = lbound(data,DIM=1); inds(2,2) = ubound(data,DIM=1)
    end if

    do i = inds(2,1), inds(2,2)
       if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
       do k = inds(1,1), inds(1,2)
          do j = 1, nbands
             lnL = lnL - 0.5d0*((data(i,k,j)-model(i,k,j))/rms(i,k,j))**2
          end do
       end do
    end do

  end function evaluate_lnL
  
  subroutine fit_band_gain(ddata, map_n, band)
    
    implicit none
    type(dang_data)                     :: ddata
    integer(i4b),            intent(in) :: map_n
    integer(i4b),            intent(in) :: band
    real(dp), allocatable, dimension(:) :: map1, map2, mask, noise, N_inv
    real(dp)                            :: norm, gain
    
    allocate(map1(0:ddata%npix-1))
    allocate(map2(0:ddata%npix-1))
    allocate(mask(0:ddata%npix-1))
    allocate(noise(0:ddata%npix-1))
    allocate(N_inv(0:ddata%npix-1))
    
    map1 = 0.d0
    map2 = 0.d0
    mask = 0.d0
    
    mask = ddata%masks(:,1)

    ! Make sure our mask doesn't have any missvals, as we'll be multiplying by it
    do i = 0, ddata%npix-1
       if (mask(i) == missval) then
          mask(i) = 0.d0
       end if
    end do
    
    ! map1 is the map we calibrate against here, being the full sky model.
    map1 = ddata%sky_model(:,map_n,band)
    
    map2  = ddata%sig_map(:,map_n,band)-ddata%offset(band)
    noise = ddata%rms_map(:,map_n,band)
    N_inv = 1.d0/(noise**2)
    
    ! Super simple - find the multiplicative factor by finding the maximum likelihood
    ! solution through a linear fit to the foreground map.
    gain = sum(mask*map1*map2)/sum(mask*map1*map1)
    
    ! Save to data type variable corresponding to the band.
    ddata%gain(band) = gain
    
  end subroutine fit_band_gain
  
  subroutine fit_band_offset(ddata, map_n, band)
    type(dang_data)                     :: ddata
    integer(i4b),            intent(in) :: map_n
    integer(i4b),            intent(in) :: band
    real(dp), allocatable, dimension(:) :: map1, map2, mask, noise, N_inv
    real(dp)                            :: norm, offset
    
    allocate(map1(0:ddata%npix-1))
    allocate(map2(0:ddata%npix-1))
    allocate(mask(0:ddata%npix-1))
    allocate(noise(0:ddata%npix-1))
    allocate(N_inv(0:ddata%npix-1))
    
    mask = ddata%masks(:,1)

    ! Make sure our mask doesn't have any missvals, as we'll be multiplying by it
    do i = 0, ddata%npix-1
       if (mask(i) == missval) then
          mask(i) = 0.d0
       end if
    end do
    
    ! map1 is the map we calibrate against here, being the full sky model.
    map1 = ddata%sky_model(:,map_n,band)
    map2 = ddata%sig_map(:,map_n,band)
    
    noise = ddata%rms_map(:,map_n,band)
    N_inv = 1.d0/(noise**2)
    
    offset = sum(mask(:)*(map2(:) - ddata%gain(band)*map1(:)))/sum(mask(:))
        
    ddata%offset(band) = offset
    
  end subroutine fit_band_offset

  ! NEEDS A FULL REWRITE FOR HOW COMPONENTS WORK NOW
  ! function eval_jeffreys_prior(dpar, dat, comp, map_n, ind, val, pixel) result(prob)
  !   implicit none

  !   class(dang_params)                    :: dpar
  !   type(dang_comps),       intent(inout) :: comp
  !   type(dang_data)                       :: dat
  !   real(dp),               intent(in)    :: val
  !   integer(i4b),           intent(in)    :: map_n, ind
  !   integer(i4b), optional, intent(in)    :: pixel
  !   real(dp)                              :: prob, sum, ss
  !   integer(i4b)                          :: i, j, k
    
  !   prob = 0.d0
  !   sum = 0.d0

  !   if (trim(dpar%fg_label(ind)) == 'synch') then
  !      ! Is this evaluated for a single pixel?
  !      if (present(pixel)) then
  !         ! If map_n = -1, then sum over poltypes
  !         if (map_n == -1) then
  !            do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
  !               do j = 1, nbands
  !                  ss  = dat%fg_map(pixel,k,0,ind)*(dpar%band_nu(j)/dpar%fg_nu_ref(ind))**val
  !                  sum = sum + (((1.0/dat%rms_map(pixel,k,j))**2)*(ss/dat%fg_map(pixel,k,0,ind))*&
  !                       & log(dpar%band_nu(j)/dpar%fg_nu_ref(ind)))**2.0
  !               end do
  !            end do
  !         else
  !            do j = 1, nbands
  !               ss  = dat%fg_map(pixel,map_n,0,ind)*(dpar%band_nu(j)/dpar%fg_nu_ref(ind))**val
  !               sum = sum + (((1.0/dat%rms_map(pixel,map_n,j))**2)*(ss/dat%fg_map(pixel,map_n,0,ind))*&
  !                    & log(dpar%band_nu(j)/dpar%fg_nu_ref(ind)))**2.0
  !            end do
  !         end if
  !      else
  !         ! If map_n = -1, then sum over poltypes
  !         if (map_n == -1) then
  !            do k = dpar%pol_type(1), dpar%pol_type(size(dpar%pol_type))
  !               do i = 0, npix-1
  !                  if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) cycle
  !                  do j = 1, nbands
  !                     ss  = dat%fg_map(i,k,0,ind)*(dpar%band_nu(j)/dpar%fg_nu_ref(ind))**val
  !                     sum = sum + (((1.0/dat%rms_map(i,k,j))**2)*(ss/dat%fg_map(i,k,0,ind))*&
  !                          & log(dpar%band_nu(j)/dpar%fg_nu_ref(ind)))**2.0
  !                  end do
  !               end do
  !            end do
  !         else
  !            do i = 0, npix-1
  !               if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) cycle
  !               do j = 1, nbands
  !                  ss  = dat%fg_map(i,k,0,ind)*(dpar%band_nu(j)/dpar%fg_nu_ref(ind))**val
  !                  sum = sum + (((1.0/dat%rms_map(i,k,j))**2)*(ss/dat%fg_map(i,k,0,ind))*& 
  !                       & log(dpar%band_nu(j)/dpar%fg_nu_ref(ind)))**2.0
  !               end do
  !            end do
  !         end if
  !      end if
  !   end if
  !   prob = sqrt(sum)

  ! end function eval_jeffreys_prior
  
end module dang_sample_mod
