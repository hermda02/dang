module dang_lnl_mod
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

  function evaluate_model_lnL(c,data,rms,model,map_inds,pixel,mask,sample,nind) result(lnL_out)
    implicit none
    type(dang_comps),    pointer, intent(in) :: c
    real(dp), dimension(0:,:,:),  intent(in) :: data, rms, model
    integer(i4b),   dimension(2), intent(in) :: map_inds
    integer(i4b),   optional,     intent(in) :: pixel
    real(dp), dimension(0:),      intent(in) :: mask
    integer(i4b),                 intent(in) :: nind
    real(dp),       dimension(:), intent(in) :: sample

    real(dp)                                 :: lnl_temp, lnl_out

    ! Evaluate the lnL (already includes the -0.5 out front)
    if (c%lnl_type(nind) == 'chisq') then
       lnl_temp = evaluate_lnL(data,rms,model,map_inds,pixel,mask)
    else if (c%lnl_type(nind) == 'marginal') then
       lnl_temp = evaluate_marginal_lnL(data,rms,model,map_inds,pixel,mask)
    end if
    
    if (trim(c%prior_type(nind)) == 'gaussian') then
       lnl_out = lnl_temp + log(eval_normal_prior(sample(nind),c%gauss_prior(nind,1),c%gauss_prior(nind,2)))
    else if (trim(c%prior_type(nind)) == 'uniform') then
       lnl_out = lnl_temp
    end if

  end function evaluate_model_lnL

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

          lnL    = lnL - 0.5d0*TNd*invTNT*TNd
       end do
    end do

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
    real(dp)                                 :: lnL, lnL_local

    ! Initialize the result to null
    lnL = 0.d0
    lnl_local = 0.d0

    ! Initialize the inds array to condense the following lines:
    inds(1,:) = map_inds

    if (pixel > -1) then
       ! For a specific pixel
       inds(2,:) = pixel
    else
       ! For all pixels
       inds(2,1) = lbound(data,DIM=1); inds(2,2) = ubound(data,DIM=1)
    end if

    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(static)
    do i = inds(2,1), inds(2,2)
       if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
       do k = inds(1,1), inds(1,2)
          do j = 1, nbands
             lnL_local = lnL_local - 0.5d0*((data(i,k,j)-model(i,k,j))/rms(i,k,j))**2
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    lnL = lnL + lnL_local

  end function evaluate_lnL

  function evaluate_lnL_fullsky(data,rms,model,map_inds,pixel,mask) result(lnL)
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
    real(dp)                                 :: lnL, lnL_local

    ! Initialize the result to null
    lnL = 0.d0
    lnl_local = 0.d0

    ! Initialize the inds array to condense the following lines:
    inds(1,:) = map_inds

    if (pixel > -1) then
       ! For a specific pixel
       inds(2,:) = pixel
    else
       ! For all pixels
       inds(2,1) = lbound(data,DIM=1); inds(2,2) = ubound(data,DIM=1)
    end if

    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(static)
    do i = inds(2,1), inds(2,2)
       if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
       do k = inds(1,1), inds(1,2)
          do j = 1, nbands
             lnL_local = lnL_local - 0.5d0*((data(i,k,j)-model(i,k,j))/rms(i,k,j))**2
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    lnL = lnL + lnL_local

  end function evaluate_lnL_fullsky

  function eval_jeffreys_prior(c,data,rms,model,map_inds,pixel,mask,val) result(prob)
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

    type(dang_comps),    pointer, intent(in) :: c
    real(dp), dimension(0:,:,:),  intent(in) :: data, rms, model
    real(dp), dimension(0:),      intent(in) :: mask
    integer(i4b),   dimension(2), intent(in) :: map_inds
    integer(i4b),                 intent(in) :: pixel
    real(dp),                     intent(in) :: val

    integer(i4b),   dimension(2,2)           :: inds
    real(dp), dimension(2)                   :: theta
    real(dp)                                 :: prob, sum, ss
    integer(i4b)                             :: i, j, k, l

    prob = 0.d0
    sum = 0.d0

    theta(1) = val

    ! Initialize the inds array to condense the following lines:
    inds(1,:) = map_inds

    if (pixel > -1) then
       ! For a specific pixel
       inds(2,:) = pixel
    else
       ! For all pixels
       inds(2,1) = lbound(data,DIM=1); inds(2,2) = ubound(data,DIM=1)
    end if

    if (trim(c%label) == 'synch') then
       ! Is this evaluated for a single pixel?
       do i = inds(2,1), inds(2,2)
          if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
          do k = inds(1,1), inds(1,2)
             do j = 1, nbands
                ss  = c%eval_signal(j,i,k,theta)
                sum = sum + (((1.0/rms(i,k,j))**2)*(ss/c%amplitude(i,k))*&
                     & log(bp(j)%nu_c/c%nu_ref))**2.0
             end do
          end do
       end do
    end if
    prob = sqrt(sum)

  end function eval_jeffreys_prior



end module dang_lnl_mod