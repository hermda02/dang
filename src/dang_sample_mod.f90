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
  implicit none
  
  private :: i, j, k, l
  integer(i4b) :: i, j, k, l
  
contains
  
  subroutine sample_joint_amp(self, dat, compo, map_n, method)
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
    type(params)                               :: self
    type(data),                intent(inout)   :: dat
    type(component)                            :: compo
    integer(i4b),              intent(in)      :: map_n
    character(len=*),          intent(in)      :: method
    real(dp), allocatable, dimension(:,:,:)    :: map2fit
    real(dp), allocatable, dimension(:)        :: b, c, d, q
    character(len=256)                         :: title
    integer(i4b)                               :: x, y, z, w, l, m, n
    integer(i4b)                               :: nfit1, nfit2, nfit3, nfit4
    integer(i4b)                               :: info
    real(dp)                                   :: t6, t7
    
    if (rank == master) then
       !write(*,fmt='(a)') 'Starting joint sampling for synch and dust_template.'
       !write(*,fmt='(a)') 'Pol_type = ', trim(poltype(:))
       t1 = mpi_wtime()
    end if
    
    ! Load which components to jointly fit for arary allocation
    ! vv These will not change based off of components
    x = dat%npix
    y = 0
    z = self%numband
    
    allocate(map2fit(0:npix-1,nmaps,nbands))
    
    map2fit = dat%sig_map
    
    do n = 1, self%ncomp
       if (ANY(self%joint_comp == trim(self%fg_label(n)))) then
          if (self%joint_pol) then
             y = y + 2*x
          else
             y = y + x
          end if
       else
          write(*,*) "remove for foreground ", trim(self%fg_label(n))
          map2fit(:,:,:) = map2fit(:,:,:) - dat%fg_map(:,:,:,n)
       end if
    end do
    do n = 1, self%ntemp
       if (ANY(self%joint_comp == trim(self%temp_label(n)))) then
          y = y + self%temp_nfit(n)
       else
          write(*,*) "remove for template ", trim(self%temp_label(n))
          map2fit(:,:,:) = map2fit(:,:,:) - dat%fg_map(:,:,:,self%ncomp+n)
       end if
    end do
    allocate(b(y),c(y))
    
    ! Initialize arrays
    b(:)              = 0.d0
    c(:)              = 0.d0
    
    write(*,*) 'Compute RHS of matrix eqn.'
    ! Computing the LHS and RHS of the linear equation
    ! RHS
    w = 0 
    do m = 1, size(self%joint_comp)
       do n = 1, self%ncomp
          if (trim(self%joint_comp(m)) == trim(self%fg_label(n))) then
             if (.not. self%joint_pol) then
                do j = 1, z
                   do i = 1, x
                      if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) then
                         c(i) = 0.d0
                         cycle
                      else
                         c(i) = c(i) +  (map2fit(i-1,map_n,j)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n))/(dat%rms_map(i-1,map_n,j)**2.d0)
                      end if
                   end do
                end do
                w = w + x
             else if (self%joint_pol) then
                do j = 1, z
                   do i = 1, x
                      if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) then
                         c(i)   = 0.d0
                         c(x+i) = 0.d0
                         cycle
                      else                           
                         c(i)   = c(i)   + (map2fit(i-1,map_n,j)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n))/(dat%rms_map(i-1,map_n,j)**2.d0)
                         c(x+i) = c(x+i) + (map2fit(i-1,map_n+1,j)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n+1))/(dat%rms_map(i-1,map_n+1,j)**2.d0)
                      end if
                   end do
                end do
                w = w + 2*x
             end if
          end if
       end do
       do n = 1, self%ntemp
          if (trim(self%joint_comp(m)) == trim(self%temp_label(n))) then
             if (.not. self%joint_pol) then
                l = 1
                do j = 1, z
                   if (self%temp_corr(n,j)) then
                      do i = 1, x
                         if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                         c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*map2fit(i-1,map_n,j)*&
                              dat%temps(i-1,map_n,n)
                      end do
                      l = l + 1
                   end if
                end do
                w = w + self%temp_nfit(n)
             else if (self%joint_pol) then
                ! If sampling Q and U jointly
                l = 1
                do j = 1, z
                   if (self%temp_corr(n,j)) then
                      do i = 1, x
                         if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                         c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*map2fit(i-1,map_n,j)*&
                              dat%temps(i-1,map_n,n)
                         c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n+1,j)**2.d0)*map2fit(i-1,map_n+1,j)*&
                              dat%temps(i-1,map_n+1,n)
                      end do
                      l = l + 1
                   end if
                end do
                w = w + self%temp_nfit(n)
             end if
          end if
       end do
    end do

    ! write(*,*) c
    
    ! Computation
    if (trim(method) == 'cholesky') then
       write(*,*) 'method: cholesky - Currently deprecated.'
       stop
       !mat_u(:,:)        = 0.d0
       !if (rank == master) write(*,fmt='(a)') 'Joint sampling using Cholesky Decomposition.'
       !call cholesky_decomp(A,mat_l)
       !mat_u  = transpose(mat_l)
       !call forward_sub(mat_l,d,c)
       !call backward_sub(mat_u,b,d)
    else if (trim(method) == 'cg') then
       if (rank == master) write(*,*) 'Joint sampling using CG.'
       call sample_cg_vec(b,c,self,dat,compo,map_n)
    else if (trim(method) == 'lu') then
       write(*,*) 'method: lu - Currently deprecated.'
       stop
       !mat_u(:,:)        = 0.d0
       if (rank == master) write(*,*) 'Joint sampling using LU Decomp'
       !call LUDecomp(A,mat_l,mat_u,y)
       !call forward_sub(mat_l,d,c)
       !call backward_sub(mat_u,b,d)
    end if
    
    ! Solver returns a vector - first filled with component amplitudes, then template amplitudes
    ! So unpack in order
    
    ! Output amplitudes to the appropriate variables
    if (.not. self%joint_pol) then
       w = 0
       do m = 1, size(self%joint_comp)
          do n = 1, self%ncomp
             if (trim(self%joint_comp(m)) == trim(self%fg_label(n))) then
                do i = 1, x
                   dat%fg_map(i-1,map_n,self%fg_ref_loc(n),n) = b(w+i)
                end do
                w = w + x
             end if
          end do
          do n = 1, self%ntemp
             if (trim(self%joint_comp(m)) == trim(self%temp_label(n))) then
                l = 1
                do while (l .lt. self%temp_nfit(n))
                   do j= 1, z
                      if (self%temp_corr(n,j)) then
                         dat%temp_amps(j,map_n,n) = b(w+l)
                         l = l + 1
                      else
                         dat%temp_amps(j,map_n,n) = 0.d0
                      end if
                   end do
                end do
                w = w + l -1
             end if
          end do
       end do
    else if (self%joint_pol) then
       write(*,*) 'output amplitudes'
       w = 0
       do m = 1, size(self%joint_comp)
          do n = 1, self%ncomp
             if (trim(self%joint_comp(m)) == trim(self%fg_label(n))) then
                do i = 1, x
                   dat%fg_map(i-1,map_n,self%fg_ref_loc(n),n) = b(w+i)
                end do
                w = w + x
                do i = 1, x
                   dat%fg_map(i-1,map_n+1,self%fg_ref_loc(n),n) = b(w+i)
                end do
                w = w + x
             end if
          end do
          do n = 1, self%ntemp
             if (trim(self%joint_comp(m)) == trim(self%temp_label(n))) then
                l = 1
                do while (l .lt. self%temp_nfit(n))
                   do j = 1, z
                      if (self%temp_corr(n,j)) then
                         dat%temp_amps(j,map_n,n)   = b(w+l)
                         dat%temp_amps(j,map_n+1,n) = b(w+l)
                         l = l + 1
                      end if
                   end do
                end do
                if (self%temp_nfit(n) == 1) then
                   do j = 1, z
                      if (self%temp_corr(n,j)) then
                         dat%temp_amps(j,map_n,n)   = b(w+l)
                         dat%temp_amps(j,map_n+1,n) = b(w+l)
                         l = l + 1
                      end if
                   end do
                end if
                w = w + l - 1
             end if
          end do
       end do
    end if
    
    ! Also checked and looks good ^^
    
    if (rank == master) then
       t3 = mpi_wtime()
       write(*,fmt='(a,f10.3,a)') 'Joint Sampler completed in ', t3-t1, 's.'
    end if
    
    write(*,*) ''
    
    ! Sure to deallocate all arrays here to free up memory
    deallocate(b)
    deallocate(c)
  end subroutine sample_joint_amp
  
  subroutine template_fit(self, dat, comp, map_n, temp_num)
    !------------------------------------------------------------------------
    ! Simple linear fit of a template to data map with a sampling term
    !------------------------------------------------------------------------
    implicit none
    
    type(params)                                           :: self
    type(component)                                        :: comp
    type(data)                                             :: dat
    !real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: dat, noise
    real(dp), dimension(0:npix-1,nbands)                   :: cov, nos, map
    real(dp), allocatable, dimension(:,:,:)                :: map2fit
    integer(i4b),                               intent(in) :: map_n
    integer(i4b), optional,                     intent(in) :: temp_num
    real(dp)                                               :: temp, sum1, sum2, norm
    integer(i4b)                                           :: i, j, k, n

    nos = dat%rms_map(:,map_n,:)
    cov = nos**2.d0
    allocate(map2fit(0:npix-1,nmaps,nbands))
    
    map2fit = dat%sig_map
    
    if (trim(self%mode) == 'comp_sep') then
       write(*,*) "Sampling for template "//trim(self%temp_label(temp_num))//", pol = "//trim(tqu(map_n))//"."
       do j = 1, self%numinc
          sum1 = 0.d0
          sum2 = 0.d0
          norm = 0.d0
          temp = 0.d0
          if (self%temp_corr(temp_num,j)) then
             ! Remove the other foregrounds from the input data before fitting
             do n = 1, self%ncomp
                ! write(*,*) "remove for foreground ", trim(self%fg_label(n)), j
                map2fit(:,:,j) = map2fit(:,:,j) - dat%fg_map(:,:,j,n)
             end do
             do n = 1, self%ntemp
                if (n == temp_num) cycle
                ! write(*,*) "remove for template ", trim(self%temp_label(n))
                map2fit(:,:,j) = map2fit(:,:,j) - dat%fg_map(:,:,j,self%ncomp+n)
             end do
             ! Calculate template amplitude
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) cycle
                sum1 = sum1 + ((map2fit(i,map_n,j)*dat%temps(i,map_n,temp_num))/cov(i,j))
                sum2 = sum2 + (dat%temps(i,map_n,temp_num)**2.d0)/cov(i,j)
                norm = norm + (dat%temps(i,map_n,temp_num)**2.d0)/cov(i,j)
             end do
             if (trim(self%ml_mode) == 'sample') then
                dat%temp_amps(j,map_n,temp_num) = sum1/sum2 + rand_normal(0.d0,1.d0)/sqrt(norm)
             else if (trim(self%ml_mode) == 'optimize') then
                dat%temp_amps(j,map_n,temp_num) = sum1/sum2
             end if
          end if
       end do
    else if (trim(self%mode) == 'hi_fit') then
       do j = 1, self%numinc
          sum1 = 0.d0
          sum2 = 0.d0
          norm = 0.d0
          temp = 0.d0
          if (self%temp_corr(1,j)) then
             do i = 0, npix-1
                if (comp%HI(i,1) > self%thresh) cycle
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) cycle
                temp = comp%HI(i,1)*planck(self%dat_nu(j)*1d9,comp%T_d(i,1))
                sum1 = sum1 + (((dat%sig_map(i,map_n,j)-dat%offset(j))/dat%gain(j))*temp)/cov(i,j)
                sum2 = sum2 + (temp)**2.d0/cov(i,j)
                norm = norm + (temp)**2.d0/cov(i,j)
             end do
          end if
          if (trim(self%ml_mode) == 'sample') then
             comp%HI_amps(j) = sum1/sum2 + rand_normal(0.d0,1.d0)/sqrt(norm)
          else if (trim(self%ml_mode) == 'optimize') then
             comp%HI_amps(j) = sum1/sum2
             end if
      end do
    end if
    
  end subroutine template_fit
  
  function sample_spec_amp(self, dat, comp, ind, map_n)
    !------------------------------------------------------------------------
    ! Samples spectral amplitude (per pixel), following the spectrum of foreground 'type'. Returns a full map of amplitudes.
    !------------------------------------------------------------------------
    implicit none
    
    class(params)                                          :: self
    type(component)                                        :: comp
    type(data)                                             :: dat
    integer(i4b),                               intent(in) :: ind
    integer(i4b),                               intent(in) :: map_n
    integer(i4b)                                           :: f
    real(dp)                                               :: sum1, sum2, spec
    real(dp)                                               :: amp, num, t, sam
    real(dp), dimension(0:npix-1)                          :: norm
    real(dp), dimension(0:npix-1,nbands)                   :: cov, map2fit
    real(dp), dimension(0:npix-1)                          :: sample_spec_amp
    
    map2fit = dat%sig_map(:,map_n,:)
    
    do f = 1, nfgs
       if (f /= ind) then
          map2fit(:,:) = map2fit(:,:) - dat%fg_map(:,map_n,:,f)
       end if
    end do
    
    cov = dat%rms_map(:,map_n,:)**2.d0
    
    do i = 0, npix-1
       sum1    = 0.0d0
       sum2    = 0.0d0
       norm(i) = 0.d0
       do j = 1, nbands
          spec           = compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n)
          sum1           = sum1 + (map2fit(i,j)*spec)/cov(i,j)
          sum2           = sum2 + (spec)**2.d0/cov(i,j)
          norm(i)        = norm(i) + ((spec)**2.d0)/cov(i,j)
       end do
       amp                = sum1/sum2
       if (trim(self%ml_mode) == 'sample') then
          sample_spec_amp(i) = amp + rand_normal(0.d0,1.d0)/sqrt(norm(i))
       else if (trim(self%ml_mode) == 'optimize') then
          sample_spec_amp(i) = amp
       end if
    end do
  end function sample_spec_amp
  
  subroutine sample_index(self, dat, comp, ind, map_n) 
    !------------------------------------------------------------------------
    ! Warning -- not set up for foregrounds with multiple spectral parameters yet
    ! I.e. will only sample beta_d, and not T_d
    !------------------------------------------------------------------------
    implicit none
    
    class(params)                                          :: self
    type(component),                         intent(inout) :: comp
    type(data)                                             :: dat
    integer(i4b),                               intent(in) :: map_n, ind
    integer(i4b)                                           :: nside1, nside2, npix2, f
    real(dp), dimension(0:npix-1,nmaps,nbands)             :: map2fit, cov 
    real(dp), dimension(0:npix-1,nmaps)                    :: indx
    real(dp), dimension(0:npix-1,nmaps)                    :: mask
    real(dp), dimension(0:npix-1)                          :: indx_sample
    real(dp), allocatable, dimension(:,:,:)                :: fg_map_high
    real(dp), allocatable, dimension(:,:,:)                :: data_low, fg_map_low, cov_low
    real(dp), allocatable, dimension(:,:)                  :: indx_low, mask_low
    real(dp), allocatable, dimension(:)                    :: indx_sample_low
    real(dp)                                               :: a, b, c, num, sam, t, p, sol
    real(dp)                                               :: time1, time2, diff
    real(dp)                                               :: like_old, like_new, ratio, s
    real(dp)                                               :: naccept, paccept, local_a, local_b
    logical                                                :: pix_samp
    logical                                                :: exist
    
    real(dp), allocatable, dimension(:)                    :: beta_grid, like_grid
    
    character(len=128) ::title
    
    real(dp), allocatable, dimension(:,:)                  :: chisq_map
    character(len=3)  :: l_str
    
    !------------------------------------------------------------------------
    ! Spectral index sampler, using the Metropolis approach.
    !------------------------------------------------------------------------    
    allocate(fg_map_high(0:npix-1,nmaps,nbands))

    do i = 0, npix-1
       do j = 1, nbands
          do k = self%pol_type(1), self%pol_type(size(self%pol_type))
             cov(i,k,j)     = dat%rms_map(i,k,j)**2.d0
             map2fit(i,k,j) = dat%sig_map(i,k,j)
          end do
       end do
       mask(i,:) = dat%masks(i,:)
       if (mask(i,self%pol_type(1)) == 0.d0) then
          mask(i,:) = missval
       end if
    end do
    
    do f = 1, nfgs
       if (f /= ind) then
          map2fit(:,:,:) = map2fit(:,:,:) - dat%fg_map(:,:,:,f)
       end if
    end do
    fg_map_high(:,:,:) = dat%fg_map(:,:,:,ind)

    !------------------------------------------------------------------------
    ! Load priors for the appropriate spectrum
    !------------------------------------------------------------------------
    if (trim(self%fg_label(ind)) == 'synch') then 
       write(*,*) "Fitting for synchrotron"
       indx     = comp%beta_s
    else if (trim(self%fg_label(ind)) == 'dust') then 
       write(*,*) "Fitting for thermal dust"
       indx     = comp%beta_d
    end if
    
    if (index(self%fg_ind_region(ind,1),'pix') /= 0) then
       write(*,*) 'single pixel region sampling'
       pix_samp = .true.
       nside2   = self%fg_samp_nside(ind,1)
    else if (index(self%fg_ind_region(ind,1),'full') /= 0) then
       write(*,*) 'Sampling fullsky'
       pix_samp = .false.
    else
       write(*,*) 'ERROR: Sampling region for component '//trim(self%fg_label(ind))//' not recognized.'
       write(*,*) '   Check COMP_BETA_REGION parameter in parameter file.'
       stop
    end if
    
    !-------------------|
    ! Per pixel sampler |
    !-------------------|
    if (pix_samp) then
       if (rank == master) then
          write(*,fmt='(a,i4)') 'Sampling ' // trim(self%fg_label(ind)) // ' beta at nside', nside2
       end if
       
       time1 = mpi_wtime()
       
       !------------------------------------------------------------------------
       ! Check to see if the data nside is the same as the sampling nside
       ! If not equal, downgrade the data before sampling
       !
       ! VITAL: Make dummy varables for each map before ud_grading
       !
       !------------------------------------------------------------------------
       nside1 = npix2nside(npix)
       if (nside1 == nside2) then
          npix2 = npix
       else
          npix2 = nside2npix(nside2)
       end if
       allocate(data_low(0:npix2-1,nmaps,nbands),fg_map_low(0:npix2-1,nmaps,nbands))
       allocate(indx_low(0:npix2-1,nmaps),cov_low(0:npix2-1,nmaps,nbands))
       allocate(indx_sample_low(0:npix2-1))
       allocate(mask_low(0:npix2-1,nmaps))
       
       if (nside1 /= nside2) then 
          if (ordering == 1) then
             call udgrade_ring(mask,nside1,mask_low,nside2,fmissval=missval,PESSIMISTIC=.false.)
             call udgrade_ring(indx,nside1,indx_low,nside2)
             do j = 1, nbands
                call udgrade_ring(map2fit(:,:,j),nside1,data_low(:,:,j),nside2)
                call udgrade_ring(fg_map_high(:,:,j),nside1,fg_map_low(:,:,j),nside2)
                call udgrade_ring(cov(:,:,j),nside1,cov_low(:,:,j),nside2)
             end do
          else
             call udgrade_nest(indx,nside1,indx_low,nside2)
             call udgrade_nest(mask,nside1,mask_low,nside2,fmissval=missval,PESSIMISTIC=.false.)
             do j = 1, nbands
                call udgrade_nest(map2fit(:,:,j),nside1,data_low(:,:,j),nside2)
                call udgrade_nest(dat%fg_map(:,:,j,1),nside1,fg_map_low(:,:,j),nside2)
                call udgrade_nest(cov(:,:,j),nside1,cov_low(:,:,j),nside2)
             end do
          end if
          cov_low = sqrt(cov_low / (npix/npix2))
          do i = 0, npix2-1
             if (mask_low(i,1) .lt. 0.50) then
                mask_low(i,:) = 0.d0
             else
                mask_low(i,:) = 1.d0
             end if
          end do
       else 
          do j = 1, nbands
             data_low(:,:,j)   = map2fit(:,:,j)
             fg_map_low(:,:,j) = dat%fg_map(:,:,j,1)
             cov_low(:,:,j)    = cov(:,:,j)
          end do
          indx_low = indx
          mask_low = mask
       end if
       
       !------------------------------------------------------------------------
       ! Metropolis algorithm:
       ! Sampling portion. Determine the log-likelihood, and accept based off of
       ! the improvement in the fit.
       !------------------------------------------------------------------------
       
       !---------------------------|
       ! Sample for Q and U jointly|
       !---------------------------|
       if (map_n == -1) then
          indx_sample_low = indx_low(:,self%pol_type(1))
          do i = 0, npix2-1
             naccept   = 0.d0
             paccept   = 0.d0
             a         = 0.d0
             sol       = indx_low(i,self%pol_type(1))
             sam       = sol
             if (mask_low(i,1) == 0.d0 .or. mask_low(i,1) == missval) then
                cycle
             end if
             
             ! First evaluate likelihood from previous sample
             !-----------------------------------------------
             
             !$OMP PARALLEL PRIVATE(j,k,local_a)
             local_a = 0.d0
             !$OMP DO SCHEDULE(static)
             do j = 1, nbands
                do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                   local_a = local_a + (((fg_map_low(i,k,self%fg_ref_loc(ind)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,sol)) &
                        - data_low(i,k,j))**2.d0)/cov_low(i,k,j)
                end do
             end do
             !$OMP END DO
             !$OMP CRITICAL
             a = a + local_a
             !$OMP END CRITICAL
             !$OMP END PARALLEL
             c = a
             
             if (self%fg_prior_type(ind,1) == 'gaussian') then
                like_old = -0.5d0*c + log(eval_normal_prior(sam,self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2)))
             else if (self%fg_prior_type(ind,1) == 'uniform') then
                like_old = -0.5d0*c
             else if (self%fg_prior_type(ind,1) == 'jeffreys') then
                like_old = -0.5d0*c + log(eval_jeffreys_prior(self, dat, comp, map_n, ind, sam, i))
             else 
                write(*,*) "error in param%fg_prior_type"
                write(*,*) "      only allowed: 'gaussian', 'uniform', 'jeffreys'"
                stop
             end if

             s = 0.5d0
             do l = 1, self%nsample
                ! Every 50 samples check acceptance rate, and adjust scaling factor
                if (mod(l,50) == 0) then
                   if (paccept > 0.6d0) then
                      s = s*2.0
                   else if (paccept < 0.4d0) then
                      s = s/2.0
                   end if
                end if
                t      = sam + rand_normal(0.d0, self%fg_gauss(ind,1,2))*s
                
                ! If sampled value is outside of uniform bounds, cycle
                if (t .gt. self%fg_uni(ind,1,2) .or. t .lt. self%fg_uni(ind,1,1)) then
                   paccept = naccept/l
                   cycle
                end if
                
                ! Evaluate likelihood given this sample
                !--------------------------------------
                b         = 0.d0
                !$OMP PARALLEL PRIVATE(j,k,local_b)
                local_b   = 0.d0
                !$OMP DO SCHEDULE(static)
                do j = 1, nbands
                   do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                      local_b = local_b + ((fg_map_low(i,k,self%fg_ref_loc(ind))*compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,t) &
                           -data_low(i,k,j))**2.d0)/cov_low(i,k,j)
                   end do
                end do
                !$OMP END DO
                !$OMP CRITICAL
                b = b + local_b
                !$OMP END CRITICAL
                !$OMP END PARALLEL
                
                if (self%fg_prior_type(ind,1) == 'gaussian') then
                   like_new = -0.5d0*b + log(eval_normal_prior(t,self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2)))
                else if (self%fg_prior_type(ind,1) == 'uniform') then
                   like_new = -0.5d0*b
                else if (self%fg_prior_type(ind,1) == 'jeffreys') then
                   like_new = -0.5d0*b + log(eval_jeffreys_prior(self, dat, comp, map_n, ind, t, i))
                else 
                   write(*,*) "error in param%fg_prior_type"
                   write(*,*) "      only allowed: 'gaussian', 'uniform', 'jeffreys'"
                   stop
                end if

                diff = like_new - like_old
                ratio = exp(diff)
                if (trim(self%ml_mode) == 'optimize') then
                   if (ratio > 1.d0) then
                      sam      = t
                      c        = b
                      like_old = like_new
                      naccept  = naccept + 1.0
                   end if
                else if (trim(self%ml_mode) == 'sample') then
                   call RANDOM_NUMBER(num)
                   if (ratio > num) then
                      sam      = t
                      c        = b
                      like_old = like_new
                      naccept  = naccept + 1.0
                   end if
                end if
                
                paccept = naccept/l
                sol = sam
                indx_sample_low(i) = sol
             end do
          end do
          
       !----------------------------|
       ! Sample for a single poltype|
       !----------------------------|
       else
          do i = 0, npix2-1
             naccept   = 0.d0
             paccept   = 0.d0
             a         = 0.d0
             sol       = indx_low(i,self%pol_type(1))
             sam       = sol
             if (mask_low(i,1) == 0.d0 .or. mask_low(i,1) == missval) then
                cycle
             end if

             ! First evaluate likelihood from previous sample
             !-----------------------------------------------
             
             !$OMP PARALLEL PRIVATE(j,k,local_a)
             local_a = 0.d0
             !$OMP DO SCHEDULE(static)
             do j = 1, nbands
                local_a = local_a + (((fg_map_low(i,map_n,self%fg_ref_loc(ind)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n,sol)) &
                     - data_low(i,map_n,j))**2.d0)/cov_low(i,map_n,j)
             end do
             !$OMP END DO
             !$OMP CRITICAL
             a = a + local_a
             !$OMP END CRITICAL
             !$OMP END PARALLEL
             c = a

             if (self%fg_prior_type(ind,1) == 'gaussian') then
                like_old = -0.5d0*c + log(eval_normal_prior(sam,self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2)))
             else if (self%fg_prior_type(ind,1) == 'uniform') then
                like_old = -0.5d0*c
             else if (self%fg_prior_type(ind,1) == 'jeffreys') then
                like_old = -0.5d0*c + log(eval_jeffreys_prior(self, dat, comp, map_n, ind, sam, i))
             else 
                write(*,*) "error in param%fg_prior_type"
                write(*,*) "      only allowed: 'gaussian', 'uniform', 'jeffreys'"
                stop
             end if
             
             s = 0.5d0
             do l = 1, self%nsample
                ! Every 50 samples check acceptance rate, and adjust scaling factor
                if (mod(l,50) == 0) then
                   if (paccept > 0.6d0) then
                      s = s*2.0
                   else if (paccept < 0.4d0) then
                      s = s/2.0
                   end if
                end if
                t      = sam + rand_normal(0.d0, self%fg_gauss(ind,1,2))*s
                
                ! If sampled value is outside of uniform bounds, cycle
                if (t .gt. self%fg_uni(ind,1,2) .or. t .lt. self%fg_uni(ind,1,1)) then
                   paccept = naccept/l
                   cycle
                end if
                
                ! Evaluate likelihood given this sample
                !--------------------------------------
                b         = 0.d0
                !$OMP PARALLEL PRIVATE(j,k,local_b)
                local_b   = 0.d0
                !$OMP DO SCHEDULE(static)
                do j = 1, nbands
                   local_b = local_b + ((fg_map_low(i,map_n,self%fg_ref_loc(ind))*compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n,t) &
                        -data_low(i,map_n,j))**2.d0)/cov_low(i,map_n,j)
                end do
                !$OMP END DO
                !$OMP CRITICAL
                b = b + local_b
                !$OMP END CRITICAL
                !$OMP END PARALLEL
                
                like_new = -0.5d0*b + log(eval_normal_prior(t,self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2)))

                if (self%fg_prior_type(ind,1) == 'gaussian') then
                   like_new = -0.5d0*b + log(eval_normal_prior(t,self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2)))
                else if (self%fg_prior_type(ind,1) == 'uniform') then
                   like_new = -0.5d0*b
                else if (self%fg_prior_type(ind,1) == 'jeffreys') then
                   like_new = -0.5d0*b + log(eval_jeffreys_prior(self, dat, comp, map_n, ind, t, i))
                else 
                   write(*,*) "error in param%fg_prior_type"
                   write(*,*) "      only allowed: 'gaussian', 'uniform', 'jeffreys'"
                   stop
                end if


                diff = like_new - like_old
                ratio = exp(diff)
                if (trim(self%ml_mode) == 'optimize') then
                   if (ratio > 1.d0) then
                      sam      = t
                      c        = b
                      like_old = like_new
                      naccept  = naccept + 1.0
                   end if
                else if (trim(self%ml_mode) == 'sample') then
                   call RANDOM_NUMBER(num)
                   if (ratio > num) then
                      sam      = t
                      c        = b
                      like_old = like_new
                      naccept  = naccept + 1.0
                   end if
                end if
                
                paccept = naccept/l
                sol = sam
                indx_sample_low(i) = sol
             end do
          end do
       end if
       if (nside1 /= nside2) then
          if (ordering == 1) then
             call udgrade_ring(indx_sample_low,nside2,indx_sample,nside1,fmissval=missval)
          else
             call udgrade_nest(indx_sample_low,nside2,indx_sample,nside1,fmissval=missval)
          end if
       else
          indx_sample = indx_sample_low
       end if
      
    !------------------|
    ! Full sky sampler |
    !------------------|
    else
       
       time1 = mpi_wtime()
       
       if (.false.) then
          
          write(*,*) iter
          write(iter_str, '(i0.5)') iter
          
          allocate(beta_grid(100))
          allocate(like_grid(100))
          
          do l = 1, 100
             beta_grid(l) = -3.5 + (-2.0+3.5)*(l-1)/100
             write(*,*) l, beta_grid(l)
             sol       = beta_grid(l)
             !!$OMP PARALLEL PRIVATE(i,j,k,local_a), SHARED(a)
             a         = 0.d0
             local_a   = 0.d0
             !!$OMP DO SCHEDULE(static)
             do i = 0, npix-1
                if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) then
                   cycle
                end if
                
                ! Chi-square from the most recent Gibbs chain update
                do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                   do j = 1, nbands
                      a = a + (((fg_map_high(i,k,self%fg_ref_loc(ind)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,sol)) &
                           - map2fit(i,k,j))**2.d0)/cov(i,k,j)
                      ! local_a = local_a + (((fg_map_high(i,k,self%fg_ref_loc(ind)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,sol)) &
                      !      - map2fit(i,k,j))**2.d0)/cov(i,k,j)
                   end do
                end do
             end do
             !!$OMP END DO
             !!$OMP CRITICAL
             ! a = a + local_a
             !!$OMP END CRITICAL
             !!$OMP END PARALLEL
             
             title = trim(self%outdir) // 'beta_grid_likelihood_k'//trim(iter_str)//'.dat'
             inquire(file=title,exist=exist)
             if (exist) then
                open(12,file=title, status="old",position="append", action="write")
             else
                open(12,file=title, status="new", action="write")
                write(12,*)  
             endif
             write(12,*) beta_grid(l), -0.5d0*a
             close(12)
          end do
       end if
       indx_sample(:) = -3.1d0
       
       !----------------------------|
       ! Sample for Q and U jointly |
       !----------------------------|
       if (map_n == -1) then
          naccept     = 0.d0
          paccept     = 0.d0
          chisq_map   = 0.d0
          indx_sample = indx(:,self%pol_type(1))
          sol         = indx(0,self%pol_type(1))
          sam         = sol
          
          ! First evaluate likelihood from previous sample
          !-----------------------------------------------
          
          ! !$OMP PARALLEL PRIVATE(i,j,k,local_a), SHARED(a)
          ! a           = 0.d0
          ! local_a     = 0.d0
          ! !$OMP DO SCHEDULE(static)
          ! do i = 0, npix-1
          !    if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
          !    do j = 1, nbands
          !       do k = self%pol_type(1), self%pol_type(size(self%pol_type))
          !          local_a = local_a + (((fg_map_high(i,k,self%fg_ref_loc(1)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,sol)) &
          !               - map2fit(i,k,j))**2.d0)/cov(i,k,j)
          !       end do
          !    end do
          ! end do
          ! !$OMP END DO
          ! !$OMP CRITICAL
          ! a = a + local_a
          ! !$OMP END CRITICAL
          ! !$OMP END PARALLEL
          ! c = a

          a           = 0.d0
          do i = 0, npix-1
             if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
             do j = 1, nbands
                do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                   sam = indx(i,k)
                   a   = a + (((fg_map_high(i,k,self%fg_ref_loc(ind))*compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,sam)) &
                        - map2fit(i,k,j))**2.d0)/cov(i,k,j)
                end do
             end do
          end do
          c = a

          if (self%fg_prior_type(ind,1) == 'gaussian') then
             like_old = -0.5d0*c + log(eval_normal_prior(sam,self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2)))
          else if (self%fg_prior_type(ind,1) == 'uniform') then
             like_old = -0.5d0*c
          else if (self%fg_prior_type(ind,1) == 'jeffreys') then
             like_old = -0.5d0*c + log(eval_jeffreys_prior(self, dat, comp, map_n, ind, sam))
          else 
             write(*,*) "error in param%fg_prior_type"
             write(*,*) "      only allowed: 'gaussian', 'uniform', 'jeffreys'"
             stop
          end if

          s = 0.5d0
          do l = 1, self%nsample
             ! Every 50 samples check acceptance rate, and adjust scaling factor
             if (mod(l,50) == 0) then
                if (paccept > 0.6d0) then
                   write(*,*) paccept,l,'s = s*2.0'
                   s = s*2.0
                else if (paccept < 0.4d0) then
                   write(*,*) paccept,l,'s = s/2.0'
                   s = s/2.0
                end if
             end if
             t      = sam + rand_normal(0.d0, self%fg_gauss(ind,1,2))*s
             
             ! If sampled value is outside of uniform bounds, cycle
             if (t .gt. self%fg_uni(ind,1,2) .or. t .lt. self%fg_uni(ind,1,1)) then
                paccept = naccept/l
                cycle
             end if
             ! ! Evaluate likelihood given this sample
             ! !--------------------------------------
             ! b         = 0.d0
             ! !$OMP PARALLEL PRIVATE(i,j,k,local_b)
             ! local_b   = 0.d0
             ! !$OMP DO SCHEDULE(static)
             ! do i = 0, npix-1
             !    if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
             !    do j = 1, nbands
             !       do k = self%pol_type(1), self%pol_type(size(self%pol_type))
             !          local_b = local_b + (((fg_map_high(i,k,self%fg_ref_loc(ind)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,t)) &
             !               - map2fit(i,k,j))**2.d0)/cov(i,k,j)
             !       end do
             !    end do
             ! end do
             ! !$OMP END DO
             ! !$OMP CRITICAL
             ! b = b + local_b
             ! !$OMP END CRITICAL
             ! !$OMP END PARALLEL

             b         = 0.d0
             do i = 0, npix-1
                if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
                do j = 1, nbands
                   do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                      b = b + (((fg_map_high(i,k,self%fg_ref_loc(1)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,t)) &
                           - map2fit(i,k,j))**2.d0)/cov(i,k,j)
                   end do
                end do
             end do

             if (self%fg_prior_type(ind,1) == 'gaussian') then
                like_new = -0.5d0*b + log(eval_normal_prior(t,self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2)))
             else if (self%fg_prior_type(ind,1) == 'uniform') then
                like_new = -0.5d0*b
             else if (self%fg_prior_type(ind,1) == 'jeffreys') then
                like_new = -0.5d0*b + log(eval_jeffreys_prior(self, dat, comp, map_n, ind, t))
             else 
                write(*,*) "error in param%fg_prior_type"
                write(*,*) "      only allowed: 'gaussian', 'uniform', 'jeffreys'"
                stop
             end if

             diff = like_new - like_old
             ratio = exp(diff)
             if (trim(self%ml_mode) == 'optimize') then
                if (ratio > 1.d0) then
                   sam      = t
                   c        = b
                   like_old = like_new
                   naccept  = naccept + 1.0
                end if
             else if (trim(self%ml_mode) == 'sample') then
                call RANDOM_NUMBER(num)
                if (ratio > num) then
                   sam      = t
                   c        = b
                   like_old = like_new
                   naccept  = naccept + 1.0
                end if
             end if
             
             paccept = naccept/l
             ! Adjust the step-size based on acceptance probability and document
          end do
          title = trim(self%outdir) // 'synch_beta_paccept.dat'
          inquire(file=title,exist=exist)
          if (exist) then
             open(12,file=title, status="old",position="append", action="write")
          else
             open(12,file=title, status="new", action="write")
             write(12,*)  
          endif
          write(12,*) paccept, s
          close(12)
          sol = sam
          indx_sample(:) = sol
          
       !---------------------------------------|
       ! Sample for a single poltype - fullsky |
       !---------------------------------------|
       else

          sol         = indx(0,map_n)
          sam         = sol
          naccept     = 0.d0
          paccept     = 0.d0

          a           = 0.d0
          do i = 0, npix-1
             if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle

             do j = 1, nbands
                a = a + (((fg_map_high(i,map_n,self%fg_ref_loc(ind)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n,sam)) &
                     - map2fit(i,map_n,j))**2.d0)/cov(i,map_n,j)
             end do
          end do

          ! ! First evaluate likelihood from previous sample
          ! !-----------------------------------------------
          ! !$OMP PARALLEL PRIVATE(i,j,k,local_a), SHARED(a)
          ! a         = 0.d0
          ! local_a   = 0.d0
          ! !$OMP DO SCHEDULE(static)
          ! do i = 0, npix-1
          !    if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) then
          !       cycle
          !    end if
             
          !    do j = 1, nbands
          !       local_a = local_a + (((fg_map_high(i,map_n,self%fg_ref_loc(1)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n,sol)) &
          !            - map2fit(i,map_n,j))**2.d0)/cov(i,map_n,j)
          !    end do
          ! end do
          ! !$OMP END DO
          ! !$OMP CRITICAL
          ! a = a + local_a
          ! !$OMP END CRITICAL
          ! !$OMP END PARALLEL

          c = a
          
          if (self%fg_prior_type(ind,1) == 'gaussian') then
             like_old = -0.5d0*c + log(eval_normal_prior(sam,self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2)))
          else if (self%fg_prior_type(ind,1) == 'uniform') then
             like_old = -0.5d0*c
          else if (self%fg_prior_type(ind,1) == 'jeffreys') then
             like_old = -0.5d0*c + log(eval_jeffreys_prior(self, dat, comp, map_n, ind, sam))
          else 
             write(*,*) "error in param%fg_prior_type"
             write(*,*) "      only allowed: 'gaussian', 'uniform', 'jeffreys'"
             stop
          end if
          
          s = 0.5d0
          do l = 1, self%nsample
             ! Every 50 samples check acceptance rate, and adjust scaling factor
             if (mod(l,50) == 0) then
                if (paccept > 0.6d0) then
                   write(*,*) paccept,l,'s = s*2.0'
                   s = s*2.0
                else if (paccept < 0.4d0) then
                   write(*,*) paccept,l,'s = s/2.0'
                   s = s/2.0
                end if
             end if
             t      = sam + rand_normal(0.d0, self%fg_gauss(ind,1,2))*s
             
             ! If sampled value is outside of uniform bounds, cycle
             if (t .gt. self%fg_uni(ind,1,2) .or. t .lt. self%fg_uni(ind,1,1)) then
                paccept = naccept/l
                cycle
             end if
             ! Evaluate likelihood given this sample
             !--------------------------------------
             b         = 0.d0

             do i = 0, npix-1
                if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
                do j = 1, nbands
                   b = b + (((fg_map_high(i,map_n,self%fg_ref_loc(ind)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n,t)) &
                        - map2fit(i,map_n,j))**2.d0)/cov(i,map_n,j)
                end do
             end do

             ! !$OMP PARALLEL PRIVATE(i,j,k,local_b)
             ! local_b   = 0.d0
             ! !$OMP DO SCHEDULE(static)
             ! do i = 0, npix-1
             !    if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
             !    do j = 1, nbands
             !       local_b = local_b + (((fg_map_high(i,map_n,self%fg_ref_loc(1)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n,t)) &
             !            - map2fit(i,map_n,j))**2.d0)/cov(i,map_n,j)
             !    end do
             ! end do
             ! !$OMP END DO
             ! !$OMP CRITICAL
             ! b = b + local_b
             ! !$OMP END CRITICAL
             ! !$OMP END PARALLEL
             
             if (self%fg_prior_type(ind,1) == 'gaussian') then
                like_new = -0.5d0*b + log(eval_normal_prior(t,self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2)))
             else if (self%fg_prior_type(ind,1) == 'uniform') then
                like_new = -0.5d0*b
             else if (self%fg_prior_type(ind,1) == 'jeffreys') then
                like_new = -0.5d0*b + log(eval_jeffreys_prior(self, dat, comp, map_n, ind, t))
             else 
                write(*,*) "error in param%fg_prior_type"
                write(*,*) "      only allowed: 'gaussian', 'uniform', 'jeffreys'"
                stop
             end if
             
             diff = like_new - like_old
             ratio = exp(diff)
             call RANDOM_NUMBER(num)
             if (trim(self%ml_mode) == 'optimize') then
                if (ratio > 1.d0) then
                   sam      = t
                   c        = b
                   like_old = like_new
                   naccept  = naccept + 1.0
                end if
             else if (trim(self%ml_mode) == 'sample') then
                call RANDOM_NUMBER(num)
                if (ratio > num) then
                   sam      = t
                   c        = b
                   like_old = like_new
                   naccept  = naccept + 1.0
                end if
             end if
             
             paccept = naccept/l
             ! Adjust the step-size based on acceptance probability and document
          end do
          sol = sam
          indx_sample(:) = sol
       end if
    end if
    
    if (map_n == -1) then
       do k = self%pol_type(1), self%pol_type(size(self%pol_type))
          if (trim(self%fg_label(ind)) == 'synch') then 
             comp%beta_s(:,k) = indx_sample(:)
          else if (trim(self%fg_label(ind)) == 'dust') then 
             comp%beta_d(:,k) = indx_sample(:)
          end if
       end do
    else
       if (trim(self%fg_label(ind)) == 'synch') then 
          comp%beta_s(:,map_n) = indx_sample(:)
       else if (trim(self%fg_label(ind)) == 'dust') then 
          comp%beta_d(:,map_n) = indx_sample(:)
       end if
    end if
    
    time2 = mpi_wtime()
    
    write(*,*) ''
    if (rank == master) then
       time2 = mpi_wtime()
       write(*,fmt='(a,f10.3,a)') 'Spectral index sampler completed in ', time2-time1, 's.'
    end if
    
    write(*,*) ''
    
  end subroutine sample_index
  
  ! This architecture of this function has not been verified yet
  subroutine sample_HI_T(self, dat, comp, map_n)
    implicit none
    
    class(params)                                          :: self
    type(data)                                             :: dat
    type(component)                                        :: comp
    integer(i4b), intent(in)                               :: map_n
    integer(i4b)                                           :: nside1, npix2, nside2
    real(dp), dimension(0:npix-1,nmaps,nbands)             :: cov
    real(dp), dimension(0:npix-1,nmaps)                    :: te
    real(dp), dimension(0:npix-1)                          :: te_sample
    real(dp), allocatable, dimension(:,:,:)                :: maps_low, cov_low
    real(dp), allocatable, dimension(:,:)                  :: T_low
    real(dp), allocatable, dimension(:)                    :: sample_T_low
    real(dp)                                               :: a, b, c, num, sam
    real(dp)                                               :: t, p, sol
    real(dp)                                               :: like_old, like_new
    real(dp)                                               :: diff, ratio
    real(dp)                                               :: paccept, naccept, s
    
    te      = comp%T_d
    cov     = dat%rms_map*dat%rms_map
    
    sam = 0.0
    
    nside2 = 16!64
    
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
    !             call convert_nest2ring(nside1,dat%rms_map(:,:,j))
    !         else
    !             call udgrade_nest(map2fit(:,:,j),nside1,band_low(:,:,j),nside2)
    !             call udgrade_nest(dat%rms_map(:,:,j),nside1,cov_low(:,:,j),nside2)
    !         end if
    !     end do
    !     cov_low = sqrt(cov_low / (npix/npix2))
    ! else
    do j = 1, nbands
       maps_low(:,:,j)   = dat%sig_map(:,:,j)
       cov_low(:,:,j)    = cov(:,:,j)
    end do
    T_low = te
    ! end if
    
    ! Metropolis algorithm
    
    do i = 0, npix2-1
       naccept = 0.d0
       paccept = 0.d0
       s       = 0.d0
       if (dat%masks(i,1) == 0.d0  .or. dat%masks(i,1) == missval) then
          sample_T_low(i) = missval
          cycle
       else
          a   = 0.d0
          sol = T_low(i,map_n)
          sam = sol
          ! Chi-square from the most recent Gibbs chain update
          do j = 1, nbands
             a = a + (((comp%HI_amps(j)*comp%HI(i,1)*planck(self%dat_nu(j)*1.d9,sol)) &
                  - (maps_low(i,map_n,j)-dat%offset(j))/dat%gain(j))**2.d0)/cov_low(i,map_n,j)
          end do
          c   = a
          
          like_old = -0.5d0*c + log(eval_normal_prior(sam,self%HI_Td_mean, self%HI_Td_std))
          
          do l = 1, self%nsample
             if (mod(l,50) == 0) then
                if (paccept > 0.6d0) then
                   s = s*2.0
                else if (paccept < 0.4d0) then
                   s = s/2.0
                end if
             end if
             ! Begin sampling from the prior
             t = rand_normal(sol,self%HI_Td_std)
             b = 0.d0
             do j = 1, nbands
                b = b + (((comp%HI_amps(j)*comp%HI(i,1)*planck(self%dat_nu(j)*1.d9,t)) &
                     - (maps_low(i,map_n,j)-dat%offset(j))/dat%gain(j))**2.d0)/cov_low(i,map_n,j)
             end do
             b = b
             
             like_new = -0.5d0*b + log(eval_normal_prior(sam,self%HI_Td_mean, self%HI_Td_std))
             diff = like_new - like_old
             ratio = exp(diff)
             if (trim(self%ml_mode) == 'optimize') then
                if (ratio > 1.d0) then
                   sam      = t
                   c        = b
                   like_old = like_new
                   naccept  = naccept + 1.0
                end if
             else if (trim(self%ml_mode) == 'sample') then
                call RANDOM_NUMBER(num)
                if (ratio > num) then
                   sam      = t
                   c        = b
                   like_old = like_new
                   naccept  = naccept + 1.0
                end if
             end if
             paccept = naccept/l
          end do
          if (sam == 0.d0) then
             write(*,*) 'Error: T_d = 0.d0 accepted!'
             stop
          end if
          
          sol             = sam
          sample_T_low(i) = sol
       end if
    end do
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
    comp%T_d(:,1) = te_sample
    
    deallocate(maps_low)
    deallocate(T_low)
    deallocate(cov_low)
    deallocate(sample_T_low)
    
  end subroutine sample_HI_T
  ! ------------------------------------------------------------
  
  subroutine sample_band_gain(self, dat, comp, map_n, band, fg, sample)
    class(params)                                          :: self
    type(data)                                             :: dat
    type(component)                                        :: comp
    integer(i4b),                               intent(in) :: map_n
    integer(i4b),                               intent(in) :: band
    integer(i4b),                               intent(in) :: fg
    integer(i4b), optional,                     intent(in) :: sample
    real(dp), allocatable, dimension(:)                    :: map1, map2, mask, noise, N_inv
    real(dp)                                               :: norm, gain
    
    allocate(map1(0:dat%npix-1))
    allocate(map2(0:dat%npix-1))
    allocate(mask(0:dat%npix-1))
    allocate(noise(0:dat%npix-1))
    allocate(N_inv(0:dat%npix-1))
    
    map1 = 0.d0
    map2 = 0.d0
    mask = 0.d0
    
    mask = dat%masks(:,1)
    
    ! map1 is the map we calibrate against here, being a component map.
    ! Must ensure the map is calculated prior to gain fitting
    ! for the HI fit, we use the foreground model to fit
    
    if (trim(self%mode) == 'hi_fit') then
       do i = 0, dat%npix-1
          if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
          map1(i) = comp%HI_amps(band)*comp%HI(i,1)*planck(self%dat_nu(band)*1d9,comp%T_d(i,1))
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
    
    norm = sqrt(sum(mask*map1*N_inv*map1))
    
    ! Sample variable can be any number. If it's present, we sample!
    if (present(sample)) then
       gain = gain + rand_normal(0.d0,1.d0)/norm
    end if
    
    ! Save to data type variable corresponding to the band.
    dat%gain(band) = gain
    
  end subroutine sample_band_gain
  
  subroutine sample_band_offset(self, dat, comp, map_n, band, fg, sample)
    class(params)                                          :: self
    type(data)                                             :: dat
    type(component)                                        :: comp
    integer(i4b),                               intent(in) :: map_n
    integer(i4b),                               intent(in) :: band
    integer(i4b),                               intent(in) :: fg
    integer(i4b), optional,                     intent(in) :: sample
    real(dp), allocatable, dimension(:)                    :: map1, map2, mask
    real(dp)                                               :: norm, offset
    real(dp)                                               :: n, x, x2, y, y2, xy
    
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
    
    ! map1 is the map we calibrate against here, being a component map.
    ! Must ensure the map is calculated prior to gain fitting
    ! for the HI fit, we use the foreground model to fit (amplitude*HI*B_nu(T))
    
    if (trim(self%mode) == 'hi_fit') then
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
    
    ! stop
    
    ! norm   = sum(mask(:))*sum(mask(:)*(map1(:)**2))-(sum(mask(:)*map1(:))**2)
    
    ! if (present(sample)) then
    !    offset = offset + rand_normal(0.d0,1.d0)/norm
    ! end if
    
    dat%offset(band) = offset
    
  end subroutine sample_band_offset
  
  subroutine calc_hi_gain_offset(self, dat, comp, map_n, fg, band)
    class(params)                                          :: self
    type(data)                                             :: dat
    type(component)                                        :: comp
    integer(i4b),                               intent(in) :: map_n
    integer(i4b),                               intent(in) :: band
    integer(i4b),                               intent(in) :: fg
    real(dp), allocatable, dimension(:)                    :: map1, map2, mask
    real(dp), allocatable, dimension(:)                    :: gain, offset
    
    
    allocate(map1(0:dat%npix-1))
    allocate(map2(0:dat%npix-1))
    allocate(mask(0:dat%npix-1))
    
    allocate(gain(0:self%nsample))
    allocate(offset(0:self%nsample))
    
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
    
    if (trim(self%mode) == 'hi_fit') then
       do i = 0, dat%npix-1
          if (mask(i) == 0.d0 .or. mask(i) == missval) cycle
          map1(i) = comp%HI_amps(band)*comp%HI(i,1)*planck(self%dat_nu(band)*1d9,comp%T_d(i,1))
       end do
    else
       map1 = dat%fg_map(:,map_n,band,fg)
    end if
    
    map2 = dat%sig_map(:,map_n,band)
    
    do i = 1, self%nsample-1
       
       ! Fit gain to the SED first (HI in the HI fit case)
       
       ! if (trim(self%mode) == 'hi_fit') then
       !    offset(i) = sum(mask(:)*(map2(:) - gain(i-1)*comp%HI(:,1)))/sum(mask(:))
       ! else
       offset(i) = sum(mask(:)*(map2(:) - gain(i-1)*map1(:)))/sum(mask(:))
       ! end if                 
       ! Calculate the offset using the HI map from the calibrated band
       
       gain(i)   = sum(mask(:)*(map1(:)*(map2(:)-offset(i))))/sum(mask(:)*(map1(:)**2))
       
       if (i > self%nsample) then
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

  function eval_jeffreys_prior(self, dat, comp, map_n, ind, val, pixel) result(prob)
    implicit none

    class(params)                                          :: self
    type(component),                         intent(inout) :: comp
    type(data)                                             :: dat
    real(dp),                                   intent(in) :: val
    integer(i4b),                               intent(in) :: map_n, ind
    integer(i4b), optional,                     intent(in) :: pixel
    real(dp)                                               :: prob, sum, ss

    prob = 0.d0
    sum = 0.d0

    if (trim(self%fg_label(ind)) == 'synch') then
       ! Is this evaluated for a single pixel?
       if (present(pixel)) then
          ! If map_n = -1, then sum over poltypes
          if (map_n == -1) then
             do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                do j = 1, nbands
                   ss  = dat%fg_map(pixel,k,self%fg_ref_loc(ind),ind)*(self%dat_nu(j)/self%fg_nu_ref(ind))**val
                   sum = sum + (((1.0/dat%rms_map(pixel,k,j))**2)*(ss/dat%fg_map(pixel,k,self%fg_ref_loc(ind),ind))*log(self%dat_nu(j)/self%fg_nu_ref(ind)))**2.0
                end do
             end do
          else
             do j = 1, nbands
                ss  = dat%fg_map(pixel,map_n,self%fg_ref_loc(ind),ind)*(self%dat_nu(j)/self%fg_nu_ref(ind))**val
                sum = sum + (((1.0/dat%rms_map(pixel,map_n,j))**2)*(ss/dat%fg_map(pixel,map_n,self%fg_ref_loc(ind),ind))*log(self%dat_nu(j)/self%fg_nu_ref(ind)))**2.0
             end do
          end if
       else
          ! If map_n = -1, then sum over poltypes
          if (map_n == -1) then
             do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                do i = 0, npix-1
                   if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) cycle
                   do j = 1, nbands
                      ss  = dat%fg_map(i,k,self%fg_ref_loc(ind),ind)*(self%dat_nu(j)/self%fg_nu_ref(ind))**val
                      sum = sum + (((1.0/dat%rms_map(i,k,j))**2)*(ss/dat%fg_map(i,k,self%fg_ref_loc(ind),ind))*log(self%dat_nu(j)/self%fg_nu_ref(ind)))**2.0
                   end do
                end do
             end do
          else
             do i = 0, npix-1
                if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) cycle
                do j = 1, nbands
                   ss  = dat%fg_map(i,k,self%fg_ref_loc(ind),ind)*(self%dat_nu(j)/self%fg_nu_ref(ind))**val
                   sum = sum + (((1.0/dat%rms_map(i,k,j))**2)*(ss/dat%fg_map(i,k,self%fg_ref_loc(ind),ind))*log(self%dat_nu(j)/self%fg_nu_ref(ind)))**2.0
                end do
             end do
          end if
       end if
    end if
    prob = sqrt(sum)

  end function eval_jeffreys_prior
  
end module dang_sample_mod
