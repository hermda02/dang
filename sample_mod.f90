module sample_mod
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
    use utility_mod
    use dang_param_mod
    use linalg_mod
    use dang_component_mod
    use dang_data_mod
    implicit none

    private :: i, j, k
    integer(i4b) :: i, j ,k
    
    logical(lgt) :: exist
    

contains

    !function temp_fit(data,template,noise,t_mask,freq)
    !    implicit none
    
    !    real(dp), dimension(0:npix-1), intent(in) :: data, template, noise, t_mask
    !    real(dp), optional, intent(in)            :: freq
    !    real(dp), dimension(0:npix-1)             :: cov
    !    real(dp)                                  :: temp_fit, norm, sam, p
    !    real(dp)                                  :: amp, old, sum1, sum2, num
    
    !    cov = noise*2

        ! Uncertainty, used for sampling.
    !    norm = 0.d0
    !    sum1 = 0.d0
    !    sum2 = 0.d0
    
    !    if (present(freq)) then
    !        do i=0,npix-1
    !            if (HI(i,1) > par%thresh) cycle
    !            p      = template(i)*planck(freq*1.d9,T_d(i,1))
    !            sum1   = sum1 + (data(i)*p)/cov(i)*t_mask(i)
    !            sum2   = sum2 + (p)**2.d0/cov(i)*t_mask(i)
    !            norm   = norm + (p)**2.d0/cov(i)*t_mask(i)
    !        end do
    !        norm = norm/n2fit
    !    else
    !        norm  = sum((template(:)**2.d0)/cov(:))/sum(t_mask)
    !        do i=0,npix-1
    !            sum1   = sum1 + (data(i)*template(i))/cov(i)*t_mask(i)
    !            sum2   = sum2 + template(i)**2.d0/cov(i)*t_mask(i)
    !        end do
    !    end if

    !    ! Don't allow negative amplitudes, if negative, set to 0.d0.
    !    if (sum1 < 0.d0) then
    !        amp = 0.d0
    !    else
    !        amp   = sum1/sum2! + rand_normal(0.d0,1.d0)/sqrt(norm)
    !    end if
    !    temp_fit  = amp
  
    !end function temp_fit

   function sample_spec_amp(self, comp, data, noise, ind, mapn)
        !------------------------------------------------------------------------
        ! Samples spectral amplitude (per pixel), following the spectrum of foreground 'type'. Returns a full map of amplitudes.
        !------------------------------------------------------------------------
        implicit none
  
        class(params)                                          :: self
        type(component)                                        :: comp
        real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: data, noise
        integer(i4b),                               intent(in) :: ind
        integer(i4b),                               intent(in) :: mapn
        real(dp)                                               :: sum1, sum2, spec
        real(dp)                                               :: chi, chi_0, chi_00, p
        real(dp)                                               :: amp, num, t, sam
        real(dp), dimension(2)                                 :: pars
        real(dp), dimension(nbands)                            :: tmp
        real(dp), dimension(0:npix-1)                          :: norm
        real(dp), dimension(0:npix-1,nbands)                   :: cov, nos
        real(dp), dimension(0:npix-1)                          :: sample_spec_amp

        nos = noise(:,mapn,:)
        cov = nos**2.d0

        do i = 0, npix-1
            sum1    = 0.0d0
            sum2    = 0.0d0
            norm(i) = 0.d0
            do j = 1, nbands
                spec           = compute_spectrum(self,comp,ind,self%dat_nu(j),i,mapn)
                sum1           = sum1 + (data(i,mapn,j)*spec)/cov(i,j)
                sum2           = sum2 + (spec)**2.d0/cov(i,j)
                norm(i)        = norm(i) + ((spec)**2.d0)/cov(i,j)
            end do
            norm(i)            = norm(i)
            amp                = sum1/sum2
            sample_spec_amp(i) = amp + rand_normal(0.d0,1.d0)/sqrt(norm(i))
        end do

    end function sample_spec_amp

    subroutine sample_index(self, comp, data, nside2, ind, map_n)
        implicit none
  
        class(params)                                          :: self
        type(component),                         intent(inout) :: comp
        integer(i4b),                               intent(in) :: map_n, nside2, ind
        real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: data
        integer(i4b)                                           :: nside1, npix2
        real(dp), dimension(0:npix-1,nmaps,nbands)             :: map2fit, cov 
        real(dp), dimension(0:npix-1,nmaps)                    :: indx
        real(dp), dimension(0:npix-1)                          :: indx_sample
        real(dp), allocatable, dimension(:,:,:)                :: data_low, fg_map_low, rms_low
        real(dp), allocatable, dimension(:,:)                  :: indx_low
        real(dp), allocatable, dimension(:)                    :: indx_sample_low
        real(dp), dimension(nbands)                            :: signal, tmp
        real(dp), dimension(2)                                 :: x
        real(dp)                                               :: a, b, c, num, sam, t, p, sol

        real(dp)                                               :: naccept   
        logical                                                :: exist

        !------------------------------------------------------------------------
        ! Spectral index sampler, using the Metropolis approach.
        !------------------------------------------------------------------------

        map2fit = data
        cov     = dang_data%rms_map*dang_data%rms_map

        !------------------------------------------------------------------------
        ! Load priors for the appropriate spectrum
        !------------------------------------------------------------------------
        if (trim(self%fg_label(ind)) == 'synch') then 
            indx     = comp%beta_s
        else if (trim(self%fg_label(ind)) == 'dust') then 
            indx     = comp%beta_d
        end if
        
        if (rank == master) then
            if (mod(iter,output_iter) .EQ. 0) then
                write(*,fmt='(a,i4)') 'Sampling ' // trim(self%fg_label(ind)) // ' beta at nside', nside2
            end if
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
        allocate(data_low(0:npix2-1,nmaps,nbands),fg_map_low(0:npix2-1,nmaps,nbands))
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
                    call udgrade_ring(dang_data%fg_map(:,:,j,1),nside1,fg_map_low(:,:,j),nside2)
                    call convert_nest2ring(nside1,dang_data%fg_map(:,:,j,1))
                    call udgrade_ring(cov(:,:,j),nside1,rms_low(:,:,j),nside2)
                    call convert_nest2ring(nside2,rms_low(:,:,j))
                else
                    call udgrade_nest(map2fit(:,:,j),nside1,data_low(:,:,j),nside2)
                    call udgrade_nest(dang_data%fg_map(:,:,j,1),nside1,fg_map_low(:,:,j),nside2)
                    call udgrade_nest(dang_data%rms_map(:,:,j),nside1,rms_low(:,:,j),nside2)
               end if
            end do
            rms_low = sqrt(rms_low / (npix/npix2))
        else 
            do j = 1, nbands
                data_low(:,:,j)   = data(:,:,j)
                fg_map_low(:,:,j) = dang_data%fg_map(:,:,j,1)
                rms_low(:,:,j)    = dang_data%rms_map(:,:,j)
            end do
            indx_low = indx
        end if

        x(1) = 1.d0           
        !------------------------------------------------------------------------
        ! Sampling portion. Determine the log-likelihood, and accept based off of
        ! the improvement in the fit.
        !------------------------------------------------------------------------
        if (map_n == -1) then

           do i = 0, npix2-1
              a         = 0.d0
              sol       = indx_low(i,self%pol_type(1))
              
              ! Chi-square from the most recent Gibbs chain update
              do j = 1, nbands
                 do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                    a = a + (((fg_map_low(i,k,par%fg_ref_loc(1)) * compute_spectrum(self,comp,ind,par%dat_nu(j),i,k,sol)) &
                         - data_low(i,k,j))**2.d0)/rms_low(i,k,j)**2.d0
                 end do
              end do
              c = a

              do l = 1, iterations

                 ! Sampling from the prior
                 t         = rand_normal(sol, self%fg_gauss(ind,1,2))
                 !t         = rand_normal(self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2))
                 b         = 0.d0
                 
                 do j = 1, nbands
                    do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                       tmp(j) = fg_map_low(i,k,par%fg_ref_loc(1))*compute_spectrum(self,comp,ind,par%dat_nu(j),i,k,t)
                       b      = b + ((tmp(j)-data_low(i,k,j))**2.d0)/rms_low(i,k,j)**2.d0
                    end do
                 end do
                 b = b

                 if (b < c .and. t .lt. self%fg_uni(ind,1,2) .and. t .gt. self%fg_uni(ind,1,1)) then
                    sam = t
                    c   = b
                 else
                    x(2) = exp(0.5d0*(c-b))
                    p = minval(x)
                    call RANDOM_NUMBER(num)
                    if (num < p) then
                       if (t .lt. self%fg_uni(ind,1,2) .and. t .gt. self%fg_uni(ind,1,1)) then
                          sam = t
                          c   = b
                       end if
                    end if
                 end if
              end do
              sol = sam
              indx_sample_low(i) = sol
           end do
        else
           do i = 0, npix2-1
              a         = 0.d0
              sol       = indx_low(i,map_n)
              
              ! Chi-square from the most recent Gibbs chain update
              do j = 1, nbands
                 a = a + (((fg_map_low(i,map_n,par%fg_ref_loc(1)) * compute_spectrum(self,comp,ind,par%dat_nu(j),i,map_n,sol)) &
                      - data_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
              end do
              c = a
              
              do l = 1, iterations
                 
                 ! Sampling from the prior
                 t         = rand_normal(sol, self%fg_gauss(ind,1,2))
                 !t         = rand_normal(self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2))
                 b         = 0.d0
                 
                 do j = 1, nbands
                    tmp(j) = fg_map_low(i,map_n,par%fg_ref_loc(1))*compute_spectrum(self,comp,ind,par%dat_nu(j),i,map_n,t)
                    b      = b + ((tmp(j)-data_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
                 end do
                 b = b
                 
                 if (b < c .and. t .lt. self%fg_uni(ind,1,2) .and. t .gt. self%fg_uni(ind,1,1)) then
                    sam = t
                    c   = b
                 else
                    x(2) = exp(0.5d0*(c-b))
                    p = minval(x)
                    call RANDOM_NUMBER(num)
                    if (num < p) then
                       if (t .lt. self%fg_uni(ind,1,2) .and. t .gt. self%fg_uni(ind,1,1)) then
                          sam = t
                          c   = b
                       end if
                    end if
                 end if
              end do
              sol = sam
              indx_sample_low(i) = sol
           end do
        end if

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

        if (map_n == -1) then
           do k = self%pol_type(1), self%pol_type(size(self%pol_type))
              if (trim(self%fg_label(ind)) == 'synch') then 
                 beta_s(:,k) = indx_sample
              else if (trim(self%fg_label(ind)) == 'dust') then 
                 beta_d(:,k) = indx_sample
              end if
           end do
        else
           if (trim(self%fg_label(ind)) == 'synch') then 
              comp%beta_s(:,k) = indx_sample
           else if (trim(self%fg_label(ind)) == 'dust') then 
              comp%beta_d(:,k) = indx_sample
           end if
        end if
        deallocate(data_low)
        deallocate(fg_map_low)
        deallocate(indx_low)
        deallocate(rms_low)
        deallocate(indx_sample_low)

    end subroutine sample_index
 

end module sample_mod
