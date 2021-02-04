module dang_sample_mod
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

    private :: i, j, k, l
    integer(i4b) :: i, j, k, l
    
    logical(lgt) :: exist
    

contains

    subroutine sample_joint_amp(self, dat, compo, map_n, method, poltype)
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
        character(len=*), dimension(:), intent(in) :: poltype
        real(dp), allocatable, dimension(:,:)      :: A, val
        real(dp), allocatable, dimension(:)        :: b, c
        character(len=256)                         :: title
        integer(i4b)                               :: x, y, z, w, l, m, n
        integer(i4b)                               :: nfit1, nfit2, nfit3, nfit4
        integer(i4b)                               :: info


        real(dp), allocatable, dimension(:)        :: damps
        real(dp), allocatable, dimension(:,:)      :: synch
        real(dp)                                   :: q, t6, t7

        real(dp), allocatable, dimension(:,:,:)    :: covar, T_nu, T_nu_T, A_2, A_1
        real(dp), allocatable, dimension(:,:)      :: A_3


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

        do i = 1, size(self%joint_comp)
           if (self%joint_comp(i) == 'synch') then
              do l = 1, size(poltype)
                 y = y + dat%npix
              end do
           else if (self%joint_comp(i) == 'dust') then
              y = y + dat%npix
           else if (self%joint_comp(i) == 'template01') then
              y = y + self%temp_nfit(1)
           else if (self%joint_comp(i) == 'template02') then
              y = y + self%temp_nfit(2)
           else if (self%joint_comp(i) == 'template03') then
              y = y + self%temp_nfit(3)
           else if (self%joint_comp(i) == 'template04') then
              y = y + self%temp_nfit(4)
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
            if (self%joint_comp(m) == 'synch') then
               if (size(poltype) == 1) then
                  do j=1, z
                     do i=1, x
                        if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) then
                           c(i) = 0.d0
                           cycle
                        else
                           c(i) = c(i) + 1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*compute_spectrum(self,compo,1,self%dat_nu(j),i-1,map_n)
                        end if
                     end do
                  end do
                  w = w + x
               else if (size(poltype) == 2) then
                  do j=1, z
                     do i=1, x
                        if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) then
                           c(i)   = 0.d0
                           c(x+i) = 0.d0
                           cycle
                        else                           
                           c(i)   = c(i)   + 1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*compute_spectrum(self,compo,1,self%dat_nu(j),i-1,map_n)
                           c(x+i) = c(x+i) + 1.d0/(dat%rms_map(i-1,map_n+1,j)**2.d0)*dat%sig_map(i-1,map_n+1,j)*compute_spectrum(self,compo,1,self%dat_nu(j),i-1,map_n+1)
                        end if
                     end do
                  end do
                  w = w + 2*x
               end if
            else if (self%joint_comp(m) == 'dust') then
               do j=1, z
                  do i=1, x
                     if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) c(i) = 0.d0
                     c(i) = c(i) + 1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*compute_spectrum(self,compo,2,self%dat_nu(j),i-1,map_n)
                  end do
               end do
               w = w + x
            end if
        end do
        do m = 1, size(self%joint_comp)
           if (self%joint_comp(m) == 'template01') then
           ! Template 1
              if (size(poltype) == 1) then
                 l = 1
                 do j = 1, z
                    if (self%temp_corr(1,j)) then
                       do i = 1, x
                          if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                               dat%temps(i-1,map_n,1)
                       end do
                       l = l + 1
                    end if
                 end do
                 w = w + self%temp_nfit(1)
              else if (size(poltype) == 2) then
              ! If sampling Q and U jointly
                 l = 1
                 do j = 1, z
                    if (self%temp_corr(1,j)) then
                       do i = 1, x
                          if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                               dat%temps(i-1,map_n,1)
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n+1,j)**2.d0)*dat%sig_map(i-1,map_n+1,j)*&
                               dat%temps(i-1,map_n+1,1)
                       end do
                       l = l + 1
                    end if
                 end do
                 w = w + self%temp_nfit(1)
              end if

           else if (self%joint_comp(m) == 'template02') then
              ! Template 2
              if (size(poltype) == 1) then
                 l = 1
                 do j = 1, z
                    if (self%temp_corr(2,j)) then
                       do i = 1, x
                          if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                               dat%temps(i-1,map_n,2)
                       end do
                       l = l + 1
                    end if
                 end do
                 w = w + self%temp_nfit(2)
              else if (size(poltype) == 2) then
              ! If sampling Q and U jointly
                 l = 1
                 do j = 1, z
                    if (self%temp_corr(2,j)) then
                       do i = 1, x
                          if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                               dat%temps(i-1,map_n,2)
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n+1,j)**2.d0)*dat%sig_map(i-1,map_n+1,j)*&
                               dat%temps(i-1,map_n+1,2)
                       end do
                       l = l + 1
                    end if
                 end do
                 w = w + self%temp_nfit(2)
              end if

           else if (self%joint_comp(m) == 'template03') then
              ! Template 3
              if (size(poltype) == 1) then
                 l = 1
                 do j = 1, z
                    if (self%temp_corr(3,j)) then
                       do i = 1, x
                          if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                               dat%temps(i-1,map_n,3)
                       end do
                       l = l + 1
                    end if
                 end do
                 w = w + self%temp_nfit(3)
              else if (size(poltype) == 2) then
              ! If sampling Q and U jointly
                 l = 1
                 do j = 1, z
                    if (self%temp_corr(3,j)) then
                       do i = 1, x
                          if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                               dat%temps(i-1,map_n,3)
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n+1,j)**2.d0)*dat%sig_map(i-1,map_n+1,j)*&
                               dat%temps(i-1,map_n+1,3)
                       end do
                       l = l + 1
                    end if
                 end do
                 w = w + self%temp_nfit(3)
              end if

           else if (self%joint_comp(m) == 'template04') then
              ! Template 4
              if (size(poltype) == 1) then
                 l = 1
                 do j = 1, z
                    if (self%temp_corr(4,j)) then
                       do i = 1, x
                          if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                               dat%temps(i-1,map_n,4)
                       end do
                       l = l + 1
                    end if
                 end do
                 w = w + self%temp_nfit(4)
              else if (size(poltype) == 2) then
              ! If sampling Q and U jointly
                 l = 1
                 do j = 1, z
                    if (self%temp_corr(4,j)) then
                       do i = 1, x
                          if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                               dat%temps(i-1,map_n,4)
                          c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n+1,j)**2.d0)*dat%sig_map(i-1,map_n+1,j)*&
                               dat%temps(i-1,map_n+1,4)
                       end do
                       l = l + 1
                    end if
                 end do
                 w = w + self%temp_nfit(4)
              end if
           end if
        end do

        ! Computation
        if (trim(method) == 'cholesky') then
           write(*,*) 'Currently deprecated.'
           stop
           !mat_u(:,:)        = 0.d0
           !if (rank == master) write(*,fmt='(a)') 'Joint sampling using Cholesky Decomposition.'
           !call cholesky_decomp(A,mat_l)
           !mat_u  = transpose(mat_l)
           !call forward_sub(mat_l,d,c)
           !call backward_sub(mat_u,b,d)
        else if (trim(method) == 'cg') then
           if (rank == master) write(*,*) 'Joint sampling using CG.'
           ! Optimize
           !call compute_cg(A,b,c,y,nnz_a,self%cg_iter,self%cg_converge)
           !call compute_cg_precond(A,b,c,y,nnz_a,self%cg_iter,self%cg_converge)
           !call compute_cg_vec(b,c,y,self,dat,compo)
           ! Sample
           !call sample_cg(A,b,c,2*npix,nnz_a,self%cg_iter,self%cg_converge,self,dat,compo)
           call sample_cg_vec(b,c,y,self,dat,compo)
        else if (trim(method) == 'lu') then
           write(*,*) 'Currently deprecated.'
           stop
           !mat_u(:,:)        = 0.d0
           if (rank == master) write(*,*) 'Joint sampling using LU Decomp'
           !call LUDecomp(A,mat_l,mat_u,y)
           !call forward_sub(mat_l,d,c)
           !call backward_sub(mat_u,b,d)
        end if

        ! Output amplitudes to the appropriate variables
        if (size(poltype) == 1) then
           w = 0
           do m = 1, size(self%joint_comp)
              if (self%joint_comp(m) == 'synch') then
                 do i = 1, x
                    dat%fg_map(i-1,map_n,self%fg_ref_loc(1),1) = b(w+i)
                 end do
                 w = w + x
              else if (self%joint_comp(m) == 'dust') then
                 do i = 1, x
                    dat%fg_map(i-1,map_n,self%fg_ref_loc(2),3) = b(w+i)
                 end do
                 w = w + x
              end if
           end do
           do m = 1, size(self%joint_comp)
              if (self%joint_comp(m) == 'template01') then
                 l = 1
                 do while (l .lt. self%temp_nfit(1))
                    do j= 1, z
                       if (self%temp_corr(1,j)) then
                          dat%temp_amps(j,map_n,1) = b(w+l)
                          l = l + 1
                       else
                          dat%temp_amps(j,map_n,1) = 0.d0
                       end if
                    end do
                 end do
                 w = w + l -1
              else if (self%joint_comp(m) == 'template02') then
                 l = 1
                 do while (l .lt. self%temp_nfit(2))
                    do j= 1, z
                       if (self%temp_corr(2,j)) then
                          dat%temp_amps(j,map_n,2) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end do
              end if
           end do
        else if (size(poltype) == 2) then
           w = 0
           do m = 1, size(self%joint_comp)
              if (self%joint_comp(m) == 'synch') then
                 do i = 1, x
                    dat%fg_map(i-1,map_n,self%fg_ref_loc(1),1) = b(w+i)
                 end do
                 w = w + x
                 do i = 1, x
                    dat%fg_map(i-1,map_n+1,self%fg_ref_loc(1),1) = b(w+i)
                 end do
                 w = w + x
              else if (self%joint_comp(m) == 'dust') then
                 do i = 1, x
                    dat%fg_map(i-1,map_n,self%fg_ref_loc(2),3) = b(w+i)
                 end do
                 w = w + x
              end if
           end do
           do m = 1, size(self%joint_comp)
              if (self%joint_comp(m) == 'template01') then
                 l = 1
                 do while (l .lt. self%temp_nfit(1))
                    do j= 1, z
                       if (self%temp_corr(1,j)) then
                          dat%temp_amps(j,map_n,1)   = b(w+l)
                          dat%temp_amps(j,map_n+1,1) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end do
                 if (self%temp_nfit(1) == 1) then
                    do j= 1, z
                       if (self%temp_corr(1,j)) then
                          dat%temp_amps(j,map_n,1)   = b(w+l)
                          dat%temp_amps(j,map_n+1,1) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end if
                 w = w + l - 1
              else if (self%joint_comp(m) == 'template02') then
                 l = 1
                 do while (l .lt. self%temp_nfit(2))
                    do j= 1, z
                       if (self%temp_corr(2,j)) then
                          dat%temp_amps(j,map_n,2)   = b(w+l)
                          dat%temp_amps(j,map_n+1,2) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end do
                 if (self%temp_nfit(2) == 1) then
                    do j= 1, z
                       if (self%temp_corr(2,j)) then
                          dat%temp_amps(j,map_n,2)   = b(w+l)
                          dat%temp_amps(j,map_n+1,2) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end if
                 w = w + l - 1
              else if (self%joint_comp(m) == 'template03') then
                 l = 1
                 do while (l .lt. self%temp_nfit(3))
                    do j= 1, z
                       if (self%temp_corr(3,j)) then
                          dat%temp_amps(j,map_n,3)   = b(w+l)
                          dat%temp_amps(j,map_n+1,3) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end do
                 if (self%temp_nfit(3) == 1) then
                    do j= 1, z
                       if (self%temp_corr(3,j)) then
                          dat%temp_amps(j,map_n,3)   = b(w+l)
                          dat%temp_amps(j,map_n+1,3) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end if
                 w = w + l - 1
              else if (self%joint_comp(m) == 'template04') then
                 l = 1
                 do while (l .lt. self%temp_nfit(4))
                    do j= 1, z
                       if (self%temp_corr(4,j)) then
                          dat%temp_amps(j,map_n,4)   = b(w+l)
                          dat%temp_amps(j,map_n+1,4) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end do
                 if (self%temp_nfit(4) == 1) then
                    do j= 1, z
                       if (self%temp_corr(4,j)) then
                          dat%temp_amps(j,map_n,4)   = b(w+l)
                          dat%temp_amps(j,map_n+1,4) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end if
                 w = w + l - 1
              end if
           end do
        end if

        if (rank == master) then
           t3 = mpi_wtime()
           write(*,fmt='(a,f10.3,a)') 'Joint Sampler completed in ', t3-t1, 's.'
        end if

        write(*,*) ''

        ! Sure to deallocate all arrays here to free up memory
        deallocate(b)
        deallocate(c)
    end subroutine sample_joint_amp

    subroutine template_fit(self, dat, comp,map_n, temp_num)
      !------------------------------------------------------------------------
      ! Simple linear fit of a template to data map with a sampling term
      !------------------------------------------------------------------------
      implicit none
      
      type(params)                                           :: self
      type(component)                                        :: comp
      type(data)                                             :: dat
      !real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: dat, noise
      real(dp), dimension(0:npix-1,nbands)                   :: cov, nos, map
      integer(i4b),                               intent(in) :: map_n
      integer(i4b), optional,                     intent(in) :: temp_num
      real(dp)                                               :: temp, sum1, sum2, norm

      nos = dat%rms_map(:,map_n,:)
      cov = nos**2.d0

      if (trim(self%mode) == 'comp_sep') then

         write(*,*) "Template fitting NOT coded for general case yet!"
         stop
         !do j = 1, self%numband
         !   if (self%temp_corr(temp_num,j)) then
            
         !   end if
         !end do
      else if (trim(self%mode) == 'hi_fit') then
         do j = 1, self%numband
            sum1 = 0.d0
            sum2 = 0.d0
            norm = 0.d0
            temp = 0.d0
            if (self%temp_corr(1,j)) then
               do i = 0, npix-1
                  if (comp%HI(i,1) > self%thresh) cycle
                  temp = comp%HI(i,1)*planck(self%dat_nu(j)*1d9,comp%T_d(i,1))
                  sum1 = sum1 + (((dat%sig_map(i,map_n,j)-dat%offset(j))/dat%gain(j))*temp)/cov(i,j)
                  sum2 = sum2 + (temp)**2.d0/cov(i,j)
                  norm = norm + (temp)**2.d0/cov(i,j)
               end do
            end if
            comp%HI_amps(j) = sum1/sum2 + rand_normal(0.d0,1.d0)/sqrt(norm)
         end do
      end if

    end subroutine template_fit


   function sample_spec_amp(self, dat, comp, noise, ind, map_n)
        !------------------------------------------------------------------------
        ! Samples spectral amplitude (per pixel), following the spectrum of foreground 'type'. Returns a full map of amplitudes.
        !------------------------------------------------------------------------
        implicit none
  
        class(params)                                          :: self
        type(component)                                        :: comp
        real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: dat, noise
        integer(i4b),                               intent(in) :: ind
        integer(i4b),                               intent(in) :: map_n
        real(dp)                                               :: sum1, sum2, spec
        real(dp)                                               :: chi, chi_0, chi_00, p
        real(dp)                                               :: amp, num, t, sam
        real(dp), dimension(2)                                 :: pars
        real(dp), dimension(nbands)                            :: tmp
        real(dp), dimension(0:npix-1)                          :: norm
        real(dp), dimension(0:npix-1,nbands)                   :: cov, nos
        real(dp), dimension(0:npix-1)                          :: sample_spec_amp

        nos = noise(:,map_n,:)
        cov = nos**2.d0

        do i = 0, npix-1
            sum1    = 0.0d0
            sum2    = 0.0d0
            norm(i) = 0.d0
            do j = 1, nbands
                spec           = compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n)
                sum1           = sum1 + (dat(i,map_n,j)*spec)/cov(i,j)
                sum2           = sum2 + (spec)**2.d0/cov(i,j)
                norm(i)        = norm(i) + ((spec)**2.d0)/cov(i,j)
            end do
            norm(i)            = norm(i)
            amp                = sum1/sum2
            sample_spec_amp(i) = amp + rand_normal(0.d0,1.d0)/sqrt(norm(i))
        end do
    end function sample_spec_amp

    subroutine sample_index(self, dat, comp, duta, nside2, ind, map_n)
        implicit none
  
        class(params)                                          :: self
        type(component),                         intent(inout) :: comp
        type(data)                                             :: dat
        integer(i4b),                               intent(in) :: map_n, nside2, ind
        real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: duta
        integer(i4b)                                           :: nside1, npix2
        real(dp), dimension(0:npix-1,nmaps,nbands)             :: map2fit, cov 
        real(dp), dimension(0:npix-1,nmaps)                    :: indx
        real(dp), dimension(0:npix-1,nmaps)                    :: mask
        real(dp), dimension(0:npix-1)                          :: indx_sample
        real(dp), allocatable, dimension(:,:,:)                :: data_low, fg_map_low, rms_low
        real(dp), allocatable, dimension(:,:)                  :: indx_low, mask_low
        real(dp), allocatable, dimension(:)                    :: indx_sample_low
        real(dp), allocatable, dimension(:,:,:)                :: fg_map_high
        real(dp), dimension(nbands)                            :: signal, tmp
        real(dp), dimension(2)                                 :: x
        real(dp)                                               :: a, b, c, num, sam, t, p, sol
        real(dp)                                               :: time1, time2

        real(dp)                                               :: naccept   
        logical                                                :: exist

        real(dp), allocatable, dimension(:,:)                  :: map, map2
        character(len=128) ::title
 

        !------------------------------------------------------------------------
        ! Spectral index sampler, using the Metropolis approach.
        !------------------------------------------------------------------------

        map2fit = duta
        cov     = dat%rms_map*dat%rms_map
        mask    = dat%masks

        !------------------------------------------------------------------------
        ! Load priors for the appropriate spectrum
        !------------------------------------------------------------------------
        if (trim(self%fg_label(ind)) == 'synch') then 
            indx     = comp%beta_s
        else if (trim(self%fg_label(ind)) == 'dust') then 
            indx     = comp%beta_d
        end if
        
        if (rank == master) then
            if (mod(iter,self%iter_out) .EQ. 0) then
                write(*,fmt='(a,i4)') 'Sampling ' // trim(self%fg_label(ind)) // ' beta at nside', nside2
            end if
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
        allocate(indx_low(0:npix2-1,nmaps),rms_low(0:npix2-1,nmaps,nbands))
        allocate(indx_sample_low(0:npix2-1))
        allocate(mask_low(0:npix2-1,nmaps))
        allocate(fg_map_high(0:npix-1,nmaps,nbands))

        do i = 0, npix-1
           if (mask(i,self%pol_type(1)) == 0.d0) then
              mask(i,:) = missval
           end if
        end do

        fg_map_high(:,:,:) = dat%fg_map(:,:,:,1)

        if (nside1 /= nside2) then 
            if (ordering == 1) then
               call udgrade_ring(mask,nside1,mask_low,nside2,fmissval=missval,PESSIMISTIC=.false.)
               call udgrade_ring(indx,nside1,indx_low,nside2)
               do j = 1, nbands
                  call udgrade_ring(map2fit(:,:,j),nside1,data_low(:,:,j),nside2)
                  call udgrade_ring(fg_map_high(:,:,j),nside1,fg_map_low(:,:,j),nside2)
                  call udgrade_ring(cov(:,:,j),nside1,rms_low(:,:,j),nside2)
                  ! title = trim(self%outdir) // trim(self%dat_label(j)) // '_synch_amplitude_Q_n0004.fits'
                  ! map(:,1)   = fg_map_low(:,2,j)
                  ! do i = 0, npix2-1
                  !    if (mask_low(i,1) == 0.d0 .or. mask_low(i,1) == missval) then
                  !       map(i,1) = missval
                  !    end if
                  ! end do
                  ! call write_result_map(trim(title), nside2, ordering, header, map)
                  ! title = trim(self%outdir) // trim(self%dat_label(j)) // '_synch_amplitude_Q_n0064.fits'
                  ! map2(:,1)   = dat%fg_map(:,2,j,1)
                  ! do i = 0, npix-1
                  !    if (dat%masks(i,1) == 0.d0 .or. dat%masks(i,1) == missval) then
                  !       map2(i,1) = missval
                  !    end if
                  ! end do
                  ! call write_result_map(trim(title), nside1, ordering, header, map2)
               end do
            else
               call udgrade_nest(indx,nside1,indx_low,nside2)
               call udgrade_nest(mask,nside1,mask_low,nside2,fmissval=missval,PESSIMISTIC=.false.)
               do j = 1, nbands
                  call udgrade_nest(map2fit(:,:,j),nside1,data_low(:,:,j),nside2)
                  call udgrade_nest(dat%fg_map(:,:,j,1),nside1,fg_map_low(:,:,j),nside2)
                  call udgrade_nest(dat%rms_map(:,:,j),nside1,rms_low(:,:,j),nside2)
               end do
            end if
            rms_low = sqrt(rms_low / (npix/npix2))
            do i = 0, npix2-1
               if (mask_low(i,1) .lt. 0.50) then
                  mask_low(i,:) = 0.d0
               else
                  mask_low(i,:) = 1.d0
               end if
            end do
        else 
            do j = 1, nbands
                data_low(:,:,j)   = duta(:,:,j)
                fg_map_low(:,:,j) = dat%fg_map(:,:,j,1)
                rms_low(:,:,j)    = dat%rms_map(:,:,j)
            end do
            indx_low = indx
            mask_low = mask
        end if

        x(1) = 1.d0           
        !------------------------------------------------------------------------
        ! Metropolis algorithm:
        ! Sampling portion. Determine the log-likelihood, and accept based off of
        ! the improvement in the fit.
        !------------------------------------------------------------------------
        if (map_n == -1) then
           indx_sample_low = indx_low(:,self%pol_type(1))
           do i = 0, npix2-1
              a         = 0.d0
              sol       = indx_low(i,self%pol_type(1))
              sam       = sol
              if (mask_low(i,1) == 0.d0 .or. mask_low(i,1) == missval) then
                 cycle
              end if

              ! Chi-square from the most recent Gibbs chain update
              do j = 1, nbands
                 do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                    a = a + (((fg_map_low(i,k,self%fg_ref_loc(1)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,sol)) &
                         - data_low(i,k,j))**2.d0)/rms_low(i,k,j)**2.d0
                 end do
              end do
              c = a

              do l = 1, self%nsample

                 ! Sampling from the prior
                 if (self%fg_spec_like(ind,1)) then
                    t      = rand_normal(sol, self%fg_gauss(ind,1,2))
                 else if (.not. self%fg_spec_like(ind,1)) then
                    t      = rand_normal(self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2))
                 end if

                 ! Check if sample is outside the bounds
                 if (t .gt. self%fg_uni(ind,1,2) .or. t .lt. self%fg_uni(ind,1,1)) then
                    cycle
                 end if
                 if (sam .gt. self%fg_uni(ind,1,2) .or. sam .lt. self%fg_uni(ind,1,1)) then
                    cycle
                 end if

                 b         = 0.d0
                 
                 do j = 1, nbands
                    do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                       tmp(j) = fg_map_low(i,k,self%fg_ref_loc(1))*compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,t)
                       b      = b + ((tmp(j)-data_low(i,k,j))**2.d0)/rms_low(i,k,j)**2.d0
                    end do
                 end do

                 if (b < c .and. t .lt. self%fg_uni(ind,1,2) .and. t .gt. self%fg_uni(ind,1,1)) then
                    sam = t
                    c   = b
                 else
                    x(2) = exp(0.5d0*(c-b))
                    p = minval(x)
                    call RANDOM_NUMBER(num)
                    if (num < p .and. t .lt. self%fg_uni(ind,1,2) .and. t .gt. self%fg_uni(ind,1,1)) then
                       sam = t
                       c   = b
                    end if
                 end if
                 if (sam .gt. self%fg_uni(ind,1,2)) then
                    write(*,*) 'Sample outside of index bounds -- error!'
                    stop
                 end if
                 sol = sam
                 indx_sample_low(i) = sol
              end do              
           end do
        else
           do i = 0, npix2-1
              a         = 0.d0
              sol       = indx_low(i,map_n)
              sam       = sol
              
              if (mask_low(i,1) == 0.d0 .or. mask_low(i,1) == missval) then
                 cycle
              end if

              ! Chi-square from the most recent Gibbs chain update
              do j = 1, nbands
                 a = a + (((fg_map_low(i,map_n,self%fg_ref_loc(1)) * compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n,sol)) &
                      - data_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
              end do
              c = a
              
              do l = 1, self%nsample
                 
                 ! Sampling from the prior
                 if (self%fg_spec_like(ind,1)) then
                    t      = rand_normal(sol, self%fg_gauss(ind,1,2))
                  else if (.not. self%fg_spec_like(ind,1)) then
                    t      = rand_normal(self%fg_gauss(ind,1,1), self%fg_gauss(ind,1,2))
                 end if
                 b         = 0.d0
                 
                 if (t .gt. self%fg_uni(ind,1,2) .or. t .lt. self%fg_uni(ind,1,1)) cycle
                 
                 do j = 1, nbands
                    tmp(j) = fg_map_low(i,map_n,self%fg_ref_loc(1))*compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n,t)
                    b      = b + ((tmp(j)-data_low(i,map_n,j))**2.d0)/rms_low(i,map_n,j)**2.d0
                 end do
                 b = b
                 
                 if (b < c) then ! .and. t .lt. self%fg_uni(ind,1,2) .and. t .gt. self%fg_uni(ind,1,1)
                    sam = t
                    c   = b
                 else
                    x(2) = exp(0.5d0*(c-b))
                    p = minval(x)
                    call RANDOM_NUMBER(num)
                    if (num < p) then
                       sam = t
                       c   = b
                    end if
                 end if
              end do
              sol = sam
              indx_sample_low(i) = sol
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
              comp%beta_s(:,k) = indx_sample(:)
           else if (trim(self%fg_label(ind)) == 'dust') then 
              comp%beta_d(:,k) = indx_sample(:)
           end if
        end if

        time2 = mpi_wtime()

        write(*,*) ''
        if (rank == master) then
           time2 = mpi_wtime()
           write(*,fmt='(a,f10.3,a)') 'Spectral index sampler completed in ', time2-time1, 's.'
        end if

        write(*,*) ''

        deallocate(data_low)
        deallocate(fg_map_low)
        deallocate(indx_low)
        deallocate(rms_low)
        deallocate(indx_sample_low)
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
      real(dp), dimension(nbands)                            :: signal, tmp
      real(dp), dimension(2)                                 :: x
      real(dp)                                               :: a, b, c, num, sam, t, p, sol, naccept_t_d
      
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
      
      x(1) = 1.d0
      do i = 0, npix2-1
         if (dat%masks(i,1) == 0.d0  .or. dat%masks(i,1) == missval) then
            sample_T_low(i) = missval
            cycle
         else
            a   = 0.d0
            sol = T_low(i,map_n)
            
            ! Chi-square from the most recent Gibbs chain update
            do j = 1, nbands
               a = a + (((comp%HI_amps(j)*comp%HI(i,1)*planck(self%dat_nu(j)*1.d9,sol)) &
                    - (maps_low(i,map_n,j)-dat%offset(j))/dat%gain(j))**2.d0)/cov_low(i,map_n,j)
            end do
            c   = a
            
            do l = 1, self%nsample
               ! Begin sampling from the prior
               t = rand_normal(sol,self%HI_Td_std)
               !t = rand_normal(self%HI_Td_mean,self%HI_Td_std)
               b = 0.d0
               do j = 1, nbands
                  tmp(j) = comp%HI_amps(j)*comp%HI(i,1)*planck(self%dat_nu(j)*1.d9,t)
                  b      = b + ((tmp(j)-(maps_low(i,map_n,j)-dat%offset(j))/dat%gain(j))**2.d0)/cov_low(i,map_n,j)
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
            if (sam == 0.d0) then
               write(*,*) 'Error: T_d = 0.d0 accepted!'
               stop
            end if

            sol             = sam
            sample_T_low(i) = sol
         end if
      end do
      if (nside1 /= nside2) then
         if (ordering == 1) then
            call udgrade_ring(sample_T_low, nside2, te_sample, nside1)
            !call convert_nest2ring(nside2, sample_T_low)
         else
            call udgrade_nest(sample_T_low, nside2, te_sample, nside1)
         end if
      else
         te_sample =  sample_T_low
      end if
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

end module dang_sample_mod
