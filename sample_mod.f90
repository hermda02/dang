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

    private :: i, j, k, l
    integer(i4b) :: i, j, k, l
    
    logical(lgt) :: exist
    

contains

    subroutine sample_joint_amp(para, dat, compo, map_n, method, poltype)
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
        type(params)                               :: para
        type(data),                intent(inout)   :: dat
        type(component)                            :: compo
        integer(i4b),              intent(in)      :: map_n
        character(len=*),          intent(in)      :: method
        character(len=*), dimension(:), intent(in) :: poltype
        real(dp), allocatable, dimension(:,:)      :: A, val
        integer(i4b), allocatable, dimension(:,:)  :: col_ptr, row_ind
        real(dp), allocatable, dimension(:,:)      :: mat_l, mat_u, unc_a_s
        real(dp), allocatable, dimension(:)        :: b, c, d, rand, samp, unc_a_d
        character(len=256)                         :: title
        integer(i4b)                               :: x, y, z, nfit1, nfit2, w, l, m, n
        integer(i4b)                               :: vi, ci, ri, co, nnz, nnz_a
        integer(i4b)                               :: info

        real(dp)                                   :: q, t6, t7

        real(dp), allocatable, dimension(:,:,:)    :: covar, T_nu, T_nu_T, A_2, A_1
        real(dp), allocatable, dimension(:,:)      :: A_3

        if (rank == master) then
           write(*,fmt='(a)') 'Starting joint sampling for synch and dust_template.'
           !write(*,fmt='(a)') 'Pol_type = ', trim(poltype(:))
           t1 = mpi_wtime()
        end if

        ! Load which components to jointly fit for arary allocation
        ! vv These will not change based off of components
        x = dat%npix
        y = 0
        z = nbands

        do i = 1, size(compo%joint)
           if (compo%joint(i) == 'synch') then
              do l = 1, size(poltype)
                 ! write(*,*) trim(poltype(l))
                 y = y + dat%npix
              end do
           else if (compo%joint(i) == 'dust') then
              y = y + dat%npix
           else if (compo%joint(i) == 'template01') then
              ! Must have at most nbands -1 templates to avoid degeneracy
              ! Count how many bands are not being fit
              nfit1 = 0
              do j = 1, nbands
                 if (para%temp_corr(1,j)) then
                    nfit1 = nfit1 + 1
                 end if
              end do
              y = y + nfit1
           else if (compo%joint(i) == 'template02') then
              ! Must have at most nbands -1 templates to avoid degeneracy
              ! Count how many bands are not being fit
              nfit2 = 0
              do j = 1, nbands
                 if (para%temp_corr(2,j)) then
                    nfit2 = nfit2 + 1
                 end if
              end do
              y = y + nfit2
           end if
        end do
        
        !------------------------------------------------------------------------------|
        ! Since A is sparse (mostly 0s), we'll try to save time and memory by          |
        ! putting all of the values, with 'pointers' into a smaller array              |
        !                                                                              |
        ! We will use the CSC (Compressed Sparse Column) scheme to save time here.     |
        ! CSC is used instead of CSR due to the nature of the sparse matrix calc.      |
        !                                                                              |
        ! There are three arrays that will be used here, namely col_ptr, row_ind, and  |
        ! val. col_ptr is size y+1 (int), row_ind is nnz (num non-zero entries) (int), |
        ! and val is also size nnz (dp).                                               |
        !                                                                              |
        ! so for each band, we have nnz = 2*npix                                       |
        !------------------------------------------------------------------------------|

        nnz = 4*x

        allocate(b(y),c(y),d(y))

        ! Initialize arrays
        b(:)              = 0.d0
        c(:)              = 0.d0

        write(*,*) 'Compute RHS of matrix eqn.'
        ! Computing the LHS and RHS of the linear equation
        ! RHS
        w = 0 
        do m = 1, size(compo%joint)
            if (compo%joint(m) == 'synch') then
               if (size(poltype) == 1) then
                  do j=1, z
                     do i=1, x
                        ! if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
                        c(i) = c(i) + 1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n)
                     end do
                  end do
                  w = w + x
               else if (size(poltype) == 2) then
                  do j=1, z
                     do i=1, x
                        c(i) = c(i) + 1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n)
                     end do
                     do i=1, x
                        c(x+i) = c(x+i) + 1.d0/(dat%rms_map(i-1,map_n+1,j)**2.d0)*dat%sig_map(i-1,map_n+1,j)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n+1)
                     end do
                  end do
                  w = w + 2*x
               end if
            else if (compo%joint(m) == 'dust') then
               do j=1, z
                  do i=1, x
!                     if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
                     c(i) = c(i) + 1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*compute_spectrum(para,compo,2,para%dat_nu(j),i-1,map_n)
                  end do
               end do
               w = w + x
            end if
        end do
        do m = 1, size(compo%joint)
            if (compo%joint(m) == 'template01') then
               if (size(poltype) == 1) then
                  l = 1
                  do j = 1, z
                     if (para%temp_corr(1,j)) then
                        do i = 1, x
                           ! if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
                           c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                                    dat%temps(i-1,map_n,1)
                        end do
                        l = l + 1
                     end if
                  end do
                  w = w + l
               else if (size(poltype) == 2) then
                  l = 1
                  do j = 1, z
                     if (para%temp_corr(1,j)) then
                        do i = 1, x
                            c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                                    dat%temps(i-1,map_n,1)
                            c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n+1,j)**2.d0)*dat%sig_map(i-1,map_n+1,j)*&
                                    dat%temps(i-1,map_n+1,1)

                        end do
                        l = l + 1
                     end if
                  end do
                  w = w + l
               end if
            else if (compo%joint(m) == 'template02') then
                l = 1
                do j = 1, z
                    if (para%temp_corr(2,j)) then
                        do i = 1, x
                           ! if (mask(i,1) == 0.d0 .or. mask(i,1) == missval) cycle
                           c(w+l) = c(w+l)+1.d0/(dat%rms_map(i-1,map_n,j)**2.d0)*dat%sig_map(i-1,map_n,j)*&
                                    dat%temps(i-1,map_n,2)
                        end do
                        l = l + 1
                    end if
                end do
                w = w + l
            end if
        end do

        t3 = mpi_wtime()

        write(*,fmt='(a,E12.4,a)') 'Sparse matrix multiply: ', t3-t2, 's.'

        ! Computation
        if (trim(method) == 'cholesky') then
           !mat_u(:,:)        = 0.d0
           !if (rank == master) write(*,fmt='(a)') 'Joint sampling using Cholesky Decomposition.'
           !call cholesky_decomp(A,mat_l)
           !mat_u  = transpose(mat_l)
           !call forward_sub(mat_l,d,c)
           !call backward_sub(mat_u,b,d)
        else if (trim(method) == 'cg') then
           if (rank == master) write(*,*) 'Joint sampling using CG.'
           ! Optimize
           !call compute_cg(A,b,c,y,nnz_a,para%cg_iter,para%cg_converge)
           !call compute_cg_precond(A,b,c,y,nnz_a,para%cg_iter,para%cg_converge)
           !call compute_cg_vec(b,c,y,para,dat,compo)
           ! Sample
           !call sample_cg(A,b,c,2*npix,nnz_a,para%cg_iter,para%cg_converge,para,dat,compo)
           call sample_cg_vec(b,c,y,para,dat,compo)
        else if (trim(method) == 'lu') then
           !write(*,*) 'Currently deprecated -- replace with LAPACK'
           !stop
           !mat_u(:,:)        = 0.d0
           if (rank == master) write(*,*) 'Joint sampling using LU Decomp'
           !call LUDecomp(A,mat_l,mat_u,y)
           !call forward_sub(mat_l,d,c)
           !call backward_sub(mat_u,b,d)
        end if
        ! Draw a sample by cholesky decompsing A^-1, and multiplying 
        ! the subsequent lower triangular by a vector of random numbers
        !stop

        ! Output amplitudes to the appropriate variables
        if (size(poltype) == 1) then
           w = 0
           do m = 1, size(compo%joint)
              if (compo%joint(m) == 'synch') then
                 do i = 1, x
                    dat%fg_map(i-1,map_n,para%fg_ref_loc(1),1) = b(w+i)
                 end do
                 w = w + x
              else if (compo%joint(m) == 'dust') then
                 do i = 1, x
                    dat%fg_map(i-1,map_n,para%fg_ref_loc(2),3) = b(w+i)
                 end do
                 w = w + x
              end if
           end do
           do m = 1, size(compo%joint)
              if (compo%joint(m) == 'template01') then
                 l = 1
                 do while (l .lt. (nfit1))
                    do j= 1, z
                       if (para%temp_corr(1,j)) then
                          dat%temp_amps(j,map_n,1) = b(w+l)
                          l = l + 1
                       else
                          dat%temp_amps(j,map_n,1) = 0.d0
                       end if
                    end do
                 end do
              else if (compo%joint(m) == 'template02') then
                 l = 1
                 do while (l .lt. (nfit2))
                    do j= 1, z
                       if (para%temp_corr(2,j)) then
                          dat%temp_amps(j,map_n,2) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end do
              end if
           end do
        else if (size(poltype) == 2) then
           w = 0
           do m = 1, size(compo%joint)
              if (compo%joint(m) == 'synch') then
                 do i = 1, x
                    dat%fg_map(i-1,map_n,para%fg_ref_loc(1),1) = b(w+i)
                 end do
                 w = w + x
                 do i = 1, x
                    dat%fg_map(i-1,map_n+1,para%fg_ref_loc(1),1) = b(w+i)
                 end do
                 w = w + x
              else if (compo%joint(m) == 'dust') then
                 do i = 1, x
                    dat%fg_map(i-1,map_n,para%fg_ref_loc(2),3) = b(w+i)
                 end do
                 w = w + x
              end if
           end do
           do m = 1, size(compo%joint)
              if (compo%joint(m) == 'template01') then
                 l = 1
                 do while (l .lt. (nfit1))
                    do j= 1, z
                       if (para%temp_corr(1,j)) then
                          dat%temp_amps(j,map_n,1) = b(w+l)
                          dat%temp_amps(j,map_n+1,1) = b(w+l)
                          l = l + 1
                       else
                          dat%temp_amps(j,map_n,1) = 0.d0
                          dat%temp_amps(j,map_n+1,1) = 0.d0
                       end if
                    end do
                 end do
              else if (compo%joint(m) == 'template02') then
                 l = 1
                 do while (l .lt. (nfit2))
                    do j= 1, z
                       if (para%temp_corr(2,j)) then
                          dat%temp_amps(j,map_n,2) = b(w+l)
                          l = l + 1
                       end if
                    end do
                 end do
              end if
           end do
        end if

        if (rank == master) then
           t3 = mpi_wtime()
           write(*,fmt='(a,f10.3,a)') 'Joint Sampler completed in ', t3-t1, 's.'
        end if

        !if (para%output_unc .and. iter == niter) then
        !   allocate(unc_a_d(nfit1))
        !   allocate(unc_a_s(0:dat%npix-1,1))!

        !   write(*,*) 'Dust amplitude uncertainties: '
        
        !   call invert_matrix_dp(A,.true.)

        !   write(*,*) 'Done inverting.'

           
        !   if (size(poltype) == 1) then
        !      do j = 1, x
        !         unc_a_s(j-1,1) = sqrt(A(j,j))
        !      end do
        !      do j = 1, nfit1
        !         unc_a_d(j) = sqrt(A(x+j,x+j))
        !      end do
        !      inquire(file=trim(para%outdir) // 'dust_' // trim(tqu(k)) // '_uncertainties.dat',exist=exist)
        !      if (exist) then
        !         open(40,file = trim(para%outdir) // 'dust_' // trim(tqu(k)) // '_uncertainties.dat', status="old", &
        !              position="append", action="write")
        !      else
        !         open(40,file = trim(para%outdir) // 'dust_' // trim(tqu(k)) // '_uncertainties.dat', status="new", action="write")
        !      endif
        !      write(40,'(6(E17.8))') unc_a_d
        !      close(40)

        !      title = trim(para%outdir) // 'a_synch_uncertainty_'// trim(tqu(k)) // '.fits'
        !      call write_bintab(unc_a_s,dat%npix,1, header, nlheader, trim(title))

        !   else if (size(poltype) == 2) then
        !      do j = 1, x
        !         unc_a_s(j-1,1) = sqrt(A(j,j))
        !      end do
        !      do j = 1, nfit1
        !         unc_a_d(j) = sqrt(A(2*x+j,2*x+j))
        !      end do
              
        !      inquire(file=trim(para%outdir) // 'dust_QU_uncertainties.dat',exist=exist)
        !      if (exist) then
        !         open(40,file = trim(para%outdir) // 'dust_QU_uncertainties.dat', status="old", &
        !              position="append", action="write")
        !      else
        !         open(40,file = trim(para%outdir) // 'dust_QU_uncertainties.dat', status="new", action="write")
        !      endif
        !      write(40,'(6(E17.8))') unc_a_d
        !      close(40)

        !   end if
                  
        !end if

        write(*,*) 'Exit joint_sampler'

        ! Sure to deallocate all arrays here to free up memory
        !deallocate(A)
        deallocate(b)
        deallocate(c)
        deallocate(d)
        !deallocate(mat_l)
        !deallocate(mat_u)
        !deallocate(rand)

    end subroutine sample_joint_amp


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

    subroutine sample_index(self, comp, dat, duta, nside2, ind, map_n)
        implicit none
  
        class(params)                                          :: self
        type(component),                         intent(inout) :: comp
        type(data)                                             :: dat
        integer(i4b),                               intent(in) :: map_n, nside2, ind
        real(dp), dimension(0:npix-1,nmaps,nbands), intent(in) :: duta
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

        map2fit = duta
        cov     = dat%rms_map*dat%rms_map

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
                    call udgrade_ring(dat%fg_map(:,:,j,1),nside1,fg_map_low(:,:,j),nside2)
                    call convert_nest2ring(nside1,dat%fg_map(:,:,j,1))
                    call udgrade_ring(cov(:,:,j),nside1,rms_low(:,:,j),nside2)
                    call convert_nest2ring(nside2,rms_low(:,:,j))
                else
                    call udgrade_nest(map2fit(:,:,j),nside1,data_low(:,:,j),nside2)
                    call udgrade_nest(dat%fg_map(:,:,j,1),nside1,fg_map_low(:,:,j),nside2)
                    call udgrade_nest(dat%rms_map(:,:,j),nside1,rms_low(:,:,j),nside2)
               end if
            end do
            rms_low = sqrt(rms_low / (npix/npix2))
        else 
            do j = 1, nbands
                data_low(:,:,j)   = duta(:,:,j)
                fg_map_low(:,:,j) = dat%fg_map(:,:,j,1)
                rms_low(:,:,j)    = dat%rms_map(:,:,j)
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
                 b         = 0.d0
                 
                 do j = 1, nbands
                    do k = self%pol_type(1), self%pol_type(size(self%pol_type))
                       tmp(j) = fg_map_low(i,k,self%fg_ref_loc(1))*compute_spectrum(self,comp,ind,self%dat_nu(j),i,k,t)
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
                 
                 do j = 1, nbands
                    tmp(j) = fg_map_low(i,map_n,self%fg_ref_loc(1))*compute_spectrum(self,comp,ind,self%dat_nu(j),i,map_n,t)
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
                 comp%beta_s(:,k) = indx_sample
              else if (trim(self%fg_label(ind)) == 'dust') then 
                 comp%beta_d(:,k) = indx_sample
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
