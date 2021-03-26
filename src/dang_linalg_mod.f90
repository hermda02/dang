module dang_linalg_mod
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
    use dang_util_mod
    use dang_param_mod
    use dang_component_mod
    use dang_data_mod
    implicit none

contains

    subroutine invert_matrix_dp(matrix, cholesky, status, ln_det)
      implicit none
      
      real(dp), dimension(1:,1:), intent(inout)         :: matrix
      logical(lgt),               intent(in),  optional :: cholesky
      integer(i4b),               intent(out), optional :: status
      real(dp),                   intent(out), optional :: ln_det
      
      integer(i4b)     :: i, j, n, lda, info, lwork
      logical(lgt)     :: use_cholesky
      character(len=1) :: uplo
      integer(i4b), allocatable, dimension(:)   :: ipiv
      real(dp),     allocatable, dimension(:)   :: work
      
      if(present(status)) status = 0
      use_cholesky = .false.; if (present(cholesky)) use_cholesky = cholesky
      n     = size(matrix(1,:))
      lda   = n
      lwork = n
      info  = 0
      uplo  = 'l'
      allocate(ipiv(n))
      allocate(work(n))
      
      if (use_cholesky) then
         call DPOTRF(uplo, n, matrix, lda, info)
         if (present(ln_det)) then
            ln_det = 0.d0
            do i = 1, n
               if (matrix(i,i) > 0.d0) then
                  ln_det = ln_det + 2.d0*log(matrix(i,i))
               else
                  ln_det = -1.d30
                  exit
               end if
            end do
         end if
      else
         call DGETRF(n, n, matrix, lda, ipiv, info)
      end if
      if (info /= 0) then
         if(present(status)) then
            status = info
            return
         end if
         write(*,*) 'DGETRF: Factorization failed. Info = ', info
         stop
      else
         
         if (use_cholesky) then
            call DPOTRI(uplo, n, matrix, lda, info)
         else
            call DGETRI(n, matrix, lda, ipiv, work, lwork, info)
         end if
         
         if (info /= 0) then
            if(present(status)) then
               status = info
               return
            end if
            write(*,*) 'DGETRI: Inversion failed. Info = ', info
            stop
         end if
         
      end if

      if (use_cholesky) then
         do i = 1, n
            do j = i+1, n
               matrix(i,j) = matrix(j,i)
            end do
         end do
      end if
      
      deallocate(work)
      deallocate(ipiv)

    end subroutine invert_matrix_dp

    subroutine invert_tri(tri,uplo,inv)
      implicit none
      real(dp), dimension(:,:), intent(in)  :: tri
      character(len=1),         intent(in)  :: uplo
      real(dp), dimension(:,:), intent(out) :: inv
      integer(i4b)                          :: n, info

      n   = size(tri(1,:))

      inv = tri

      call dtrtri(uplo,'N',n,inv,n,info)

    end subroutine invert_tri

    subroutine cholesky_decomp(mat,low)
        implicit none
        real(dp), dimension(:,:), intent(in)  :: mat
        real(dp), dimension(:,:), intent(out) :: low
        integer(i4b)                          :: n, info

        n   = size(mat(1,:))

        low = mat

        call dpotrf('L',n,low,n,info)

    end subroutine cholesky_decomp

    subroutine LUDecomp(A,L,U,n)
        
        ! Using Doolittle's method for LU decomposition of matrix A
        ! L and U can then be used to solve the matrix equation Ax = b by solving
        ! Lv = b (using forward substitution), Ux = v (backwards substitution).
        
        implicit none        
        real(dp), dimension(:,:), intent(in)  :: A
        real(dp), dimension(:,:), intent(out) :: L
        real(dp), dimension(:,:), intent(out) :: U
        integer(i4b), intent(in)              :: n
        integer(i4b)                          :: i, j, m
        real(dp)                              :: sum

        L = 0.d0
        U = 0.d0

        do i = 1, n
            do m = i, n
                sum = 0
                do j = 1, i
                   sum = sum + (L(i,j) * U(j,m))
                   U(i,m) = A(i,m) - sum
                end do
            end do
            do m = i, n
                if (i == m) then
                    L(i,i) = 1
                else
                    sum = 0
                    do j = 1, i
                        sum = sum + (L(m,j)*U(j,i))
                    end do
                    L(m,i) = (A(m,i) - sum)/U(i,i)
                end if
            end do
        end do
    end subroutine LUDecomp

    subroutine forward_sub(L,xf,bf)
        ! Forward substitution to solve the matrix equation Lx=b
        implicit none
        real(dp), dimension(:,:), intent(in)  :: L
        real(dp), dimension(:),   intent(in)  :: bf
        real(dp), dimension(:),   intent(out) :: xf
        integer(i4b)                          :: i, j, n

        n = size(bf)

        do i = 1, n
            xf(i) = bf(i)
            !$OMP PARALLEL PRIVATE(j)
            !$OMP DO SCHEDULE(static)
            do j = 1, i-1
                xf(i) = xf(i)-l(i,j)*xf(j) 
            end do
            !$OMP END DO
            !$OMP END PARALLEL
            xf(i) = xf(i)/L(i,i)
        end do

    end subroutine forward_sub

    subroutine backward_sub(U,xb,bb)
        ! Backward substitution to solve the matrix equation Ux=b
        implicit none
        real(dp), dimension(:,:), intent(in)  :: U
        real(dp), dimension(:),   intent(in)  :: bb
        real(dp), dimension(:),   intent(out) :: xb
        integer(i4b)                          :: i, j, n

        n = size(bb)
        xb(n) = bb(n)/U(n,n)

        do i = n-1, 1, -1
            xb(i) = bb(i)
            !$OMP PARALLEL PRIVATE(i)
            !$OMP DO SCHEDULE(static)
            do j = n, i+1, -1
                xb(i) = xb(i) - U(i,j)*xb(j)
            end do
            !$OMP END DO
            !$OMP END PARALLEL
            xb(i) = xb(i)/U(i,i)
        end do

    end subroutine backward_sub

    subroutine sample_cg_vec(x,b,param,dat,compos, map_n)
        
      ! Implementation of the canned algorithm (B2) outlined in Jonathan Richard Shewuck (1994)
      ! "An introduction to the Conjugate Gradient Method Without the Agonizing Pain"

      implicit none
      type(params)                          :: param
      type(data)                            :: dat
      type(component)                       :: compos
      real(dp), dimension(:), intent(in)    :: b
      real(dp), dimension(:), intent(out)   :: x
      integer(i4b),           intent(in)    :: map_n
      real(dp), allocatable, dimension(:)   :: r, q, d, s, b2, eta
      real(dp)                              :: converge
      real(dp)                              :: alpha, beta, delta_0
      real(dp)                              :: delta_old, delta_new, t3, t4, t5, t6
      integer(i4b)                          :: i_max, i, j, n, m

      m = size(b)

      t5 = mpi_wtime()
      
      allocate(r(m),q(m),d(m))
      allocate(b2(m))

      ! Determine length of univariate with
      ! is also the size of the square matrix N
      if (param%joint_pol) then
         n = dat%npix*2
      else
         n = dat%npix
      end if

      allocate(eta(n))

      ! Draw a univariate
      do i = 1, n
         eta(i) = rand_normal(0.d0,1.d0)
      end do

      if (trim(param%ml_mode) == 'sample') then
         b2 = b + compute_sample_vec(param, dat, compos, eta, param%numband, map_n, b)
      else if (trim(param%ml_mode) == 'optimize') then
         b2 = b
      end if

      ! Initial condition parameters
      !-----------------------------
      x(:)     = 0.0d0
      i_max    = param%cg_iter
      converge = param%cg_converge
      !-----------------------------

      r  = b2 - return_Ax(param, dat, compos, x, param%numband, map_n)
      d  = r
      delta_new = sum(r*r)
      delta_0   = delta_new
      i = 0     
      do while( (i .lt. i_max) .and. (delta_new .gt. converge))
         t3    = mpi_wtime()
         q     = return_Ax(param, dat, compos, d, param%numband, map_n)
         alpha = delta_new/(sum(d*q))
         x     = x + alpha*d

         if (mod(i,50) == 0) then
            r = b2 - return_Ax(param, dat, compos, x, param%numband, map_n)
         else
            r = r - alpha*q
         end if

         delta_old = delta_new
         delta_new = sum(r*r)
         beta      = delta_new/delta_old
         d         = r + beta*d
         t4        = mpi_wtime()
         i         = i + 1

         if (delta_new .gt. converge) then
            write(*,fmt='(a,i4,a,e12.5,a,e12.5,a)') 'CG Iter: ', i, ' | delta: ', delta_new, ' | time: ', t4-t3, 's.'
         else
            write(*,fmt='(a,i4,a,e12.5,a,e12.5,a)') 'Final CG Iter: ', i, ' | delta: ', delta_new, ' | time: ', t4-t3, 's.'
         end if
      end do

      t6 = mpi_wtime()
      
      write(*,fmt='(a,e12.5,a)') 'CG Total time: ', t6-t5, 's.'
      
      deallocate(r)
      deallocate(q)
      deallocate(d)

    end subroutine sample_cg_vec

    subroutine compute_cg_precond(A,x,b,n,nnz_a,iters,converge)
        
      ! Implementation of the canned algorithm (B2) outlined in Jonathan Richard Shewuck (1994)
      ! "An introduction to the Conjugate Gradient Method Without the Agonizing Pain"

      implicit none
      real(dp), dimension(:,:), intent(in)  :: A
      real(dp), dimension(:), intent(in)    :: b
      real(dp), dimension(:), intent(out)   :: x
      integer(i4b), intent(in)              :: n,iters
      real(dp), intent(in)                  :: converge
      integer(i4b), optional, intent(in)    :: nnz_a
      real(dp), allocatable, dimension(:)   :: r, q, d, s, M 
      real(dp)                              :: epsil, alpha, beta, delta_0
      real(dp)                              :: delta_old, delta_new, t3, t4, t5, t6
      integer(i4b)                          :: i_max, i, j
        
      real(dp),allocatable,dimension(:)     :: v
      integer(i4b),allocatable,dimension(:) :: rp, ci
      
      t5 = mpi_wtime()
      allocate(r(n),q(n),d(n),M(n),s(n))
      
      if(present(nnz_a)) then
         allocate(v(nnz_a),rp(n+1),ci(nnz_a))
         call A_to_CSR(A,rp,ci,v)
      end if
      
      x(:)     = 0.0d0
      M(:)     = 1.d0
      i_max = iters
      
      i = 0

      write(*,*) 'Define preconditioner'
      do i = 1, n
         M(i) = 1.d0/A(i,i)
      end do

      i = 0

      r = b

      !$OMP PARALLEL PRIVATE(i)
      !$OMP DO SCHEDULE(static)
      do i = 1, n
         d(i) = M(i)*r(i)
      end do
      !$OMP END DO
      !$OMP END PARALLEL
      delta_new = sum(r*d)
      delta_0   = delta_new
      write(*,*) delta_0
      t3 = mpi_wtime()
      do while( (i .lt. i_max) .and. (delta_new .gt. converge))!(epsil**2)*delta_0))
         t3 = mpi_wtime()
         if (present(nnz_a)) then
            q = Ax_CSR(rp,ci,v,d)
         else
            q = matmul(A,d)
         end if
         alpha = delta_new/(sum(d*q))
         x = x + alpha*d
         if (mod(i,50) == 0) then
            if (present(nnz_a)) then
               r = b - Ax_CSR(rp,ci,v,x)
            else
               r = b - matmul(A,x)
            end if
         else
            r = r - alpha*q
         end if
         !$OMP PARALLEL PRIVATE(i)
         !$OMP DO SCHEDULE(static)
         do j = 1, n
            s(j) = M(j)*r(j)
         end do
         !$OMP END DO
         !$OMP END PARALLEL
         delta_old = delta_new
         delta_new = sum(r*s)
         beta = delta_new/delta_old
         d  = s + beta*d
         t4 = mpi_wtime()
         i = i + 1
         if (delta_new .gt. converge) then
            write(*,fmt='(a,i4,a,e12.5,a,e12.5,a)') 'CG Iter: ', i, ' | delta: ', delta_new, ' | time: ', t4-t3, 's.'
         else
            write(*,fmt='(a,i4,a,e12.5,a,e12.5,a)') 'Final CG Iter: ', i, ' | delta: ', delta_new, ' | time: ', t4-t3, 's.'
         end if
      end do
      t6 = mpi_wtime()

      write(*,fmt='(a,e12.5,a)') 'CG Total time: ', t6-t5, 's.'

      deallocate(r)
      deallocate(q)
      deallocate(d)

    end subroutine compute_cg_precond

    function compute_sample_vec(self, dat, compo, vec, nbands, map_n, rhs) result(res)
      !-----------------------------------------------------------------
      ! Function to compute the sampling vector for the CG method
      !
      ! Sampling vector has the form T^t N^-1/2 eta, where eta is a 
      ! univariate vector of length size(b).
      !-----------------------------------------------------------------
      implicit none
      type(params)                           :: self
      type(data)                             :: dat
      type(component)                        :: compo

      real(dp), dimension(:), intent(in)     :: vec
      real(dp), dimension(:), intent(in)     :: rhs
      real(dp), allocatable, dimension(:)    :: temp1, temp2, temp3, res
      real(dp)                               :: x
      integer(i4b), intent(in)               :: nbands, map_n
      integer(i4b)                           :: i, j, k, len, l, m, n, r, t, s
      integer(i4b)                           :: z, w, o

      x = dat%npix
      n = size(vec)
      m = size(rhs)
      len = 0

      do i = 1, size(self%joint_comp)
         if (self%joint_comp(i) == 'synch') then
            if (self%joint_pol) then
               len = len + 2*x
            else
               len = len + x
            end if
         end if
      end do

      ! How this works: result = sum_nu(T_nu^t N_nu^{-1/2})vec
      !--------------------------------------------------------------
      ! Suppose T_nu is of size n, m, and input vector is of size m |
      ! Then N^{-1/2} is of size n, n                               |
      ! And T_nu^t is of size m, n                                  |
      !                                                             |
      ! The returned vector, like the input vector, is of size m    |
      !                                                             |
      ! temp1 is of size n                                          |
      ! temp2 is also of size n                                     |
      ! temp3 is of size m <- the final vector                      |
      !--------------------------------------------------------------

      allocate(temp1(n))
      allocate(temp2(m))
      allocate(res(m))

      res = 0.d0

      do j = 1, nbands
         temp1 = 0.d0
         temp2 = 0.d0
         temp3 = 0.d0

         ! Step one is to solve temp1 = (N^{-1/2})(vec)
         !$OMP PARALLEL PRIVATE(i)
         !$OMP DO SCHEDULE(static)
         do i = 1, x
            if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
            temp1(i) = vec(i)/(dat%rms_map(i-1,map_n,j))
            if (self%joint_pol) temp1(x+i) = vec(x+i)/(dat%rms_map(i-1,map_n+1,j))
         end do
         !$OMP END DO
         !$OMP END PARALLEL

         ! Step three is to solve temp3 = (T_nu^t)(temp2) - this part is complex/confusing...
         do m = 1, size(self%joint_comp)
            do n = 1, self%ncomp
               if (trim(self%joint_comp(m)) == trim(self%fg_label(n))) then
                  !$OMP PARALLEL PRIVATE(i)
                  !$OMP DO SCHEDULE(static)
                  do i = 1, x
                     if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                     temp2(i) = temp1(i)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n)
                     if (self%joint_pol) temp2(x+i) = temp1(x+i)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n+1)
                  end do
                  !$OMP END DO
                  !$OMP END PARALLEL
               end if
            end do
            do n = 1, self%ntemp
               if (trim(self%joint_comp(m)) == trim(self%temp_label(n))) then
                  if (self%temp_corr(n,j)) then
                     !!$OMP PARALLEL PRIVATE(i)
                     !!$OMP DO SCHEDULE(static)
                     do i = 1, x
                        if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                        temp2(len+l) = temp2(len+l) + temp1(i)*dat%temps(i-1,map_n,n)
                        if (self%joint_pol) temp2(len+l) = temp2(len+l) + temp1(x+i)*dat%temps(i-1,map_n+1,n)
                     end do
                     !!$OMP END DO
                     !!$OMP END PARALLEL
                     l = l + 1
                  end if
               end if
            end do
         end do
         res = res + temp2
      end do

      ! do m = 1, size(self%joint_comp)
      !    do n = 1, self%ncomp
      !       if (trim(self%joint_comp(m)) == trim(self%fg_label(n))) then
      !          if (.not. self%joint_pol) then
      !             do j=1, z
      !                do i=1, x
      !                   if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) then
      !                      res(i) = 0.d0
      !                      cycle
      !                   else
      !                      res(i) = res(i) + 1.d0/(dat%rms_map(i-1,map_n,j))*vec(i)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n)
      !                   end if
      !                end do
      !             end do
      !             w = w + x
      !          else if (self%joint_pol) then
      !             ! If sampling Q and U jointly
      !              do j=1, z
      !                do i=1, x
      !                   if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) then
      !                      res(i)   = 0.d0
      !                      res(x+i) = 0.d0
      !                      cycle
      !                   else                           
      !                      res(i)   = res(i)   + 1.d0/(dat%rms_map(i-1,map_n,j))*vec(i)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n)
      !                      res(x+i) = res(x+i) + 1.d0/(dat%rms_map(i-1,map_n+1,j))*vec(i)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n+1)
      !                   end if
      !                end do
      !             end do
      !             w = w + 2*x
      !          end if
      !       end if
      !    end do
      !    do n = 1, self%ntemp
      !       if (trim(self%joint_comp(m)) == trim(self%temp_label(n))) then
      !          if (.not. self%joint_pol) then
      !             l = 1
      !             do j = 1, z
      !                if (self%temp_corr(n,j)) then
      !                   do i = 1, x
      !                      if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
      !                      res(w+l) = res(w+l)+1.d0/(dat%rms_map(i-1,map_n,j))*vec(i)*&
      !                           dat%temps(i-1,map_n,n)
      !                   end do
      !                   l = l + 1
      !                end if
      !             end do
      !             w = w + self%temp_nfit(n)
      !          else if (self%joint_pol) then
      !             ! If sampling Q and U jointly
      !             l = 1
      !             do j = 1, z
      !                if (self%temp_corr(n,j)) then
      !                   do i = 1, x
      !                      if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
      !                      res(w+l) = res(w+l)+1.d0/(dat%rms_map(i-1,map_n,j))*vec(i)*&
      !                           dat%temps(i-1,map_n,n)
      !                      res(w+l) = res(w+l)+1.d0/(dat%rms_map(i-1,map_n+1,j))*vec(i)*&
      !                           dat%temps(i-1,map_n+1,n)
      !                   end do
      !                   l = l + 1
      !                end if
      !             end do
      !             w = w + self%temp_nfit(n)
      !          end if
      !       end if
      !    end do
      ! end do

    end function compute_sample_vec

    subroutine multiply_with_A(self, dat, compo, vec, nbands, map_n)
      implicit none
      type(params)                           :: self
      type(data)                             :: dat
      type(component)                        :: compo

      real(dp), dimension(:), intent(inout)  :: vec
      real(dp), allocatable, dimension(:)    :: temp1, temp2, temp3, res
      real(dp)                               :: x
      integer(i4b), intent(in)               :: nbands, map_n
      integer(i4b)                           :: n, i, j, k, len, l, m, r, t, s

      x = dat%npix
      n = size(vec)
      len = 0

      do i = 1, size(self%joint_comp)
         if (self%joint_comp(i) == 'synch') then
            if (self%joint_pol) then
               len = len + 2*x
            else
               len = len + x
            end if
         end if
      end do

      ! How this works: result = sum_nu(T_nu^t N_nu^-1 T_nu)vec
      !--------------------------------------------------------------
      ! Suppose T_nu is of size n, m, and input vector is of size m |
      ! Then N^-1 is of size n, n                                   |
      ! And T_nu^t is of size m, n                                  |
      !                                                             |
      ! The returned vector, like the input vector, is of size m    |
      !                                                             |
      ! temp1 is of size n                                          |
      ! temp2 is also of size n                                     |
      ! temp3 is of size m <- the final vector                      |
      !--------------------------------------------------------------

      allocate(temp1(len))
      allocate(temp2(len))
      allocate(temp3(n))
      allocate(res(n))

      res     = 0.d0

      !---------------------------------
      ! No multi-template handling yet
      ! l = 1
      ! if (self%ntemp == 2) then 
      !    r = 1+self%temp_nfit(1)
      ! else if (self%ntemp == 3) then
      !    r = 1+self%temp_nfit(1)
      !    t = 1+self%temp_nfit(1)+self%temp_nfit(2)
      ! else if (self%ntemp == 4) then 
      !    r = 1+self%temp_nfit(1)
      !    t = 1+self%temp_nfit(1)+self%temp_nfit(2)
      !    s = 1+self%temp_nfit(1)+self%temp_nfit(2)+self%temp_nfit(3)
      ! end if
      !---------------------------------

      do j = 1, nbands
         temp1 = 0.d0
         temp2 = 0.d0
         temp3 = 0.d0
         ! A is composed of sum_nu(T_nu^T N_nu^-1 T_nu)

         ! Step one is to solve temp1 = (T_nu)(vec)
         do m = 1, size(self%joint_comp)
            do n = 1, self%ncomp
               if (trim(self%joint_comp(m)) == trim(self%fg_label(n))) then
                  ! First the multiply npix values of x by the npix diagonals of T_nu
                  !$OMP PARALLEL PRIVATE(i)
                  !$OMP DO SCHEDULE(static)
                  do i = 1, x
                     if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                     temp1(i) = compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n)*vec(i)
                     if (self%joint_pol) temp1(x+i) = compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n+1)*vec(x+i)
                   end do
                   !$OMP END DO
                   !$OMP END PARALLEL
                end if
            end do
            do n = 1, self%ntemp
               if (trim(self%joint_comp(m)) == trim(self%temp_label(n))) then
                  ! Then add vec(i)*template(i) to each temp1(i)
                  if (self%temp_corr(n,j)) then
                     !$OMP PARALLEL PRIVATE(i)
                     !$OMP DO SCHEDULE(static)
                     do i = 1, x
                        if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                        temp1(i) = temp1(i) + dat%temps(i-1,map_n,n)*vec(len+l)
                        if (self%joint_pol) temp1(x+i) = temp1(x+i) + dat%temps(i-1,map_n+1,n)*vec(len+l)
                     end do
                     !$OMP END DO
                     !$OMP END PARALLEL
                  end if
               end if
            end do
         end do
         ! Step two is to solve temp2 = (N^-1)(temp1)
         !$OMP PARALLEL PRIVATE(i)
         !$OMP DO SCHEDULE(static)
         do i = 1, x
            if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
            temp2(i) = temp1(i)/(dat%rms_map(i-1,map_n,j)**2.d0)
            if (self%joint_pol) temp2(x+i) = temp1(x+i)/(dat%rms_map(i-1,map_n+1,j)**2.d0)
         end do
         !$OMP END DO
         !$OMP END PARALLEL
         ! Step three is to solve temp3 = (T_nu^t)(temp2) - this part is complex/confusing...
         do m = 1, size(self%joint_comp)
            do n = 1, self%ncomp
               if (trim(self%joint_comp(m)) == trim(self%fg_label(n))) then
                  !$OMP PARALLEL PRIVATE(i)
                  !$OMP DO SCHEDULE(static)
                  do i = 1, x
                     if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                     temp3(i) = temp2(i)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n)
                     if (self%joint_pol) temp3(x+i) = temp2(x+i)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n+1)
                  end do
                  !$OMP END DO
                  !$OMP END PARALLEL
               end if
            end do
            do n = 1, self%ntemp
               if (self%temp_corr(n,j)) then
                  !$OMP PARALLEL PRIVATE(i)
                  !$OMP DO SCHEDULE(static)
                  do i = 1, x
                     if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                     temp3(len+l) = temp2(i)*dat%temps(i-1,map_n,n)
                     if (self%joint_pol) temp3(len+l) = temp3(len+l) + temp2(x+i)*dat%temps(i-1,map_n+1,n)
                  end do
                  !$OMP END DO
                  !$OMP END PARALLEL
               end if
            end do
         end do
         res = res + temp3
      end do

      vec = res

    end subroutine multiply_with_A

    function return_Ax(self, dat, compo, vec, nbands, map_n) result(res)
      implicit none
      type(params)                           :: self
      type(data)                             :: dat
      type(component)                        :: compo

      real(dp), dimension(:), intent(in)     :: vec
      real(dp), allocatable, dimension(:)    :: temp1, temp2, temp3, res
      real(dp)                               :: x
      integer(i4b), intent(in)               :: nbands, map_n
      integer(i4b)                           :: n, i, j, k, len, l, m, r, t, s

      x = dat%npix
      n = size(vec)
      len = 0

      do i = 1, size(self%joint_comp)
         if (self%joint_comp(i) == 'synch') then
            if (self%joint_pol) then
               len = len + 2*x
            else
               len = len + x
            end if
         end if
      end do

      ! How this works: result = sum_nu(T_nu^t N_nu^-1 T_nu)vec
      !--------------------------------------------------------------
      ! Suppose T_nu is of size n, m, and input vector is of size m |
      ! Then N^-1 is of size n, n                                   |
      ! And T_nu^t is of size m, n                                  |
      !                                                             |
      ! The returned vector, like the input vector, is of size m    |
      !                                                             |
      ! temp1 is of size n                                          |
      ! temp2 is also of size n                                     |
      ! temp3 is of size m <- the final vector                      |
      !--------------------------------------------------------------

      allocate(temp1(len))
      allocate(temp2(len))
      allocate(temp3(n))
      allocate(res(n))

      res     = 0.d0

      !---------------------------------
      ! No multi-template handling yet
      l = 1
      ! if (self%ntemp == 2) then 
      !    r = 1+self%temp_nfit(1)
      ! else if (self%ntemp == 3) then
      !    r = 1+self%temp_nfit(1)
      !    t = 1+self%temp_nfit(1)+self%temp_nfit(2)
      ! else if (self%ntemp == 4) then 
      !    r = 1+self%temp_nfit(1)
      !    t = 1+self%temp_nfit(1)+self%temp_nfit(2)
      !    s = 1+self%temp_nfit(1)+self%temp_nfit(2)+self%temp_nfit(3)
      ! end if
      !---------------------------------

      do j = 1, nbands
         temp1 = 0.d0
         temp2 = 0.d0
         temp3 = 0.d0
         ! A is composed of sum_nu(T_nu^T N_nu^-1 T_nu)

         ! Step one is to solve temp1 = (T_nu)(vec)
         do m = 1, size(self%joint_comp)
            do n = 1, self%ncomp
               if (trim(self%joint_comp(m)) == trim(self%fg_label(n))) then
                  ! First the multiply npix values of x by the npix diagonals of T_nu
                  !$OMP PARALLEL PRIVATE(i)
                  !$OMP DO SCHEDULE(static)
                  do i = 1, x
                     if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                     temp1(i) = compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n)*vec(i)
                     if (self%joint_pol) temp1(x+i) = compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n+1)*vec(x+i)
                  end do
                  !$OMP END DO
                  !$OMP END PARALLEL
               end if
            end do
            do n = 1, self%ntemp
               if (trim(self%joint_comp(m)) == trim(self%temp_label(n))) then
                  ! Then add vec(i)*template(i) to each temp1(i)
                  if (self%temp_corr(n,j)) then
                     !!$OMP PARALLEL PRIVATE(i)
                     !!$OMP DO SCHEDULE(static)
                     do i = 1, x
                        if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                        temp1(i) = temp1(i) + dat%temps(i-1,map_n,n)*vec(len+l)
                        if (self%joint_pol) temp1(x+i) = temp1(x+i) + dat%temps(i-1,map_n+1,n)*vec(len+l)
                     end do
                     !!$OMP END DO
                     !!$OMP END PARALLEL
                  end if
               end if
            end do
         end do

         ! Step two is to solve temp2 = (N^-1)(temp1)
         !$OMP PARALLEL PRIVATE(i)
         !$OMP DO SCHEDULE(static)
         do i = 1, x
            if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
            temp2(i) = temp1(i)/(dat%rms_map(i-1,map_n,j)**2.d0)
            if (self%joint_pol) temp2(x+i) = temp1(x+i)/(dat%rms_map(i-1,map_n+1,j)**2.d0)
         end do
         !$OMP END DO
         !$OMP END PARALLEL

         ! Step three is to solve temp3 = (T_nu^t)(temp2) - this part is complex/confusing...
         do m = 1, size(self%joint_comp)
            do n = 1, self%ncomp
               if (trim(self%joint_comp(m)) == trim(self%fg_label(n))) then
                  !$OMP PARALLEL PRIVATE(i)
                  !$OMP DO SCHEDULE(static)
                  do i = 1, x
                     if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                     temp3(i) = temp2(i)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n)
                     if (self%joint_pol) temp3(x+i) = temp2(x+i)*compute_spectrum(self,compo,n,self%dat_nu(j),i-1,map_n+1)
                  end do
                  !$OMP END DO
                  !$OMP END PARALLEL
               end if
            end do
            do n = 1, self%ntemp
               if (trim(self%joint_comp(m)) == trim(self%temp_label(n))) then
                  if (self%temp_corr(n,j)) then
                     !!$OMP PARALLEL PRIVATE(i)
                     !!$OMP DO SCHEDULE(static)
                     do i = 1, x
                        if (dat%masks(i-1,1) == 0.d0 .or. dat%masks(i-1,1) == missval) cycle
                        temp3(len+l) = temp3(len+l) + temp2(i)*dat%temps(i-1,map_n,n)
                        if (self%joint_pol) temp3(len+l) = temp3(len+l) + temp2(x+i)*dat%temps(i-1,map_n+1,n)
                     end do
                     !!$OMP END DO
                     !!$OMP END PARALLEL
                     l = l + 1
                  end if
               end if
            end do
         end do
         res = res + temp3
      end do

    end function return_Ax

    !-------------------------------------------------------------------------------------

    subroutine lower_tri_Ax(A,x,n)
      implicit none
      real(dp), dimension(:,:), intent(in)  :: A
      real(dp), dimension(:), intent(inout) :: x
      integer(i4b)                          :: i, j, k, n
      real(dp)                              :: temp
      
      do j = n,1,-1
         temp = 0.d0
         if (x(j) /= 0.d0) then
            temp = x(j)
            !$OMP PARALLEL PRIVATE(i)
            !$OMP DO SCHEDULE(static)
            do i = n, j+1, -1
               x(i) = x(i) + temp*a(i,j)
            end do
            !$OMP END DO
            !$OMP END PARALLEL
            x(j) = x(j)*a(j,j)
         end if
      end do
      
    end subroutine lower_tri_Ax
    
    function compute_ATA(mat) result(B)
      implicit none
      real(dp), dimension(:,:), intent(in)  :: mat
      integer(i4b)                          :: ii, ij, ik, m, n
      real(dp), allocatable, dimension(:,:) :: B
      real(dp)                              :: q
      
      ! Little function to compute the product of a matrix with its transpose
      ! saves time since the matrix is symmetric!
      
      m = size(mat(1,:))
      n = size(mat(:,1))
      
      allocate(B(m,m))
      ! Compute A^T A
      B(:,:) = 0.d0
      
      do ii = 1, m
         do ij = ii, m
            q = 0.d0
            do ik = 1, n
               q = q + mat(ik,ii)*mat(ik,ij)
            end do
            B(ii,ij) = q
            B(ij,ii) = q
         end do
      end do
      
    end function compute_ATA
    
    function compute_ATA_CSC(v,row_i,col_p) result(B)
      implicit none
      real(dp), dimension(:),     intent(in) :: v
      integer(i4b), dimension(:), intent(in) :: row_i, col_p
      integer(i4b)                           :: n, col, row, i, ii, k, ik, nnz
      real(dp), allocatable, dimension(:,:)  :: B
      
      n   = size(col_p)-1
      nnz = size(row_i)
      
      allocate(B(n,n))
      
      b = 0.d0
      
      !$OMP PARALLEL PRIVATE(i,k,ii,ik,row,col)
      !$OMP DO SCHEDULE(DYNAMIC)
      do i = 1, n+1
         do k = col_p(i), col_p(i+1)-1
            if (k > nnz) then
               write(*,*) 'k too big'
               !stop
               exit
            end if
            col = i
            do ii = 1, n+1
               row = ii
               do ik = col_p(ii), col_p(ii+1)-1
                  if (ik > nnz) then
                     !write(*,*) ik
                     !write(*,*) 'ik too big'
                     !stop
                  end if
                  if (row_i(k) == row_i(ik)) then
                     b(row,col) = b(row,col) + v(k)*v(ik)
                  end if
               end do
            end do
         end do
      end do
      !$OMP END DO
      !$OMP DO SCHEDULE(DYNAMIC)
      do i = 1, n
         do k = 1, n
            if (b(i,k) /= 0.d0) b(k,i) = b(i,k)
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
      
    end function compute_ATA_CSC
    
    function compute_ATA_CSR(v,col_i,row_p) result(B)
      implicit none
      real(dp), dimension(:),     intent(in) :: v
      integer(i4b), dimension(:), intent(in) :: col_i, row_p
      integer(i4b)                           :: n, col, row, i, ii, k, ik, nnz
      real(dp), allocatable, dimension(:,:)  :: B
    
      n   = size(row_p)-1
      nnz = size(col_i)
      
      allocate(B(n,n))
      
      b = 0.d0
      
      !$OMP PARALLEL PRIVATE(i,k,ii,ik,row,col)
      !$OMP DO SCHEDULE(STATIC)
      do i = 1, n+1
         do k = row_p(i), row_p(i+1)-1
            row = i
            do ii = 1, n+1
               col = ii
               do ik = row_p(ii), row_p(ii+1)
                  if (col_i(k) == col_i(ik)) then
                     b(row,col) = b(row,col) + v(k)*v(ik)
                  end if
               end do
            end do
         end do
      end do
      !$OMP END DO
      !$OMP DO SCHEDULE(STATIC)
      do i = 1, n
         do k = 1, n
            if (b(i,k) /= 0.d0) b(k,i) = b(i,k)
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
      
    end function compute_ATA_CSR
    
    subroutine A_to_CSC(A,col_p,row_i,v)
      implicit none
      real(dp), dimension(:,:),   intent(in)  :: A
      real(dp), dimension(:),     intent(out) :: v
      integer(i4b), dimension(:), intent(out) :: col_p, row_i
      integer(i4b)               :: vi, ci, ri, co, n, m, i ,j

      m = size(A(1,:))
      n = size(A(:,1))

      vi  = 1
      ci  = 1
      ri  = 1
      co  = 1
      col_p(ci) = 1
      ci  = ci + 1

      do i = 1, n
         do j = 1, m
            if (A(j,i) /= 0.d0) then
               v(vi) = A(j,i)
               row_i(ri) = j
               co = co + 1
               vi = vi + 1
               ri = ri + 1
            end if
         end do
         col_p(ci) = co
         ci = ci + 1
      end do

    end subroutine A_to_CSC

    subroutine A_to_CSR(A,row_p,col_i,v)
      implicit none
      real(dp), dimension(:,:),   intent(in)  :: A
      real(dp), dimension(:),     intent(out) :: v
      integer(i4b), dimension(:), intent(out) :: row_p, col_i
      integer(i4b)                            :: vi, ci, ri, ro, n, m, i, j

      m = size(A(1,:))
      n = size(A(:,1))

      vi  = 1
      ci  = 1
      ri  = 1
      ro  = 1
      row_p(ri) = 1
      ri  = ri + 1

      do i = 1, n
         do j = 1, m
            if (A(i,j) /= 0.d0) then
               v(vi)       = A(i,j)
               col_i(ci) = j
               ro  = ro + 1
               vi  = vi + 1
               ci  = ci + 1
            end if
         end do
         row_p(ri) = ro
         ri = ri + 1
      end do
  
    end subroutine A_to_CSR

    function Ax_csr(row_p,col_i,val,x) result(res)
      implicit none
      real(dp), dimension(:), intent(in)     :: val, x
      integer(i4b), dimension(:), intent(in) :: row_p, col_i
      real(dp), allocatable, dimension(:)    :: res

      integer(i4b)                           :: i, j, k, n, nnz, nval, mtype
      integer(i4b)                           :: nrows, offset, dest, source, end

      n   = size(row_p)-1
      nnz = size(col_i)


      allocate(res(n))
      res = 0.d0

      !$OMP PARALLEL PRIVATE(i,k)
      !$OMP DO SCHEDULE(GUIDED)
      do i = 1, n
         res(i) = 0.d0
         do k = row_p(i), row_p(i+1)-1
            res(i) = res(i) + val(k)*x(col_i(k))
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end function Ax_csr

    function Ax_csc(col_p,row_i,val,x) result(res)
      ! Does not work currently...I don't think.
      implicit none
      real(dp), dimension(:), intent(in)     :: val, x
      integer(i4b), dimension(:), intent(in) :: col_p, row_i
      real(dp), allocatable, dimension(:)    :: res

      integer(i4b)                           :: i, k, n

      n = size(col_p)-1

      allocate(res(n))

      res = 0.d0

      do i = 1, n
         do k = col_p(i), col_p(i+1)-1
            res(row_i(k)) = res(row_i(k)) + val(k)*x(i)
         end do
      end do
    end function Ax_csc

end module dang_linalg_mod
