module linalg_mod
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
    use utility_mod
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

    subroutine compute_cg(A,x,b,n,nnz_a,iters,converge)
        
        ! Implementation of the canned algorithm (B2) outlined in Jonathan Richard Shewuck (1994)
        ! "An introduction to the Conjugate Gradient Method Without the Agonizing Pain"

        implicit none
        real(dp), dimension(:,:), intent(in)  :: A
        real(dp), dimension(:), intent(in)    :: b
        real(dp), dimension(:), intent(out)   :: x
        integer(i4b), intent(in)              :: n,iters
        real(dp), intent(in)                  :: converge
        integer(i4b), optional, intent(in)    :: nnz_a
        real(dp), allocatable, dimension(:)   :: r, q, d
        real(dp)                              :: epsil, alpha, beta, delta_0
        real(dp)                              :: delta_old, delta_new, t3, t4, t5, t6
        integer(i4b)                          :: i_max, i, j

        real(dp),allocatable,dimension(:)     :: v
        integer(i4b),allocatable,dimension(:) :: rp, ci

        t5 = mpi_wtime()
        allocate(r(n),q(n),d(n))

        if(present(nnz_a)) then
           allocate(v(nnz_a),rp(n+1),ci(nnz_a))
           call A_to_CSR(A,rp,ci,v)
        end if

        x(:) = 0.0d0
        i_max = iters

        i = 0
        epsil = 1.0d-16

        if (present(nnz_a)) then
           r = b - Ax_CSR(rp,ci,v,x)
        else
           r = b - matmul(A,x)
        end if
        d = r
        delta_new = sum(r*r)
        delta_0   = delta_new
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
           delta_old = delta_new
           delta_new = sum(r*r)
           beta = delta_new/delta_old
           d = r + beta*d
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
        if(present(nnz_a)) then
           deallocate(v)
           deallocate(rp)
           deallocate(ci)
        end if

        deallocate(r)
        deallocate(q)
        deallocate(d)

    end subroutine compute_cg

    subroutine compute_cg_vec(x,b,n,param,datas,compos)
        
      ! Implementation of the canned algorithm (B2) outlined in Jonathan Richard Shewuck (1994)
      ! "An introduction to the Conjugate Gradient Method Without the Agonizing Pain"

      implicit none
      type(params)                          :: param
      type(data)                            :: datas
      type(component)                       :: compos
      real(dp), dimension(:), intent(in)    :: b
      real(dp), dimension(:), intent(out)   :: x
      integer(i4b), intent(in)              :: n
      real(dp), allocatable, dimension(:)   :: r, q, d, s
      real(dp), allocatable, dimension(:)   :: M, w
      real(dp)                              :: converge
      real(dp)                              :: alpha, beta, delta_0
      real(dp)                              :: delta_old, delta_new, t3, t4, t5, t6
      integer(i4b)                          :: i_max, i, j

      t5 = mpi_wtime()
      allocate(r(n),q(n),d(n),M(n),s(n),w(n))

      x(:) = 0.0d0
      i_max    = param%cg_iter
      converge = param%cg_converge

      i = 0
      call multiply_with_A(param, datas, compos, w, 5, 2)

      r = b - return_Ax(param, datas, compos, x, 5, 2)
      d = r
      delta_new = sum(r*r)
      delta_0   = delta_new
      do while( (i .lt. i_max) .and. (delta_new .gt. converge))!(epsil**2)*delta_0))
         t3 = mpi_wtime()
         q = return_Ax(param, datas, compos, d, 5, 2)
         alpha = delta_new/(sum(d*q))
         x = x + alpha*d
         if (mod(i,50) == 0) then
            r = b - return_Ax(param, datas, compos, x, 5, 2)
         else
            r = r - alpha*q
         end if
         delta_old = delta_new
         delta_new = sum(r*r)
         beta = delta_new/delta_old
         d = r + beta*d
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

    end subroutine compute_cg_vec


    subroutine compute_cg_precond(x,b,n,param,datas,compos)
        
      ! Implementation of the canned algorithm (B2) outlined in Jonathan Richard Shewuck (1994)
      ! "An introduction to the Conjugate Gradient Method Without the Agonizing Pain"
      
      ! HAS NOT BEEN COMPLETED YET
      
      ! Note - the preconditioner form here is different that is normally seen
      ! Normally a preconditioner matrix M is inverted and applied. M here is the 
      ! inverse of the preconditioner we want (inverse of the diagonal version of A).

      implicit none
      type(params)                          :: param
      type(data)                            :: datas
      type(component)                       :: compos
      real(dp), dimension(:), intent(in)    :: b
      real(dp), dimension(:), intent(out)   :: x
      integer(i4b), intent(in)              :: n
      real(dp), allocatable, dimension(:)   :: r, q, d, s
      real(dp), allocatable, dimension(:)   :: M, w
      real(dp)                              :: converge
      real(dp)                              :: epsil, alpha, beta, delta_0
      real(dp)                              :: delta_old, delta_new, t3, t4, t5, t6
      integer(i4b)                          :: i_max, i, j

      t5 = mpi_wtime()
      allocate(r(n),q(n),d(n),M(n),s(n),w(n))

      x(:)     = 0.0d0
      M(:)     = 1.d0
      w(:)     = 0.d0
      i_max    = param%cg_iter
      converge = param%cg_converge

      !call multiply_with_A(param, datas, compos, M, 5, 2)

      !write(*,*) 'Define preconditioner'
      !do i = 1, n
      !   M(i) = 1.d0/M(i)
      !end do

      i = 0
      epsil = 1.0d-16

      r = b

      !!$OMP PARALLEL PRIVATE(i)
      !!$OMP DO SCHEDULE(static)
      !do i = 1, n
      !   d(i) = M(i)*r(i)
      !end do
      !!$OMP END DO
      !!$OMP END PARALLEL
      !delta_new = sum(r*d)
      
      d = r
      delta_new = sum(r*r)
      delta_0   = delta_new
      t3 = mpi_wtime()
      do while( (i .lt. i_max) .and. (delta_new .gt. converge))!(epsil**2)*delta_0))
         t3 = mpi_wtime()
         q = d
         call multiply_with_A(param, datas, compos, q, 5, 2)
         alpha = delta_new/(sum(d*q))
         write(*,*) alpha
         x = x + alpha*d
         if (mod(i,50) == 0) then
            w = x
            call multiply_with_A(param, datas, compos, w, 5, 2)
            r = b - w
         else
            r = r - alpha*q
         end if
         !!$OMP PARALLEL PRIVATE(i)
         !!$OMP DO SCHEDULE(static)
         !do j = 1, n
         !   s(j) = M(j)*r(j)
         !end do
         !!$OMP END DO
         !!$OMP END PARALLEL
         delta_old = delta_new
         delta_new = sum(r*r)
         beta = delta_new/delta_old
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
          if (i == k) write(*,*) i, b(i,i)
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

    subroutine multiply_with_A(para, dat, compo, vec, nbands, map_n)
      implicit none
      type(params)                           :: para
      type(data)                             :: dat
      type(component)                        :: compo

      real(dp), dimension(:), intent(inout)  :: vec
      real(dp), allocatable, dimension(:)    :: v_temp, v_temp2, v_temp3
      real(dp), allocatable, dimension(:)    :: vech
      real(dp)                               :: x
      integer(i4b), intent(in)               :: nbands, map_n
      integer(i4b)                           :: n, i, j, k, len, l

      x = dat%npix
      n = size(vec)
      len = 2*x

      allocate(vech(n))
      allocate(v_temp(n))
      allocate(v_temp2(len))
      allocate(v_temp3(n))

      vech = vec

      v_temp3 = 0.d0

      write(*,*) n, len

      l = 1
      do j = 1, nbands
         v_temp  = 0.d0
         v_temp2 = 0.d0
         ! A is composed of sum_nu(T_nu^T N_nu^-1 T_nu)
         
         ! first multiply by T_nu
         !$OMP PARALLEL PRIVATE(i,k)
         !$OMP DO SCHEDULE(static)
         do i = 1, x
            v_temp2(i)   = vech(i)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n)
            v_temp2(x+i) = vech(x+i)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n+1)
            if (para%temp_corr(1,j)) then
               if (i == 1) write(*,*) 'truuu', j, l
               v_temp2(i)   = v_temp2(i)   + vech(len+l)*dat%temps(i-1,map_n,1)
               v_temp2(x+i) = v_temp2(x+i) + vech(len+l)*dat%temps(i-1,map_n+1,1)
            end if
            v_temp2(i)   = v_temp2(i)/(dat%rms_map(i-1,map_n,j))**2.d0
            v_temp2(x+i) = v_temp2(x+i)/(dat%rms_map(i-1,map_n+1,j))**2.d0
            v_temp(i)    = v_temp(i)   + v_temp2(i)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n)
            v_temp(x+i)  = v_temp(x+i) + v_temp2(x+i)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n+1)
         end do
         !$OMP END DO

         !do i = len+1, n
         if (para%temp_corr(1,j)) then
            !$OMP DO SCHEDULE(static)
            do k = 1, x
               v_temp(len+l) = v_temp(len+l) + dat%temps(k-1,map_n,1)*v_temp2(k)
               v_temp(len+l) = v_temp(len+l) + dat%temps(k-1,map_n+1,1)*v_temp2(k+x)
            end do
            !$OMP END DO
            l = l+1
         end if
         !end do
         !$OMP END PARALLEL
         !write(*,*) v_temp(len), v_temp2(len)
         v_temp3 = v_temp3 + v_temp
      end do

      write(*,*) v_temp2(len:)
      stop
      
      vec = v_temp3

    end subroutine multiply_with_A

    function return_Ax(para, dat, compo, vec, nbands, map_n) result(res)
      implicit none
      type(params)                           :: para
      type(data)                             :: dat
      type(component)                        :: compo

      real(dp), dimension(:), intent(in)     :: vec
      real(dp), allocatable, dimension(:)    :: v_temp, v_temp2, v_temp3
      real(dp), allocatable, dimension(:)    :: vech, res
      real(dp)                               :: x
      integer(i4b), intent(in)               :: nbands, map_n
      integer(i4b)                           :: n, i, j, k, len, l

      x = dat%npix
      n = size(vec)
      len = 2*x

      allocate(vech(n))
      allocate(v_temp(n))
      allocate(v_temp2(len))
      allocate(v_temp3(n))
      allocate(res(n))

      vech = vec

      v_temp3 = 0.d0

      l = 1
      do j = 1, nbands
         v_temp  = 0.d0
         v_temp2 = 0.d0
         ! A is composed of sum_nu(T_nu^T N_nu^-1 T_nu)
         
         ! first multiply by T_nu
         !$OMP PARALLEL PRIVATE(i,k)
         !$OMP DO SCHEDULE(static)
         do i = 1, x
            v_temp2(i)   = vech(i)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n)
            v_temp2(x+i) = vech(x+i)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n+1)
            if (para%temp_corr(1,j)) then
               v_temp2(i)   = v_temp2(i)   + vech(len+l)*dat%temps(i-1,map_n,1)
               v_temp2(x+i) = v_temp2(x+i) + vech(len+l)*dat%temps(i-1,map_n+1,1)
            end if
         end do
         !$OMP END DO

         ! Then multiply by N_nu^-1
         !$OMP DO SCHEDULE(static)
         do i = 1, x
            v_temp2(i)   = v_temp2(i)/(dat%rms_map(i-1,map_n,j))**2.d0
            v_temp2(x+i) = v_temp2(x+i)/(dat%rms_map(i-1,map_n+1,j))**2.d0
         end do
         !$OMP END DO

         ! Then multiply by T_nu^T
         !$OMP DO SCHEDULE(static)
         do i = 1, x
            v_temp(i)   = v_temp(i)   + v_temp2(i)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n)
            v_temp(x+i) = v_temp(x+i) + v_temp2(x+i)*compute_spectrum(para,compo,1,para%dat_nu(j),i-1,map_n+1)
         end do
         !$OMP END DO

         !do i = len+1, n
         if (para%temp_corr(1,j)) then
            !$OMP DO SCHEDULE(static)
            do k = 1, x
               v_temp(len+l) = v_temp(len+l) + dat%temps(k-1,map_n,1)*v_temp2(k)
               v_temp(len+l) = v_temp(len+l) + dat%temps(k-1,map_n+1,1)*v_temp2(k+x)
            end do
            !$OMP END DO
            l = l+1
         end if
         !end do
         !$OMP END PARALLEL
         v_temp3 = v_temp3 + v_temp
      end do
      
      res = v_temp3

    end function return_Ax


end module linalg_mod
