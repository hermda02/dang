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

    subroutine create_inv_N_precond(self,dat,precond,map_n)
      !--------------------------------------------------------------
      ! Subroutine which takes in the preconditioner and spits out
      ! a preconditioner vector which is just the sum (over bands) of 
      ! the inverse noise covariance matrices.
      !--------------------------------------------------------------
      implicit none
      type(dang_params)         :: self
      type(dang_data)           :: dat
      integer(i4b), intent(in)  :: map_n

      real(dp), dimension(:), intent(inout) :: precond
      integer(i4b)                          :: leng, mode, i, j

      leng = size(precond)

      if (leng < 2*dat%npix) then
         mode = 1
      else if (leng > 2*dat%npix) then
         mode = 2
      else
         write(*,*) "Something wrong with the preconditioner"
         write(*,*) leng
         stop
      end if

      precond(:) = 0.d0

      if (mode == 1) then
         do i = 1, dat%npix
            do j = 1, nbands
               precond(i) = precond(i) + 1.d0/(dat%rms_map(i-1,map_n,j)**2.0)
            end do
         end do
         precond(dat%npix:) = 1.d0
      else if (mode == 2) then
         do i = 1, dat%npix
            do j = 1, nbands
               precond(i) = precond(i) + 1.d0/(dat%rms_map(i-1,2,j)**2.0)
               precond(i+dat%npix) = precond(i) + 1.d0/(dat%rms_map(i-1,3,j)**2.0)
            end do
         end do
         precond(2*dat%npix:) = 1.d0
      end if

    end subroutine create_inv_N_precond

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
