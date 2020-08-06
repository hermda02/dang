module linalg_mod
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
    use init_mod
    use utility_mod
    implicit none

contains

    function inv(A) result(Ainv)
        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(size(A,1),size(A,2)) :: Ainv
    
        real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info
    
        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI
    
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
    
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)
    
        if (info /= 0) then
        stop 'Matrix is numerically singular!'
        end if
    
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)
    
        if (info /= 0) then
        stop 'Matrix inversion failed!'
        end if
    end function inv

    subroutine cholesky_decomp(mat,low,n)
        implicit none
        real(dp), dimension(:,:), intent(in)  :: mat
        real(dp), dimension(:,:), intent(out) :: low
        integer(i4b),             intent(in)  :: n
        integer(i4b)                          :: ip, i, j
        real(dp)                              :: s

        low(:,:)   = 0.d0

        do i = 1, n
        low(i,i) = sqrt(mat(i,i) - dot_product(low(i,1:i-1),low(i,1:i-1)) )
        do j = i+1, n
            low(j,i) = (mat(j,i) - dot_product(low(j,1:i-1),low(i,1:i-1)))/low(i,i)
        end do
        end do

        ! This code would be used for the LDU decomp:
        ! -------------------------------------------
        ! do i = 1, n
        !     mat_s(i,i) = low(i,i)
        ! end do

        ! low  = matmul(low,inv(mat_s))
        ! diag = mat_s**2

        ! deallocate(mat_s)
        ! -------------------------------------------

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
        integer(i4b)                          :: i, n, m

        n = size(bf)

        do i = 1, n
            xf(i) = bf(i)
            do m = 1, i-1
                xf(i) = xf(i)-l(i,m)*xf(m) 
            end do
            xf(i) = xf(i)/L(i,i)
        end do

    end subroutine forward_sub

    subroutine backward_sub(U,xb,bb)
        ! Backward substitution to solve the matrix equation Ux=b
        implicit none
        real(dp), dimension(:,:), intent(in)  :: U
        real(dp), dimension(:),   intent(in)  :: bb
        real(dp), dimension(:),   intent(out) :: xb
        integer(i4b)                          :: n, m, i

        n = size(bb)
        xb(n) = bb(n)/U(n,n)

        do i = n-1, 1, -1
            xb(i) = bb(i)
            do m = n, i+1, -1
                xb(i) = xb(i) - U(i,m)*xb(m)
            end do
            xb(i) = xb(i)/U(i,i)
        end do

    end subroutine backward_sub

    subroutine compute_cg(A,x,b,n,nnz_a)
        
        ! Implementation of the canned algorithm (B2) outlined in Jonathan Richard Shewuck (1994)
        ! "An introduction to the Conjugate Gradient Method Without the Agonizing Pain"

        implicit none

        real(dp), dimension(:,:), intent(in)  :: A
        real(dp), dimension(:), intent(in)    :: b
        real(dp), dimension(:), intent(out)   :: x
        integer(i4b), intent(in)              :: n
        integer(i4b), optional, intent(in)    :: nnz_a
        real(dp), allocatable, dimension(:)   :: r, q, d
        real(dp)                              :: epsil, alpha, beta, delta_0
        real(dp)                              :: delta_old, delta_new, t3, t4
        integer(i4b)                          :: i_max, i, j

        real(dp),allocatable,dimension(:)     :: v
        integer(i4b),allocatable,dimension(:) :: rp, ci


        allocate(r(n),q(n),d(n))

        if(present(nnz_a)) then
           allocate(v(nnz_a),rp(n+1),ci(nnz_a))
           call A_to_CSR(A,rp,ci,v)
        end if

        x(:) = 0.0d0
        i_max = 10

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
        t3 = mpi_wtime()
        do while( (i .lt. i_max) .and. (delta_new .gt. (epsil**2)*delta_0))
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
            i = i + 1
        end do

        deallocate(r)
        deallocate(q)
        deallocate(d)

    end subroutine compute_cg

    function mm_mpi(A,B)
      implicit none

      integer(i4b)   :: nworkers, source, dest, mtype
      integer(i4b)   :: cols, avecol, extra, offset
      integer(i4b)   :: m, n, p, ii, ij, ik
      integer status(mpi_status_size)

      real(dp), dimension(:,:), intent(in)  :: A
      real(dp), dimension(:,:), intent(in)  :: B
      real(dp), dimension(size(A(:,1)),size(B(1,:))) :: mm_mpi

      m = size(A(:,1)) !number of rows in A
      n = size(A(1,:)) !number of columns in A (rows of B)
      p = size(B(1,:)) !number of columns of B 

      if (size(A(1,:)) /= size(B(:,1))) then
         if (rank == master) then
            write(*,*) 'Matrix sizes do not conform!'
            stop
         end if
      end if

      nworkers = numprocs - 1

      !----------------------- MASTER TASK ------------------------
      if (rank == master) then

         ! Send out matrix information to the workers
         avecol = p/nworkers
         extra  = mod(p,nworkers)
         offset = 1
         mtype  = FROM_MASTER

         do dest=1,nworkers
            if (dest .le. extra) then
               cols = avecol + 1
            else
               cols = avecol
            end if
            write(*,*) '  sending ', cols, ' cols to task ', dest
            call mpi_send(offset, 1, MPI_INTEGER, dest, mtype, MPI_COMM_WORLD, ierr)
            call mpi_send(cols, 1, mpi_integer, dest, mtype, mpi_comm_world, ierr)
            call mpi_send(a, m*n, mpi_double_precision, dest, mtype, mpi_comm_world, ierr)
            call mpi_send(b(1,offset), cols*n, mpi_double_precision, dest, mtype, &
                      mpi_comm_world, ierr)
            offset = offset + cols
         end do

         ! Aaand receive results from the workers
         mtype = from_worker
         do ii=1, nworkers
            source = ii
            call mpi_recv(offset, 1, mpi_integer, source, mtype, mpi_comm_world, status, ierr)
            call mpi_recv(cols, 1, mpi_integer, source, mtype, mpi_comm_world, status, ierr)
            call mpi_recv(mm_mpi(1,offset), cols*m, mpi_double_precision, source, mtype, &
                      mpi_comm_world, ierr)
         end do
      end if

      ! ------------------- WORKER TASK ---------------------------
  
      if (rank > master) then
         ! Receive data from master
         mtype = from_master
         call mpi_recv(offset, 1, mpi_integer, master, mtype, mpi_comm_world, status, ierr)
         call mpi_recv(cols, 1, mpi_integer,  master, mtype, mpi_comm_world, status, ierr)
         call mpi_recv(a, m*n, mpi_double_precision, master, mtype, mpi_comm_world, status, ierr)
         call mpi_recv(b, cols*n, mpi_double_precision, master, mtype, mpi_comm_world, status, ierr)

         ! Multiply matrices
         do ik = 1, cols
            do ii = 1, m
               mm_mpi(ii,ik) = 0.0
               do ij=1, n
                  mm_mpi(ii,ik) = mm_mpi(ii,ik) + a(ii,ij)*b(ij,ik)
               end do
            end do
         end do
         
         mtype = from_worker
         ! and send the results back
         call mpi_send(offset, 1, mpi_integer, master, mtype, mpi_comm_world, ierr)
         call mpi_send(cols, 1, mpi_integer, master, mtype, mpi_comm_world, ierr)
         call mpi_send(mm_mpi, cols*m, mpi_double_precision, master, mtype, mpi_comm_world, ierr)
      end if
    end function mm_mpi

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
      integer(i4b)                           :: n, col, row, i, ii, k, ik
      real(dp), allocatable, dimension(:,:)  :: B

      n = size(col_p)-1

      allocate(B(n,n))

      b = 0.d0

      do i = 1, n+1
         do k = col_p(i), col_p(i+1)-1
            col = i
            do ii = 1, n+1
               row = ii
               do ik = col_p(ii), col_p(ii+1)-1
                  if (row_i(k) == row_i(ik)) then
                     b(row,col) = b(row,col) + v(k)*v(ik)
                  end if
               end do
            end do
         end do
      end do
      do i = 1, n
         do k = 1, n
            if (b(i,k) /= 0.d0) b(k,i) = b(i,k)
         end do
      end do

    end function compute_ATA_CSC

    function compute_ATA_CSR(v,col_i,row_p) result(B)
      implicit none
      real(dp), dimension(:),     intent(in) :: v
      integer(i4b), dimension(:), intent(in) :: col_i, row_p
      integer(i4b)                           :: n, col, row, i, ii, k, ik
      real(dp), allocatable, dimension(:,:)  :: B

      n = size(row_p)-1

      allocate(B(n,n))

      b = 0.d0

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
      do i = 1, n
         do k = 1, n
            if (b(i,k) /= 0.d0) b(k,i) = b(i,k)
         end do
      end do

    end function compute_ATA_CSR

    subroutine A_to_CSC(A,col_ptr,row_ind,v)
      implicit none
      real(dp), dimension(:,:),   intent(in)  :: A
      real(dp), dimension(:),     intent(out) :: v
      integer(i4b), dimension(:), intent(out) :: col_ptr, row_ind
      integer(i4b)               :: vi, ci, ri, co, n, m, i ,j

      m = size(A(1,:))
      n = size(A(:,1))

      vi  = 1
      ci  = 1
      ri  = 1
      co  = 1
      col_ptr(ci) = 1
      ci  = ci + 1

      do i = 1, n
         do j = 1, m
            if (A(j,i) /= 0.d0) then
               v(vi) = A(j,i)
               row_ind(ri) = j
               co = co + 1
               vi = vi + 1
               ri = ri + 1
            end if
         end do
         col_ptr(ci) = co
         ci = ci + 1
      end do

    end subroutine A_to_CSC

    subroutine A_to_CSR(A,row_ptr,col_ind,v)
      implicit none
      real(dp), dimension(:,:),   intent(in)  :: A
      real(dp), dimension(:),     intent(out) :: v
      integer(i4b), dimension(:), intent(out) :: row_ptr, col_ind
      integer(i4b)                            :: vi, ci, ri, ro, n, m, i, j

      m = size(A(1,:))
      n = size(A(:,1))

      vi  = 1
      ci  = 1
      ri  = 1
      ro  = 1
      row_ptr(ri) = 1
      ri  = ri + 1

      do i = 1, n
         do j = 1, m
            if (A(i,j) /= 0.d0) then
               v(vi)       = A(i,j)
               col_ind(ci) = j
               ro  = ro + 1
               vi  = vi + 1
               ci  = ci + 1
            end if
         end do
         row_ptr(ri) = ro
         ri = ri + 1
      end do
  
    end subroutine A_to_CSR

    function Ax_csr(rptr,cind,val,x) result(res)
      implicit none
      real(dp), dimension(:), intent(in)     :: val, x
      integer(i4b), dimension(:), intent(in) :: rptr, cind
      real(dp), allocatable, dimension(:)    :: res

      integer(i4b)                           :: i, k, n

      n = size(rptr)-1

      allocate(res(n))

      res = 0.d0

      do i = 1, n+1
         do k = rptr(i), rptr(i+1)-1
            res(i) = res(i) + val(k)*x(cind(k))
         end do
      end do

    end function Ax_csr

    function Ax_csc(cptr,rind,val,x) result(res)
      implicit none
      real(dp), dimension(:), intent(in)     :: val, x
      integer(i4b), dimension(:), intent(in) :: cptr, rind
      real(dp), allocatable, dimension(:)    :: res

      integer(i4b)                           :: i, k, n

      n = size(cptr)-1

      allocate(res(n))

      res = 0.d0

      do i = 1, n+1
         do k = cptr(i), cptr(i+1)-1
            res(i) = res(i) + val(k)*x(i)
         end do
      end do

    end function Ax_csc

end module linalg_mod
