module linalg_mod
    use healpix_types
    use pix_tools
    use fitstools
    use udgrade_nr
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

    subroutine compute_cg(A,x,b,n)
        
        ! Implementation of the canned algorithm (B2) outlined in Jonathan Richard Shewuck (1994)
        ! "An introduction to the Conjugate Gradient Method Without the Agonizing Pain"

        implicit none

        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(:), intent(in)   :: b
        real(dp), dimension(:), intent(out)  :: x
        integer(i4b), intent(in)             :: n
        real(dp), allocatable, dimension(:)  :: r, q, d
        real(dp)                             :: epsil, alpha, beta, delta_0
        real(dp)                             :: delta_old, delta_new
        integer(i4b)                         :: i_max, i

        allocate(r(n),q(n),d(n))

        x(:) = 0.0d0
        i_max = 10

        i = 0
        epsil = 1.0d-16

        r = b - matmul(A,x)
        d = r
        delta_new = sum(r*r)
        delta_0   = delta_new

        do while( (i .lt. i_max) .and. (delta_new .gt. (epsil**2)*delta_0))
            q = matmul(A,d)
            alpha = delta_new/(sum(d*q))
            x = x + alpha*d
            if (mod(i,50) == 0) then
                r = b - matmul(A,x)
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
end module linalg_mod
