program test_dang_linalg_mod
  use healpix_types
  use dang_linalg_mod, only: compute_ATA, LUDecomp, forward_sub, backward_sub
  use dang_linalg_mod, only: cholesky_decomp, invert_tri, A_to_CSR, Ax_csr
  use dang_linalg_mod, only: A_to_CSC, Ax_csc, lower_tri_Ax
  implicit none

  integer(i4b) :: failures

  failures = 0

  call test_compute_ata_against_definition(failures)
  call test_lu_forward_backward_solve(failures)
  call test_cholesky_factorization_identity(failures)
  call test_invert_triangular_identity(failures)
  call test_sparse_matvec_matches_dense(failures)
  call test_lower_tri_ax_matches_formula(failures)

  if (failures /= 0) then
    write(*, '(A,I0)') "test_dang_linalg_mod failures: ", failures
    stop 1
  end if

  write(*, '(A)') "test_dang_linalg_mod passed"

contains

  subroutine test_compute_ata_against_definition(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(3, 2) :: A
    real(dp), allocatable, dimension(:,:) :: B
    real(dp), dimension(2, 2) :: B_ref
    real(dp), dimension(2) :: x
    real(dp) :: quad

    A = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, -1.0_dp, 2.0_dp], [3, 2])
    B = compute_ATA(A)
    B_ref = matmul(transpose(A), A)
    call expect_close_mat(B, B_ref, 1.0e-12_dp, "compute_ATA equals A^T A", failures)

    x = [0.75_dp, -0.25_dp]
    quad = dot_product(x, matmul(B, x))
    call expect_close_scalar(quad, sum(matmul(A, x)**2), 1.0e-12_dp, "x^T(A^TA)x = ||Ax||^2", failures)
    if (quad < -1.0e-12_dp) then
      failures = failures + 1
      write(*, '(A)') "FAIL: compute_ATA produced non-PSD quadratic form"
    end if
  end subroutine test_compute_ata_against_definition

  subroutine test_lu_forward_backward_solve(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(3, 3) :: A, L, U
    real(dp), dimension(3) :: b, y, x

    A = reshape([4.0_dp, 2.0_dp, 0.0_dp, 2.0_dp, 5.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp], [3, 3])
    b = [2.0_dp, -1.0_dp, 4.0_dp]

    call LUDecomp(A, L, U, 3_i4b)
    call forward_sub(L, y, b)
    call backward_sub(U, x, y)

    call expect_close_vec(matmul(A, x), b, 1.0e-10_dp, "LU solve residual", failures)
    call expect_close_mat(matmul(L, U), A, 1.0e-10_dp, "A = L*U", failures)
  end subroutine test_lu_forward_backward_solve

  subroutine test_cholesky_factorization_identity(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(3, 3) :: A, L

    A = reshape([6.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 5.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 4.0_dp], [3, 3])
    call cholesky_decomp(A, L)
    call zero_upper_triangle(L)
    call expect_close_mat(matmul(L, transpose(L)), A, 1.0e-10_dp, "cholesky reconstruction", failures)
  end subroutine test_cholesky_factorization_identity

  subroutine test_invert_triangular_identity(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(3, 3) :: L, Linv

    L = reshape([2.0_dp, -1.0_dp, 4.0_dp, 0.0_dp, 3.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 1.5_dp], [3, 3])
    call invert_tri(L, 'L', Linv)
    call zero_upper_triangle(Linv)
    call expect_close_mat(matmul(L, Linv), eye3(), 1.0e-10_dp, "invert_tri identity", failures)
  end subroutine test_invert_triangular_identity

  subroutine test_sparse_matvec_matches_dense(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(3, 3) :: A
    real(dp), dimension(3) :: x, dense
    real(dp), allocatable, dimension(:) :: vals, y_csr, y_csc
    integer(i4b), allocatable, dimension(:) :: row_p, col_i, col_p, row_i
    integer(i4b) :: nnz

    A = reshape([2.0_dp, 0.0_dp, -1.0_dp, 1.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, -2.0_dp, 4.0_dp], [3, 3])
    x = [1.0_dp, -1.0_dp, 2.0_dp]
    dense = matmul(A, x)
    nnz = count(abs(A) > 0.0_dp)

    allocate(vals(nnz), row_p(4), col_i(nnz), col_p(4), row_i(nnz))

    call A_to_CSR(A, row_p, col_i, vals)
    y_csr = Ax_csr(row_p, col_i, vals, x)
    call expect_close_vec(y_csr, dense, 1.0e-12_dp, "CSR Ax", failures)

    call A_to_CSC(A, col_p, row_i, vals)
    y_csc = Ax_csc(col_p, row_i, vals, x)
    call expect_close_vec(y_csc, dense, 1.0e-12_dp, "CSC Ax", failures)

    deallocate(vals, row_p, col_i, col_p, row_i, y_csr, y_csc)
  end subroutine test_sparse_matvec_matches_dense

  subroutine test_lower_tri_ax_matches_formula(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(3, 3) :: L
    real(dp), dimension(3) :: x_work, x_ref

    L = reshape([2.0_dp, -1.0_dp, 4.0_dp, 0.0_dp, 3.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 1.5_dp], [3, 3])
    x_work = [0.5_dp, -2.0_dp, 1.0_dp]

    x_ref = matmul(L, x_work)
    call lower_tri_Ax(L, x_work, 3_i4b)
    call expect_close_vec(x_work, x_ref, 1.0e-12_dp, "lower_tri_Ax", failures)
  end subroutine test_lower_tri_ax_matches_formula

  function eye3() result(E)
    real(dp), dimension(3, 3) :: E
    E = 0.0_dp
    E(1, 1) = 1.0_dp
    E(2, 2) = 1.0_dp
    E(3, 3) = 1.0_dp
  end function eye3

  subroutine zero_upper_triangle(A)
    real(dp), dimension(:, :), intent(inout) :: A
    integer(i4b) :: i, j
    do i = 1, size(A, 1)
      do j = i + 1, size(A, 2)
        A(i, j) = 0.0_dp
      end do
    end do
  end subroutine zero_upper_triangle

  subroutine expect_close_scalar(actual, expected, tol, label, failures)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures
    if (abs(actual - expected) > tol) then
      failures = failures + 1
      write(*, '(A,A,A,ES12.5,A,ES12.5)') "FAIL: ", trim(label), " expected=", expected, " actual=", actual
    end if
  end subroutine expect_close_scalar

  subroutine expect_close_vec(actual, expected, tol, label, failures)
    real(dp), dimension(:), intent(in) :: actual, expected
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures
    if (any(abs(actual - expected) > tol)) then
      failures = failures + 1
      write(*, '(A,A)') "FAIL: ", trim(label)
    end if
  end subroutine expect_close_vec

  subroutine expect_close_mat(actual, expected, tol, label, failures)
    real(dp), dimension(:, :), intent(in) :: actual, expected
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures
    if (any(abs(actual - expected) > tol)) then
      failures = failures + 1
      write(*, '(A,A)') "FAIL: ", trim(label)
    end if
  end subroutine expect_close_mat

end program test_dang_linalg_mod
