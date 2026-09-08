program test_math_tools_linear
  use healpix_types
  use math_tools, only: tridag, cholesky_decompose, cholesky_solve, invert_matrix_dp
  implicit none

  integer(i4b) :: failures

  failures = 0

  call test_tridag_solver(failures)
  call test_cholesky_decompose_spd(failures)
  call test_cholesky_decompose_non_spd(failures)
  call test_cholesky_solve_system(failures)
  call test_invert_matrix_dp_general(failures)
  call test_invert_matrix_dp_cholesky(failures)
  call test_invert_matrix_dp_singular_status(failures)
  call test_invert_matrix_dp_ill_conditioned(failures)

  if (failures /= 0) then
    write(*, '(A,I0)') "test_math_tools_linear failures: ", failures
    stop 1
  end if

  write(*, '(A)') "test_math_tools_linear passed"

contains

  subroutine test_tridag_solver(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2) :: a, c
    real(dp), dimension(3) :: b, r, u

    a = [-1.0_dp, -1.0_dp]
    b = [2.0_dp, 2.0_dp, 2.0_dp]
    c = [-1.0_dp, -1.0_dp]
    r = [1.0_dp, 0.0_dp, 1.0_dp]

    call tridag(a, b, c, r, u)

    call expect_close_vec(u, [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-12_dp, "tridag", failures)
  end subroutine test_tridag_solver

  subroutine test_cholesky_decompose_spd(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A, L, recon
    integer(i4b) :: ierr

    A = reshape([4.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], [2, 2])

    call cholesky_decompose(A, L, ierr)
    call expect_equal_i4(ierr, 0_i4b, "cholesky ierr spd", failures)
    call expect_close_scalar(L(1, 2), 0.0_dp, 1.0e-12_dp, "cholesky upper(1,2)", failures)

    recon = matmul(L, transpose(L))
    call expect_close_mat(recon, A, 1.0e-12_dp, "cholesky reconstruction", failures)
  end subroutine test_cholesky_decompose_spd

  subroutine test_cholesky_decompose_non_spd(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A, L
    integer(i4b) :: ierr

    A = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2])

    call cholesky_decompose(A, L, ierr)
    if (ierr == 0) then
      failures = failures + 1
      write(*, '(A)') "FAIL: cholesky non-spd should return non-zero ierr"
    end if
  end subroutine test_cholesky_decompose_non_spd

  subroutine test_cholesky_solve_system(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A, L
    real(dp), dimension(2) :: b, x
    integer(i4b) :: ierr

    A = reshape([4.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], [2, 2])
    b = [6.0_dp, 7.0_dp]

    call cholesky_decompose(A, L, ierr)
    call expect_equal_i4(ierr, 0_i4b, "cholesky for solve ierr", failures)
    call cholesky_solve(L, b, x)

    call expect_close_vec(x, [0.5_dp, 2.0_dp], 1.0e-12_dp, "cholesky_solve", failures)
  end subroutine test_cholesky_solve_system

  subroutine test_invert_matrix_dp_general(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A, Ainv, ident
    integer(i4b) :: status

    A = reshape([4.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], [2, 2])
    Ainv = A
    call invert_matrix_dp(Ainv, status=status)
    call expect_equal_i4(status, 0_i4b, "invert_matrix_dp general status", failures)

    ident = matmul(A, Ainv)
    call expect_close_mat(ident, eye2(), 1.0e-10_dp, "invert_matrix_dp general identity", failures)
    call expect_relative_residual(ident, eye2(), 1.0e-10_dp, "invert_matrix_dp general residual", failures)
  end subroutine test_invert_matrix_dp_general

  subroutine test_invert_matrix_dp_cholesky(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A, Ainv, ident
    integer(i4b) :: status

    A = reshape([4.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], [2, 2])
    Ainv = A
    call invert_matrix_dp(Ainv, cholesky=.true., status=status)
    call expect_equal_i4(status, 0_i4b, "invert_matrix_dp cholesky status", failures)

    ident = matmul(A, Ainv)
    call expect_close_mat(ident, eye2(), 1.0e-10_dp, "invert_matrix_dp cholesky identity", failures)
    call expect_relative_residual(ident, eye2(), 1.0e-10_dp, "invert_matrix_dp cholesky residual", failures)
  end subroutine test_invert_matrix_dp_cholesky

  subroutine test_invert_matrix_dp_singular_status(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A
    integer(i4b) :: status

    A = reshape([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], [2, 2])
    call invert_matrix_dp(A, status=status)
    if (status == 0) then
      failures = failures + 1
      write(*, '(A)') "FAIL: invert_matrix_dp singular matrix should report non-zero status"
    end if
  end subroutine test_invert_matrix_dp_singular_status

  subroutine test_invert_matrix_dp_ill_conditioned(failures)
    integer(i4b), intent(inout) :: failures
    integer(i4b), parameter :: n = 5
    real(dp), dimension(n, n) :: A, Ainv_lu, Ainv_chol, ident_lu, ident_chol
    integer(i4b) :: status

    call fill_hilbert(A)

    Ainv_lu = A
    call invert_matrix_dp(Ainv_lu, status=status)
    call expect_equal_i4(status, 0_i4b, "invert_matrix_dp hilbert LU status", failures)
    ident_lu = matmul(A, Ainv_lu)
    call expect_close_mat(ident_lu, eye(n), 5.0e-8_dp, "invert_matrix_dp hilbert LU identity", failures)
    call expect_relative_residual(ident_lu, eye(n), 5.0e-8_dp, "invert_matrix_dp hilbert LU residual", failures)

    Ainv_chol = A
    call invert_matrix_dp(Ainv_chol, cholesky=.true., status=status)
    call expect_equal_i4(status, 0_i4b, "invert_matrix_dp hilbert Cholesky status", failures)
    ident_chol = matmul(A, Ainv_chol)
    call expect_close_mat(ident_chol, eye(n), 5.0e-8_dp, "invert_matrix_dp hilbert Cholesky identity", failures)
    call expect_relative_residual(ident_chol, eye(n), 5.0e-8_dp, "invert_matrix_dp hilbert Cholesky residual", failures)
  end subroutine test_invert_matrix_dp_ill_conditioned

  subroutine fill_hilbert(A)
    real(dp), dimension(:, :), intent(out) :: A
    integer(i4b) :: i, j

    do i = 1, size(A, 1)
      do j = 1, size(A, 2)
        A(i, j) = 1.0_dp / real(i + j - 1, dp)
      end do
    end do
  end subroutine fill_hilbert

  function eye2() result(E)
    real(dp), dimension(2, 2) :: E

    E = 0.0_dp
    E(1, 1) = 1.0_dp
    E(2, 2) = 1.0_dp
  end function eye2

  function eye(n) result(E)
    integer(i4b), intent(in) :: n
    real(dp), dimension(n, n) :: E
    integer(i4b) :: i

    E = 0.0_dp
    do i = 1, n
      E(i, i) = 1.0_dp
    end do
  end function eye

  subroutine expect_equal_i4(actual, expected, label, failures)
    integer(i4b), intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (actual /= expected) then
      failures = failures + 1
      write(*, '(A,A,A,I0,A,I0)') "FAIL: ", trim(label), " expected=", expected, " actual=", actual
    end if
  end subroutine expect_equal_i4

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

  subroutine expect_relative_residual(actual, expected, tol, label, failures)
    real(dp), dimension(:, :), intent(in) :: actual, expected
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures
    real(dp) :: resid, denom

    denom = max(sqrt(sum(expected**2)), 1.0e-30_dp)
    resid = sqrt(sum((actual - expected)**2)) / denom

    if (resid > tol) then
      failures = failures + 1
      write(*, '(A,A,A,ES12.5,A,ES12.5)') "FAIL: ", trim(label), " residual=", resid, " tol=", tol
    end if
  end subroutine expect_relative_residual

end program test_math_tools_linear
