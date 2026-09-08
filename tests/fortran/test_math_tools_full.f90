program test_math_tools_full
  use healpix_types
  use math_tools, only: invert_matrix_dpc, invert_matrix_sp, invert_matrix_with_mask_dpc, &
    & invert_matrix_with_mask_dp, invert_matrix_with_mask_sp, get_eigen_decomposition, &
    & get_eigenvalues, eigen_decomp, compute_hermitian_root, solve_system, invert_singular_matrix, &
    & solve_system_real, comp_normalised_Plm, convert_sigma2frac, convert_fract2sigma_dp, &
    & convert_fract2sigma_sp, cholesky_decompose_single, matmul_symm, matmul_gen, MOV3, corr_erf
  implicit none

  integer(i4b) :: failures

  failures = 0

  call test_invert_matrix_sp_ok(failures)
  call test_invert_matrix_dpc_ok(failures)
  call test_invert_matrix_with_mask_real_and_complex(failures)
  call test_get_eigenvalues_and_vectors(failures)
  call test_eigen_decomp_wrapper(failures)
  call test_compute_hermitian_root(failures)
  call test_compute_hermitian_root_non_spd(failures)
  call test_solve_systems(failures)
  call test_invert_singular_matrix_threshold(failures)
  call test_probability_conversions(failures)
  call test_cholesky_decompose_single(failures)
  call test_matmul_wrappers(failures)
  call test_mov3(failures)
  call test_comp_normalised_plm_smoke(failures)

  if (failures /= 0) then
    write(*, '(A,I0)') 'test_math_tools_full failures: ', failures
    stop 1
  end if

  write(*, '(A)') 'test_math_tools_full passed'

contains

  subroutine test_invert_matrix_sp_ok(failures)
    integer(i4b), intent(inout) :: failures
    real(sp), dimension(2, 2) :: A, Ainv, ident

    A = reshape([2.0_sp, 1.0_sp, 1.0_sp, 2.0_sp], [2, 2])
    Ainv = A
    call invert_matrix_sp(Ainv)
    ident = matmul(A, Ainv)
    call expect_close_mat_sp(ident, eye2_sp(), 2.0e-5_sp, 'invert_matrix_sp identity', failures)
  end subroutine test_invert_matrix_sp_ok

  subroutine test_invert_matrix_dpc_ok(failures)
    integer(i4b), intent(inout) :: failures
    complex(dpc), dimension(2, 2) :: A, Ainv, ident

    A = reshape([cmplx(2.0_dp, 0.0_dp, dpc), cmplx(0.0_dp, 0.0_dp, dpc), &
      & cmplx(0.0_dp, 0.0_dp, dpc), cmplx(3.0_dp, 0.0_dp, dpc)], [2, 2])
    Ainv = A
    call invert_matrix_dpc(Ainv)
    ident = matmul(A, Ainv)
    call expect_close_mat_cplx(ident, eye2_cplx(), 1.0e-12_dp, 'invert_matrix_dpc identity', failures)
  end subroutine test_invert_matrix_dpc_ok

  subroutine test_invert_matrix_with_mask_real_and_complex(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: Ar
    real(sp), dimension(2, 2) :: As
    complex(dpc), dimension(2, 2) :: Ac

    Ar = reshape([0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], [2, 2])
    call invert_matrix_with_mask_dp(Ar)
    call expect_close_scalar_dp(Ar(1, 1), 0.0_dp, 1.0e-14_dp, 'mask dp zero diag', failures)
    call expect_close_scalar_dp(Ar(2, 2), 0.25_dp, 1.0e-14_dp, 'mask dp invert diag', failures)

    As = reshape([0.0_sp, 0.0_sp, 0.0_sp, 9.0_sp], [2, 2])
    call invert_matrix_with_mask_sp(As)
    call expect_close_scalar_sp(As(1, 1), 0.0_sp, 1.0e-6_sp, 'mask sp zero diag', failures)
    call expect_close_scalar_sp(As(2, 2), 1.0_sp / 9.0_sp, 2.0e-6_sp, 'mask sp invert diag', failures)

    Ac = reshape([cmplx(0.0_dp, 0.0_dp, dpc), cmplx(0.0_dp, 0.0_dp, dpc), &
      & cmplx(0.0_dp, 0.0_dp, dpc), cmplx(4.0_dp, 0.0_dp, dpc)], [2, 2])
    call invert_matrix_with_mask_dpc(Ac)
    call expect_close_scalar_dp(real(Ac(1, 1), dp), 0.0_dp, 1.0e-14_dp, 'mask dpc zero diag', failures)
    call expect_close_scalar_dp(real(Ac(2, 2), dp), 0.25_dp, 1.0e-14_dp, 'mask dpc invert diag', failures)
  end subroutine test_invert_matrix_with_mask_real_and_complex

  subroutine test_get_eigenvalues_and_vectors(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A, V
    real(dp), dimension(2) :: W
    real(dp), dimension(3, 3) :: B
    real(dp), dimension(3) :: WB

    A = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
    call get_eigenvalues(A, W)
    call expect_close_scalar_dp(W(1), 3.0_dp, 1.0e-12_dp, 'get_eigenvalues n2 w1', failures)
    call expect_close_scalar_dp(W(2), 1.0_dp, 1.0e-12_dp, 'get_eigenvalues n2 w2', failures)

    call get_eigen_decomposition(A, W, V)
    call expect_eigen_residual(A, W, V, 1.0e-10_dp, 'get_eigen_decomposition residual', failures)
    call expect_orthonormal_columns(V, 1.0e-10_dp, 'get_eigen_decomposition orthogonality', failures)

    B = 0.0_dp
    B(1, 1) = 1.0_dp
    B(2, 2) = 3.0_dp
    B(3, 3) = 2.0_dp
    call get_eigenvalues(B, WB)
    call expect_close_scalar_dp(WB(1), 1.0_dp, 1.0e-12_dp, 'get_eigenvalues n3 w1', failures)
    call expect_close_scalar_dp(WB(2), 2.0_dp, 1.0e-12_dp, 'get_eigenvalues n3 w2', failures)
    call expect_close_scalar_dp(WB(3), 3.0_dp, 1.0e-12_dp, 'get_eigenvalues n3 w3', failures)
  end subroutine test_get_eigenvalues_and_vectors

  subroutine test_eigen_decomp_wrapper(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A, V
    real(dp), dimension(2) :: W
    integer(i4b) :: status

    A = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
    call eigen_decomp(A, W, V, status)
    if (status /= 0_i4b) then
      failures = failures + 1
      write(*, '(A,I0)') 'FAIL: eigen_decomp status=', status
    end if
    call expect_eigen_residual(A, W, V, 1.0e-10_dp, 'eigen_decomp residual', failures)
    call expect_orthonormal_columns(V, 1.0e-10_dp, 'eigen_decomp orthogonality', failures)
  end subroutine test_eigen_decomp_wrapper

  subroutine test_compute_hermitian_root(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A, Aorig

    A = reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
    Aorig = A
    call compute_hermitian_root(A, 0.5_dp)
    call expect_close_mat_dp(matmul(A, A), Aorig, 1.0e-9_dp, 'compute_hermitian_root sqrt', failures)
  end subroutine test_compute_hermitian_root

  subroutine test_compute_hermitian_root_non_spd(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A

    A = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2])
    call compute_hermitian_root(A, 0.5_dp)
    call expect_close_scalar_dp(A(1, 1), -1.0d30, 0.0_dp, 'compute_hermitian_root non-spd sentinel', failures)
  end subroutine test_compute_hermitian_root_non_spd

  subroutine test_solve_systems(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: Ar
    real(dp), dimension(2) :: br, xr
    complex(dpc), dimension(2, 2) :: Ac
    complex(dpc), dimension(2) :: bc, xc

    Ar = reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
    br = [9.0_dp, 8.0_dp]
    call solve_system_real(Ar, xr, br)
    call expect_close_vec_dp(xr, [1.7272727272727273_dp, 2.0909090909090908_dp], 1.0e-12_dp, 'solve_system_real', failures)

    Ac = reshape([cmplx(2.0_dp, 0.0_dp, dpc), cmplx(0.0_dp, 0.0_dp, dpc), &
      & cmplx(0.0_dp, 0.0_dp, dpc), cmplx(3.0_dp, 0.0_dp, dpc)], [2, 2])
    bc = [cmplx(4.0_dp, 0.0_dp, dpc), cmplx(9.0_dp, 0.0_dp, dpc)]
    call solve_system(Ac, xc, bc)
    call expect_close_scalar_dp(real(xc(1), dp), 2.0_dp, 1.0e-12_dp, 'solve_system complex x1', failures)
    call expect_close_scalar_dp(real(xc(2), dp), 3.0_dp, 1.0e-12_dp, 'solve_system complex x2', failures)
  end subroutine test_solve_systems

  subroutine test_invert_singular_matrix_threshold(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A

    A = reshape([4.0_dp, 0.0_dp, 0.0_dp, 1.0e-8_dp], [2, 2])
    call invert_singular_matrix(A, 0.1_dp)
    call expect_close_scalar_dp(A(1, 1), 0.25_dp, 1.0e-12_dp, 'invert_singular_matrix kept mode', failures)
    call expect_close_scalar_dp(A(2, 2), 0.0_dp, 1.0e-12_dp, 'invert_singular_matrix cutoff mode', failures)
  end subroutine test_invert_singular_matrix_threshold

  subroutine test_probability_conversions(failures)
    integer(i4b), intent(inout) :: failures
    real(sp) :: sigma_sp, frac_sp, sigma_sp_back
    real(dp) :: sigma_dp_back
    real(dp) :: fract_dp

    sigma_sp = 2.0_sp
    call convert_sigma2frac(frac_sp, sigma_sp)
    fract_dp = corr_erf(real(sigma_sp, dp) / sqrt(2.0_dp))
    call convert_fract2sigma_sp(sigma_sp_back, real(fract_dp, sp))
    call convert_fract2sigma_dp(sigma_dp_back, fract_dp)

    call expect_close_scalar_sp(sigma_sp_back, sigma_sp, 2.0e-4_sp, 'fract2sigma_sp roundtrip', failures)
    call expect_close_scalar_dp(sigma_dp_back, real(sigma_sp, dp), 2.0e-4_dp, 'fract2sigma_dp roundtrip', failures)
    call expect_close_scalar_sp(frac_sp, 0.97724986_sp, 5.0e-5_sp, 'convert_sigma2frac sigma2', failures)
  end subroutine test_probability_conversions

  subroutine test_cholesky_decompose_single(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A, L, recon
    integer(i4b) :: ierr

    A = reshape([4.0_dp, 2.0_dp, 2.0_dp, 3.0_dp], [2, 2])
    L = A
    call cholesky_decompose_single(L, ierr)
    call expect_equal_i4(ierr, 0_i4b, 'cholesky_decompose_single ierr', failures)
    call expect_close_scalar_dp(L(1, 2), 0.0_dp, 1.0e-12_dp, 'cholesky_decompose_single upper zero', failures)
    recon = matmul(L, transpose(L))
    call expect_close_mat_dp(recon, A, 1.0e-12_dp, 'cholesky_decompose_single recon', failures)
  end subroutine test_cholesky_decompose_single

  subroutine test_matmul_wrappers(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(2, 2) :: A, B, C1, C2

    A = reshape([2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
    B = reshape([1.0_dp, 4.0_dp, 2.0_dp, 5.0_dp], [2, 2])

    call matmul_symm(A, B, C1)
    C2 = matmul(A, B)
    call expect_close_mat_dp(C1, C2, 1.0e-12_dp, 'matmul_symm', failures)

    call matmul_gen(A, B, C1)
    call expect_close_mat_dp(C1, C2, 1.0e-12_dp, 'matmul_gen', failures)
  end subroutine test_matmul_wrappers

  subroutine test_mov3(failures)
    integer(i4b), intent(inout) :: failures
    real(dp) :: a1, b1, c1, a2, b2, c2

    a1 = 1.0_dp; b1 = 2.0_dp; c1 = 3.0_dp
    a2 = 7.0_dp; b2 = 8.0_dp; c2 = 9.0_dp
    call MOV3(a1, b1, c1, a2, b2, c2)

    call expect_close_scalar_dp(a1, 7.0_dp, 1.0e-12_dp, 'MOV3 a1', failures)
    call expect_close_scalar_dp(b1, 8.0_dp, 1.0e-12_dp, 'MOV3 b1', failures)
    call expect_close_scalar_dp(c1, 9.0_dp, 1.0e-12_dp, 'MOV3 c1', failures)
    call expect_close_scalar_dp(a2, 1.0_dp, 1.0e-12_dp, 'MOV3 a2', failures)
    call expect_close_scalar_dp(b2, 2.0_dp, 1.0e-12_dp, 'MOV3 b2', failures)
    call expect_close_scalar_dp(c2, 3.0_dp, 1.0e-12_dp, 'MOV3 c2', failures)
  end subroutine test_mov3

  subroutine test_comp_normalised_plm_smoke(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(0:0) :: p0
    real(dp), dimension(0:3) :: p3

    call comp_normalised_Plm(0_i4b, 0_i4b, 0.7_dp, p0)
    call expect_close_scalar_dp(p0(0), 1.0_dp / sqrt(4.0_dp * pi), 1.0e-12_dp, 'comp_normalised_Plm l0m0', failures)

    call comp_normalised_Plm(3_i4b, 1_i4b, 0.7_dp, p3)
    if (.not. all(abs(p3) < huge(1.0_dp))) then
      failures = failures + 1
      write(*, '(A)') 'FAIL: comp_normalised_Plm finite output'
    end if
  end subroutine test_comp_normalised_plm_smoke

  subroutine expect_eigen_residual(A, W, V, tol, label, failures)
    real(dp), dimension(:, :), intent(in) :: A, V
    real(dp), dimension(:), intent(in) :: W
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures
    integer(i4b) :: i
    real(dp), dimension(size(W)) :: res

    do i = 1, size(W)
      res = matmul(A, V(:, i)) - W(i) * V(:, i)
      if (sqrt(sum(res**2)) > tol) then
        failures = failures + 1
        write(*, '(A,A,A,I0)') 'FAIL: ', trim(label), ' vec=', i
      end if
    end do
  end subroutine expect_eigen_residual

  subroutine expect_orthonormal_columns(V, tol, label, failures)
    real(dp), dimension(:, :), intent(in) :: V
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(size(V, 2), size(V, 2)) :: gram
    real(dp), dimension(size(V, 2), size(V, 2)) :: ident

    gram = matmul(transpose(V), V)
    ident = 0.0_dp
    ident = identity(size(V, 2))

    if (any(abs(gram - ident) > tol)) then
      failures = failures + 1
      write(*, '(A,A)') 'FAIL: ', trim(label)
    end if
  end subroutine expect_orthonormal_columns

  function identity(n) result(E)
    integer(i4b), intent(in) :: n
    real(dp), dimension(n, n) :: E
    integer(i4b) :: i

    E = 0.0_dp
    do i = 1, n
      E(i, i) = 1.0_dp
    end do
  end function identity

  function eye2_sp() result(E)
    real(sp), dimension(2, 2) :: E
    E = 0.0_sp
    E(1, 1) = 1.0_sp
    E(2, 2) = 1.0_sp
  end function eye2_sp

  function eye2_cplx() result(E)
    complex(dpc), dimension(2, 2) :: E
    E = cmplx(0.0_dp, 0.0_dp, dpc)
    E(1, 1) = cmplx(1.0_dp, 0.0_dp, dpc)
    E(2, 2) = cmplx(1.0_dp, 0.0_dp, dpc)
  end function eye2_cplx

  subroutine expect_equal_i4(actual, expected, label, failures)
    integer(i4b), intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (actual /= expected) then
      failures = failures + 1
      write(*, '(A,A,A,I0,A,I0)') 'FAIL: ', trim(label), ' expected=', expected, ' actual=', actual
    end if
  end subroutine expect_equal_i4

  subroutine expect_close_scalar_dp(actual, expected, tol, label, failures)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (abs(actual - expected) > tol) then
      failures = failures + 1
      write(*, '(A,A,A,ES12.5,A,ES12.5)') 'FAIL: ', trim(label), ' expected=', expected, ' actual=', actual
    end if
  end subroutine expect_close_scalar_dp

  subroutine expect_close_scalar_sp(actual, expected, tol, label, failures)
    real(sp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (abs(actual - expected) > tol) then
      failures = failures + 1
      write(*, '(A,A,A,ES12.5,A,ES12.5)') 'FAIL: ', trim(label), ' expected=', real(expected, dp), ' actual=', real(actual, dp)
    end if
  end subroutine expect_close_scalar_sp

  subroutine expect_close_vec_dp(actual, expected, tol, label, failures)
    real(dp), dimension(:), intent(in) :: actual, expected
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (any(abs(actual - expected) > tol)) then
      failures = failures + 1
      write(*, '(A,A)') 'FAIL: ', trim(label)
    end if
  end subroutine expect_close_vec_dp

  subroutine expect_close_mat_dp(actual, expected, tol, label, failures)
    real(dp), dimension(:, :), intent(in) :: actual, expected
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (any(abs(actual - expected) > tol)) then
      failures = failures + 1
      write(*, '(A,A)') 'FAIL: ', trim(label)
    end if
  end subroutine expect_close_mat_dp

  subroutine expect_close_mat_sp(actual, expected, tol, label, failures)
    real(sp), dimension(:, :), intent(in) :: actual, expected
    real(sp), intent(in) :: tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (any(abs(actual - expected) > tol)) then
      failures = failures + 1
      write(*, '(A,A)') 'FAIL: ', trim(label)
    end if
  end subroutine expect_close_mat_sp

  subroutine expect_close_mat_cplx(actual, expected, tol, label, failures)
    complex(dpc), dimension(:, :), intent(in) :: actual, expected
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (any(abs(actual - expected) > tol)) then
      failures = failures + 1
      write(*, '(A,A)') 'FAIL: ', trim(label)
    end if
  end subroutine expect_close_mat_cplx

end program test_math_tools_full
