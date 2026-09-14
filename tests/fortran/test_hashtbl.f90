program test_hashtbl
  use hashtbl, only: hash_tbl_sll, sum_string
  implicit none

  integer :: failures
  type(hash_tbl_sll) :: tbl
  character(len=:), allocatable :: val

  failures = 0

  call test_sum_string_equation(failures)
  call test_put_get_roundtrip(failures)
  call test_collision_handling(failures)
  call test_missing_key_behavior(failures)

  if (failures /= 0) then
    write(*, '(A,I0)') "test_hashtbl failures: ", failures
    stop 1
  end if

  write(*, '(A)') "test_hashtbl passed"

contains

  subroutine test_sum_string_equation(failures)
    integer, intent(inout) :: failures
    integer :: expected

    expected = ichar('a') + ichar('B') + ichar('3')
    if (sum_string('aB3') /= expected) then
      failures = failures + 1
      write(*, '(A)') "FAIL: sum_string should equal sum of character codes"
    end if
  end subroutine test_sum_string_equation

  subroutine test_put_get_roundtrip(failures)
    integer, intent(inout) :: failures
    type(hash_tbl_sll) :: t
    character(len=:), allocatable :: out

    call t%init(32)
    call t%put('alpha', '42')
    call t%get('alpha', out)

    if (.not. allocated(out) .or. out /= '42') then
      failures = failures + 1
      write(*, '(A)') "FAIL: hashtable put/get roundtrip"
    end if

    call t%free()
  end subroutine test_put_get_roundtrip

  subroutine test_collision_handling(failures)
    integer, intent(inout) :: failures
    type(hash_tbl_sll) :: t
    character(len=:), allocatable :: out

    call t%init(1)
    call t%put('k1', 'v1')
    call t%put('k2', 'v2')
    call t%put('k3', 'v3')

    call t%get('k1', out)
    if (.not. allocated(out) .or. out /= 'v1') then
      failures = failures + 1
      write(*, '(A)') "FAIL: collision retrieval k1"
    end if

    if (allocated(out)) deallocate(out)
    call t%get('k2', out)
    if (.not. allocated(out) .or. out /= 'v2') then
      failures = failures + 1
      write(*, '(A)') "FAIL: collision retrieval k2"
    end if

    call t%free()
  end subroutine test_collision_handling

  subroutine test_missing_key_behavior(failures)
    integer, intent(inout) :: failures
    type(hash_tbl_sll) :: t
    character(len=:), allocatable :: out

    call t%init(8)
    call t%put('known', 'value')
    call t%get('unknown', out)
    if (allocated(out)) then
      failures = failures + 1
      write(*, '(A)') "FAIL: missing key should leave output unallocated"
      deallocate(out)
    end if
    call t%free()
  end subroutine test_missing_key_behavior

end program test_hashtbl
