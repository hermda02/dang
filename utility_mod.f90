module utility_mod
    use healpix_types
    implicit none
    include 'mpif.h'

    integer(i4b) :: ierr, rank, numprocs

contains 

  ! Small utility for converting an integer to a string                                              
  subroutine int2string(integer, string)
    implicit none

    integer(i4b),     intent(in)  :: integer
    character(len=*), intent(out) :: string

    integer(i4b)               :: temp_int, i, k

    temp_int = integer
    do i = 1, len(string)
       k = temp_int / 10**(len(string)-i)
       write(string(i:i),'(I1)') k
       temp_int = temp_int - k * 10**(len(string)-i)
    end do

  end subroutine int2string

    subroutine tolower(strIn)! result(strOut)
        ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
        ! Original author: Clive Page - switched to tolower
        
             implicit none
        
             character(len=*), intent(inout) :: strIn
             character(len=len(strIn)) :: strOut
             integer :: i,j
        
             do i = 1, len(strIn)
                  j = iachar(strIn(i:i))
                  if (j>= iachar("A") .and. j<=iachar("Z") ) then
                       strOut(i:i) = achar(iachar(strIn(i:i))+32)
                  else
                       strOut(i:i) = strIn(i:i)
                  end if
             end do

            strIn = strOut
        
    end subroutine tolower

end module utility_mod
