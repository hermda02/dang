module init_mod
    use healpix_types
    implicit none

    ! Constants
    real(dp)           :: k_B     = 1.3806503d-23
    real(dp)           :: h       = 1.0545726691251021d-34 * 2.d0 * pi
    real(dp)           :: c       = 2.99792458d8
    real(dp)           :: T_CMB   = 2.7255d0

    integer(i4b)       :: nbands

contains 

    function getlun()
        implicit none
        integer(i4b) :: getlun
        logical(lgt) :: exists, isopen
        getlun = 9
        do
        getlun = getlun+1
        inquire(unit=getlun,exist=exists)
        if(exists) then
            inquire(unit=getlun,opened=isopen)
            if(.not. isopen) return
        end if
        end do
    end function getlun

end module init_mod