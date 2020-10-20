module dang_component_mod
  use healpix_types
  use utility_mod
  use dang_param_mod
  implicit none

  type, public                              :: component

     character(len=64), allocatable, dimension(:) :: joint

     real(dp), allocatable, dimension(:,:)  :: beta_s, beta_d, T_d

  end type component

contains 

    function planck(fre,T)
        implicit none
        real(dp), intent(in)  :: fre
        real(dp), intent(in)  :: T
        real(dp)              :: planck
        ! Output in units of [W sr^-1 m^-2 Hz^-1]
        planck  = ((2.d0*h*fre**3.d0)/(c**2.d0))*(1.d0/(exp((h*fre)/(k_B*T))-1))
    end function planck

    function compute_spectrum(self, comp, ind, freq, pix, mapn, param)
        implicit none
        class(params)                  :: self
        type(component)                :: comp
        real(dp),           intent(in) :: freq
        integer(i4b),       intent(in) :: ind
        integer(i4b),       intent(in) :: pix
        integer(i4b),       intent(in) :: mapn
        real(dp), optional             :: param
        real(dp)                       :: y, rj_cmb ! Conversion factor from K_{RJ} -> K_{CMB}
        real(dp)                       :: z, compute_spectrum

        y = h*(freq*1.0d9) / (k_B*T_CMB)
        rj_cmb = ((exp(y)-1)**2.d0)/(y**2.d0*exp(y))

        if (trim(self%fg_label(ind)) == 'synch') then
           if (present(param)) then
              compute_spectrum = (freq/self%fg_nu_ref(ind))**param
           else 
              compute_spectrum = (freq/self%fg_nu_ref(ind))**comp%beta_s(pix,mapn)
           end if
        else if (trim(self%fg_label(ind)) == 'mbb') then
           z = h / (k_B*comp%T_d(pix,mapn))
           compute_spectrum = (exp(z*self%fg_nu_ref(ind)*1d9)-1.d0) / &
                (exp(z*freq*1d9)-1.d0) * (freq/self%fg_nu_ref(ind))**(comp%beta_d(pix,mapn)+1.d0)!*rj_cmb
        end if

    end function compute_spectrum



end module dang_component_mod
