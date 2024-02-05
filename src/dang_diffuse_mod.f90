module dang_diffuse
  use dang_component_mod
  use dang_data_mod
  use dang_bp_mod

  type :: dang_diffuse
     class(dang_comps), pointer :: comp => null()
     class(dang_bp),    pointer :: bp   => null()
     type(spline_type) :: s

  end type dang_diffuse

  integer(i4b) :: n_spline = 1000
  
contains

  procedure :: initCompSpline
  procedure :: update
  
  ! stores the component and bandpass info
  function constructor(comp, bp, pol)
    implicit none
    class(dang_comps),     intent(in), target :: comp
    class(dang_bp),        intent(in), target :: bp
    integer(i4b),          intent(in)         :: pol

    allocate(constructor)
    constructor%comp => comp
    constructor%bp   => bp

    call constructor%initCompSpline(pol=pol)
    
  end function constructor
  
  subroutine initCompSpline()
    implicit none
    class(dang_comps)             :: self

    
    if (self%comp%nindices == 0) then
       ! No need to spline
    else if (self%comp%nindices == 1) then
       



end module dang_diffuse
