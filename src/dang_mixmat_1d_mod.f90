module dang_mixmat_1d_mod
  use dang_mixmat_mod
  use dang_util_mod
  use dang_component_mod
  use dang_bp_mod
  use spline_1d_mod

  private
  public dang_mixmat_1d

  type, extends (dang_mixmat) :: dang_mixmat_1d
     type(spline_type) :: s
     class(dang_comp), pointer :: comp => null()
   contains
     procedure :: init_spline
     procedure :: eval    => evalIntegratedSignal
     procedure :: eval_dI => evalIntegratedDSignal 
  end type dang_mixmat_1d

  interface dang_mixmat_1d
     procedure constructor_mixmat
  end interface dang_mixmat_1d

  integer(i4b) :: n = 1000

contains

  function constructor_mixmat(comp, bp, pol)
    implicit none
    class(dang_comp),     intent(in), target :: comp
    class(bandinfo),      intent(in), target :: bp
    integer(i4b),         intent(in)         :: pol
    class(dang_mixmat_1d),           pointer :: constructor_mixmat

    allocate(constructor_mixmat)
    constructor_mixmat%comp  => comp
    constructor_mixmat%bp    => bp

    call constructor_mixmat%init_spline(pol)

  end function constructor_mixmat
  
  subroutine init_spline(self, pol)
    implicit none
    class(dang_mixmat_1d)                     :: self
    integer(i4b),   intent(in), optional  :: pol
    
    real(dp), allocatable, dimension(:)   :: theta
    real(dp), allocatable, dimension(:)   :: f, s

    integer(i4b)                          :: i, j, k
    integer(i4b)                          :: m

    m = self%bp%n
    
    allocate(theta(n))
    allocate(s(m))
    allocate(f(n))
    do i = 1, n
       theta(i) = self%comp%uni_prior(1,1) + (self%comp%uni_prior(1,2)-self%comp%uni_prior(1,1))/(n-1) * (i-1)
       do j = 1, m
          s(j)  = self%comp%S(nu=self%bp%nu0(j),pol=pol,theta=theta)
       end do
       f(i) = self%bp%integrate(s)
    end do
    call spline_simple(self%s, theta, f, regular=.true.)

  end subroutine init_spline

  function evalIntegratedSignal(self, theta)
    implicit none
    class(dang_mixmat_1d),                 intent(in) :: self
    real(dp),           dimension(1:), intent(in) :: theta
    real(dp)                                      :: evalIntegratedSignal
    evalIntegratedSignal = splint_simple(self%s, theta(1))
  end function evalIntegratedSignal

  function evalIntegratedDSignal(self, theta, par)
    implicit none
    class(dang_mixmat_1d),                 intent(in) :: self
    real(dp),           dimension(1:), intent(in) :: theta
    integer(i4b),                      intent(in) :: par
    real(dp)                                      :: evalIntegratedDSignal

    real(dp) :: p(2), delta = 1.d-10, f1, f2
    
    evalIntegratedDSignal = splint_deriv(self%s%x, self%s%y, self%s%y2, theta(1))
    
  end function evalIntegratedDSignal

end module
