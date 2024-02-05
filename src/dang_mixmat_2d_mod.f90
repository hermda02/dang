module dang_mixmat_2d_mod
  use dang_mixmat_mod
  use dang_util_mod
  use dang_component_mod
  use dang_bp_mod
  use spline_2d_mod

  private
  public dang_mixmat_2d

  type, extends (dang_mixmat) :: dang_mixmat_2d
     real(dp), allocatable, dimension(:)              :: x, y
     real(dp), allocatable, dimension(:,:,:,:)        :: coeff
     class(dang_comp), pointer :: comp => null()
   contains
     procedure :: init_spline
     procedure :: eval    => evalIntegratedSignal
     procedure :: eval_dI => evalIntegratedDSignal 
  end type dang_mixmat_2d

  interface dang_mixmat_2d
     procedure constructor_mixmat_2d
  end interface dang_mixmat_2d

  integer(i4b) :: n=100

contains

  function constructor_mixmat_2d(comp, bp, pol)
    implicit none
    class(dang_comp),    intent(in), target :: comp
    class(bandinfo),      intent(in), target :: bp
    integer(i4b),         intent(in)         :: pol
    class(dang_mixmat_2d),          pointer :: constructor_mixmat_2d

    allocate(constructor_mixmat_2d)
    constructor_mixmat_2d%comp  => comp
    constructor_mixmat_2d%bp    => bp

    call constructor_mixmat_2d%init_spline(pol)

  end function constructor_mixmat_2d
  
  subroutine init_spline(self, pol)
    implicit none
    class(dang_mixmat_2d)                 :: self
    integer(i4b),   intent(in), optional  :: pol
    
    real(dp), allocatable, dimension(:)   :: theta, s
    real(dp), allocatable, dimension(:,:) :: f

    integer(i4b)                          :: i, j, k
    integer(i4b)                          :: m

    m = self%bp%n
    
    allocate(f(n,n))
    allocate(s(m))
    allocate(self%x(n))
    allocate(self%y(n))
    allocate(self%coeff(4,4,n,n))
    do i = 1, n
       self%x(i) = self%comp%uni_prior(1,1) + (self%comp%uni_prior(1,2)-self%comp%uni_prior(1,1))/(n-1) * (i-1)
       self%y(i) = self%comp%uni_prior(2,1) + (self%comp%uni_prior(2,2)-self%comp%uni_prior(2,1))/(n-1) * (i-1)
    end do
    do i = 1, n
       do j = 1, n
          do k = 1, m
             s(k) = self%comp%S(nu=self%comp%nu_ref,pol=pol,theta=[self%x(i),self%y(j)])
          end do
          f(i,j) = self%bp%integrate(s)
       end do
    end do
    call splie2_full_precomp(self%x, self%y, f, self%coeff)

  end subroutine init_spline

  function evalIntegratedSignal(self, theta)
    implicit none
    class(dang_mixmat_2d),             intent(in) :: self
    real(dp),           dimension(1:), intent(in) :: theta
    real(dp)                                      :: evalIntegratedSignal
    evalIntegratedSignal = splin2_full_precomp(self%x, self%y, self%coeff, theta(1), theta(2))
  end function evalIntegratedSignal

  function evalIntegratedDSignal(self, theta, par)
    implicit none
    class(dang_mixmat_2d),             intent(in) :: self
    real(dp),           dimension(1:), intent(in) :: theta
    integer(i4b),                      intent(in) :: par
    real(dp)                                      :: evalIntegratedDSignal

    real(dp) :: p(2), delta = 1.d-10, f1, f2
    
    p      = theta
    p(par) = p(par) + delta
    f1     = self%eval(theta)
    f2     = self%eval(p)
    evalIntegratedDSignal = (f2-f1)/delta
  end function evalIntegratedDSignal

end module
