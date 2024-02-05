module dang_mixmat_mod
  use dang_util_mod
  use dang_bp_mod
  
  private
  public dang_mixmat, mixmat_ptr

  type, abstract :: dang_mixmat
     class(bandinfo), pointer :: bp
   contains
     procedure(init_spline),              deferred :: init_spline
     procedure(evalIntegratedSignal),     deferred :: eval
     procedure(evalIntegratedDSignal),    deferred :: eval_dI
  end type dang_mixmat

  abstract interface
     subroutine init_spline(self, pol)
       import :: dang_mixmat, i4b
       class(dang_mixmat)                     :: self
       integer(i4b),   intent(in), optional  :: pol
     end subroutine init_spline
     
     function evalIntegratedSignal(self, theta)
       import :: dang_mixmat, dp
       class(dang_mixmat),         intent(in) :: self
       real(dp),    dimension(1:), intent(in) :: theta
       real(dp)                               :: evalIntegratedSignal
     end function evalIntegratedSignal

     function evalIntegratedDSignal(self, theta, par)
       import :: dang_mixmat, dp, i4b
       class(dang_mixmat),      intent(in) :: self
       real(dp), dimension(1:), intent(in) :: theta
       integer(i4b),            intent(in) :: par
       real(dp)                            :: evalIntegratedDSignal
     end function evalIntegratedDSignal
  end interface
     
  type mixmat_ptr
     class(dang_mixmat), pointer :: p
  end type mixmat_ptr

end module dang_mixmat_mod
