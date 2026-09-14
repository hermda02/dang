# dang_mixmat_mod

Defines the abstract mixing-matrix interface used by component models.
It specifies deferred APIs for spline initialization, integrated signal evaluation, and derivatives, and includes pointer-wrapper types used throughout component code.

In practice, this gives component code a uniform way to ask for bandpass-integrated SED values regardless of whether the backend is 1D or 2D.
That abstraction keeps model logic clean while allowing specialized interpolation implementations underneath.

Conceptually, implementations expose a common map

$$
(b,\theta)\mapsto M_b(\theta)=\int R_b(\nu)S(\nu;\theta)\,d\nu,
$$

plus derivative access \(\partial M_b/\partial\theta_i\) for proposal/optimization logic.
