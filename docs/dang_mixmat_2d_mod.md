# dang_mixmat_2d_mod

Provides a concrete 2D mixing-matrix implementation with bicubic spline coefficients over two spectral parameters.
It enables fast interpolated bandpass-integrated evaluation and derivative estimation for multi-parameter SED components.

In practice, this module is critical for components like MBB where two spectral parameters vary during sampling.
It keeps the runtime tractable by replacing repeated 2D integrations with spline lookups and finite-difference derivatives.

It tabulates a 2D surface,

$$
M_b(\theta_1,\theta_2)=\int R_b(\nu)\,S(\nu;\theta_1,\theta_2)\,d\nu,
$$

then evaluates/interpolates \(M_b\) and numerical derivatives during sampling.
