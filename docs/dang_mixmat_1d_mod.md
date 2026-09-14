# dang_mixmat_1d_mod

Provides a concrete 1D mixing-matrix implementation that precomputes spline tables over one spectral parameter.
It supports fast bandpass-integrated interpolation of signal values and corresponding derivatives.

In practice, one-parameter components use this module to avoid expensive on-the-fly bandpass integrals at every sampler step.
It trades small precomputation cost for much faster repeated evaluations.

It precomputes and splines

$$
M_b(\theta)=\int R_b(\nu)\,S(\nu;\theta)\,d\nu,
$$

so runtime calls are mostly interpolation \(M_b(\theta)\) and \(\partial M_b/\partial\theta\).
