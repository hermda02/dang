# dang_cg_mod

Implements conjugate-gradient machinery for solving or sampling grouped component amplitudes.
It builds right-hand-side vectors, applies the implicit normal-equation operator (`compute_Ax`), runs CG iterations, and maps solved values back into component amplitude/template structures.

In practice, this is the module that solves the big linear system behind each amplitude update step without ever forming a huge dense matrix explicitly.
Given current component spectra, noise weights, and masks, it repeatedly computes matrix-vector products and converges to the best-fit (or sampled) component amplitudes used in map reconstruction.

The solved system is the weighted normal equation for one CG group and one polarization selection:

$$
\left(\sum_{\nu} T_{\nu}^{\mathsf T} N_{\nu}^{-1} T_{\nu}\right) x
=
\sum_{\nu} T_{\nu}^{\mathsf T} N_{\nu}^{-1} d_{\nu}^{\ast},
$$

where `compute_rhs` forms the right-hand side after subtracting non-group components from data (`d_\nu^*`), and `compute_Ax` applies the left-hand operator matrix-free as:

$$
x \xrightarrow{\;T_\nu\;} y_\nu \xrightarrow{\;N_\nu^{-1}\;} z_\nu \xrightarrow{\;T_\nu^{\mathsf T}\;} w_\nu,
\qquad
Ax = \sum_\nu w_\nu.
$$

For per-pixel amplitudes, entries are weighted by inverse noise variance, e.g.

$$
b_{c,p} = \sum_{\nu}\frac{d_{\nu,p}^{\ast}\,S_{c,\nu,p}}{\sigma_{\nu,p}^2},
\qquad
(Ax)_{c,p} = \sum_{\nu}\frac{S_{c,\nu,p}}{\sigma_{\nu,p}^2}\sum_{c'} S_{c',\nu,p}x_{c',p},
$$

with template/monopole/`hi_fit` terms collapsed over pixels into band-level coefficients rather than one unknown per pixel.

In sampling mode, the RHS gets a fluctuation term:

$$
b \leftarrow b + \sum_{\nu}T_{\nu}^{\mathsf T}N_{\nu}^{-1/2}\eta,
\qquad \eta \sim \mathcal N(0,I),
$$

so the CG solve returns a constrained Gaussian sample instead of only the MAP/ML amplitude solution.
