# dang_mbb_comp_mod

Defines a modified blackbody component with two spectral parameters (`BETA` and `T`).
Its evaluator computes the MBB spectrum with Planck-factor ratios for single-frequency and bandpass-integrated cases, using a 2D mixing-matrix spline backend.

In practice, this is the thermal dust workhorse component: it models dust emission amplitude, spectral index, and temperature.
The 2D spline acceleration keeps repeated bandpass evaluations fast during iterative fitting.

Its SED is the usual modified blackbody form,

$$
S(\nu;\beta,T) \propto \left(\frac{\nu}{\nu_{\mathrm{ref}}}\right)^\beta
\frac{B_\nu(T)}{B_{\nu_{\mathrm{ref}}}(T)},
$$

so map prediction is \(m_{\nu,p}=a_p\,S(\nu;\beta_p,T_p)\) (or a bandpass integral of that SED).
