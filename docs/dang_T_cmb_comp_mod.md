# dang_T_cmb_comp_mod

Implements a one-parameter CMB-temperature component with a global temperature-like spectral parameter.
Its evaluator computes a Planck-spectrum-based scaling term (or bandpass equivalent) and returns values in microkelvin units.

In practice, this module allows fitting a temperature-dependent CMB-like term beyond fixed anisotropy scaling.
It is useful when analyses need explicit control of a global CMB temperature parameter in the spectral model.

Its spectral factor is Planck-based,

$$
S(\nu;T) \propto \frac{B_\nu(T)}{B'_{\nu,\mathrm{RJ}}},
$$

so predictions take the form \(m_{\nu,p}=a_p\,S(\nu;T)\) (or bandpass-integrated equivalent).
