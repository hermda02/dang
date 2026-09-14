# dang_hi_fit_comp_mod

Implements a specialized HI-template fitting component with one sampled temperature-like index and per-band template amplitudes.
Its signal model combines template amplitude with a modified-blackbody-like frequency term, with support for priors and bandpass-integrated evaluation.

In practice, this module is used when neutral-hydrogen-correlated emission is fit as a dedicated component instead of being absorbed into broader dust terms.
It enables band-specific template coupling while still keeping a physically parameterized frequency scaling.

Its fitted signal is effectively

$$
m_{\nu,p}=A_{\nu}\,T_{\mathrm{HI},p}\,f_{\nu}(T),
$$

where \(A_\nu\) are solved template amplitudes for correlated bands and \(f_\nu(T)\) is the temperature-dependent spectral term.
