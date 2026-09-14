# dang_freefree_comp_mod

Implements a free-free emission component with one spectral parameter (`T_e`) and standard prior/sampling map handling.
Its evaluator computes normalized free-free spectral scaling at a frequency or through bandpass integration.

In practice, this module models thermal bremsstrahlung foreground emission, letting the pipeline separate that contribution from CMB and other foregrounds.
Its temperature-like parameter can be sampled while amplitudes are solved in the linear step.

The per-band model is

$$
m_{\nu,p}=a_{p}\,S_{\mathrm{ff}}(\nu;T_e),
$$

with `S_ff` computed from the free-free spectral law (including temperature-dependent terms) and normalized to the component reference frequency.
