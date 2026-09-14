# dang_powlaw_comp_mod

Implements a power-law foreground component with one spectral index (`beta`).
The constructor loads amplitudes/index priors and maps, and the evaluator computes normalized power-law scaling at direct frequency or integrated over a bandpass.

In practice, this is the standard module for simple scale-free foregrounds (such as synchrotron-like behavior in some regimes).
It is lightweight and efficient when one index parameter is enough to describe the spectrum.

The component SED is

$$
S(\nu;\beta)=\left(\frac{\nu}{\nu_{\mathrm{ref}}}\right)^\beta,
$$

used in predictions as \(m_{\nu,p}=a_p\,S(\nu;\beta_p)\) (or bandpass-integrated form).
