# dang_lognormal_comp_mod

Implements a two-parameter lognormal-like component with parameters such as `nu_p` and `w_ame`.
The evaluator computes the lognormal spectral shape (including normalization factors) at direct frequencies or via bandpass integration.

In practice, this module captures peaked foreground spectra (for example AME-like behavior) that cannot be represented well by a simple power law.
It lets the model fit both peak location and width while solving amplitudes jointly with other components.

Its spectrum is lognormal-like,

$$
S(\nu) \propto \left(\frac{\nu_{\mathrm{ref}}}{\nu}\right)^2
\exp\!\left[-\frac{(\ln(\nu/\nu_p))^2}{2w_{\mathrm{ame}}^2}\right],
$$

then used as \(m_{\nu,p}=a_p\,S(\nu)\) (or bandpass-integrated equivalent).
