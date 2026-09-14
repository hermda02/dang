# dang_template_comp_mod

Defines a template-driven component with no sampled spectral indices and per-band fitted template amplitudes.
It loads/normalizes a template map and evaluates template contributions only for bands flagged as correlated.

In practice, this module is used for externally supplied spatial templates where only channel scaling needs to be fit.
It is a convenient way to inject known morphology into the component-separation model.

The model is template-linear,

$$
m_{\nu,p} = a_\nu\,t_p,
$$

with \(t_p\) fixed from the input map and \(a_\nu\) solved only for bands marked as correlated.
