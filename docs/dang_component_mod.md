# dang_component_mod

Defines the abstract base component type `dang_comp` and shared state used by all concrete components.
It includes common fields for amplitudes, templates, spectral indices, priors, polarization flags, and per-band mixing matrices, plus shared signal/spectrum helper routines.

In practice, all component-specific modules (CMB, dust, synchrotron-like models, templates) inherit from this base so they can be handled uniformly by samplers and solvers.
It provides the common contract that lets the rest of the pipeline evaluate any component without special-case logic.

The shared signal contract is the standard linear mixture form,

$$
\hat d_{\nu,p,k} = \sum_c a_{c,p,k}\,S_{c,\nu,p,k}(\theta_c) + \sum_t a_{t,\nu,k}\,T_{t,p,k},
$$

implemented through common evaluators so CG/sampling code can treat all component subclasses consistently.
