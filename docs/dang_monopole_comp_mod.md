# dang_monopole_comp_mod

Implements a monopole template component with fixed intensity template behavior and zero polarization template.
It supports per-band fitted offsets and evaluates to active/inactive template contributions depending on band-correlation settings.

In practice, this module captures constant map offsets per band, helping absorb large-scale zero-level mismatches that are not true sky structure.
That prevents those offsets from biasing other astrophysical component amplitudes.

Its template is spatially constant in intensity, so for correlated bands,

$$
m_{\nu,p} = a_\nu\cdot 1,
$$

with \(a_\nu\) solved as a per-band offset parameter.
