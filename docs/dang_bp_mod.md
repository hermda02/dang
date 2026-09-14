# dang_bp_mod

Handles bandpass metadata and operations through `bandinfo` and the global `bp` array.
It reads and normalizes transmission curves, integrates spectra over bandpasses, and provides unit conversion helpers (such as `a2t`, `a2f`, and `f2t`) plus derivative utilities used across the codebase.

In practice, whenever a component is evaluated for a real instrument channel, this module converts a theoretical spectrum into what that channel actually measures.
It is the bridge between physics-level SED formulas and detector-level map units.

Operationally, channel response is handled as a normalized band integral,

$$
\langle S \rangle_b = \frac{\int R_b(\nu)\,S(\nu)\,d\nu}{\int R_b(\nu)\,d\nu},
$$

with conversion factors then applied between amplitude, flux, and temperature units (via `a2t`, `a2f`, `f2t`).
