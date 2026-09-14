# dang_bp_mod

Handles bandpass metadata and operations through `bandinfo` and the global `bp` array.
It reads and normalizes transmission curves, integrates spectra over bandpasses, and provides unit conversion helpers (such as `a2t`, `a2f`, and `f2t`) plus derivative utilities used across the codebase.
