# dang_cmb_comp_mod

Defines the CMB anisotropy component as a `dang_comp` extension with no sampled spectral indices.
Its evaluation returns the RJ/CMB conversion scaling (or bandpass-integrated equivalent) in microkelvin units.

In practice, this is the baseline sky component for CMB fluctuations: it contributes a fixed spectral shape while only the amplitude map is solved/sampled.
That makes it a stable anchor term in multi-component fits against foregrounds.

For each band, the component contribution is effectively

$$
m_{\nu,p} = a_{\mathrm{cmb},p}\,S_{\nu}^{\mathrm{cmb}},
$$

where `S_\nu^{cmb}` is the RJ/CMB conversion factor (single-\(\nu\) or bandpass-integrated) and only the map amplitude \(a_{\mathrm{cmb},p}\) is fitted.
