# dang_sample_mod

Implements sampling workflows for model parameter inference.
It includes Metropolis-Hastings updates for spectral parameters (full-sky or pixel-wise), calibration fitting routines, model update helpers, and MH step-size tuning logic.

In practice, this is the control loop for nonlinear parameter updates between linear amplitude solves.
It proposes new spectral/calibration parameters, evaluates acceptance using likelihood routines, and adapts proposal scales for stable chain behavior.

For MH updates it follows

$$
\alpha = \min\left(1,\exp\left[\ln P(\theta'\mid d)-\ln P(\theta\mid d)\right]\right),
$$

accepting proposals with probability \(\alpha\), while alternating with CG amplitude solves conditioned on current nonlinear parameters.
