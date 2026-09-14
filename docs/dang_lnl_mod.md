# dang_lnl_mod

Implements likelihood and prior evaluation routines used by samplers.
It provides chi-square and marginalized-amplitude log-likelihood functions, combined posterior wrappers, and specialized prior evaluators (including Jeffreys-style cases).

In practice, this module decides whether a proposed parameter move is statistically better or worse.
The sampling code calls these routines repeatedly during Metropolis-Hastings steps.

The core scoring is standard,

$$
\ln \mathcal L(\theta) = -\tfrac12\chi^2(\theta)+C,
\qquad
\chi^2=\sum\frac{(d-\hat d(\theta))^2}{\sigma^2},
$$

with posterior terms formed as \(\ln P(\theta\mid d)=\ln\mathcal L(\theta)+\ln\pi(\theta)\).
