# dang_data_mod

Provides the central `dang_data` container for maps, residuals, masks, gains, unit conversions, and I/O metadata.
Key routines initialize/read/convert datasets, update model and residual maps from components, compute chi-square diagnostics, and write outputs for analysis.

In practice, this module is the state hub of a run: observed maps come in here, model maps and residuals are updated here, and diagnostics/products are written from here.
Most high-level iteration steps read from and write to `dang_data`.

Its core bookkeeping follows

$$
r_{\nu,p,k} = d_{\nu,p,k} - \hat d_{\nu,p,k},
\qquad
\chi^2 = \sum_{\nu,p,k}\frac{r_{\nu,p,k}^2}{\sigma_{\nu,p,k}^2},
$$

which is what downstream likelihood and convergence diagnostics consume.
