# Dang Model Equations

This note ties the main equations to the current implementation so docs can be used as a consistency check against code behavior.

## 1) Forward model

For band $\nu$, pixel $p$, map/pol channel $k$:

$$
\hat d_{\nu p k} = \sum_{c \in \text{pixel-amp}} a_{c p k}\,S_{c,\nu p k}(\theta_c) + \sum_{t \in \text{template-like}} A_{t\nu k}\,T_{tpk}.
$$

Residual and weighted loss:

$$
r_{\nu p k} = d_{\nu p k} - \hat d_{\nu p k},
\qquad
\chi^2 = \sum_{\nu,p,k} \frac{r_{\nu p k}^2}{\sigma_{\nu p k}^2}.
$$

Code anchors: `src/dang_component_mod.f90`, `src/dang_data_mod.f90`.

## 2) Bandpass mapping

Channel-integrated response is represented as

$$
\langle S \rangle_b = \frac{\int R_b(\nu)S(\nu)\,d\nu}{\int R_b(\nu)\,d\nu},
$$

then converted between internal units via functions like `a2t`, `a2f`, `f2t`.

Code anchors: `src/dang_bp_mod.f90`, `src/dang_mixmat_1d_mod.f90`, `src/dang_mixmat_2d_mod.f90`.

## 3) Linear amplitude system (per CG group)

With nonlinear parameters fixed, amplitudes are solved from

$$
\left(\sum_\nu T_\nu^{\mathrm{T}}N_\nu^{-1}T_\nu\right)x = \sum_\nu T_\nu^{\mathrm{T}}N_\nu^{-1}d_\nu^*
$$

where $d_\nu^*$ is data after subtracting components not included in the active CG solve.

Matrix-free operator application in code is

$$
x \to T_\nu x \to N_\nu^{-1}(T_\nu x) \to T_\nu^{\mathrm{T}}N_\nu^{-1}T_\nu x,
\quad
Ax=\sum_\nu(\cdot).
$$

Code anchors: `src/dang_cg_mod.f90` (`compute_rhs`, `compute_Ax`, `cg_search`).

## 4) Sampling step (nonlinear parameters)

Metropolis-Hastings acceptance uses

$$
\alpha = \min\left(1,\exp\left[\ln P(\theta'\mid d)-\ln P(\theta\mid d)\right]\right),
\qquad
\ln P = -\tfrac12\chi^2 + \ln \pi(\theta) + C.
$$

In CG sampling mode, RHS is perturbed by

$$
b \leftarrow b + \sum_\nu T_\nu^{\mathrm{T}}N_\nu^{-1/2}\eta,
\qquad \eta\sim\mathcal N(0,I).
$$

Code anchors: `src/dang_sample_mod.f90`, `src/dang_lnl_mod.f90`, `src/dang_cg_mod.f90`.

## 5) Component SED patterns in use

Power law:

$$
S(\nu;\beta)=\left(\frac{\nu}{\nu_{\mathrm{ref}}}\right)^\beta.
$$

Modified blackbody (dust-like):

$$
S(\nu;\beta,T)\propto \left(\frac{\nu}{\nu_{\mathrm{ref}}}\right)^\beta \frac{B_\nu(T)}{B_{\nu_{\mathrm{ref}}}(T)}.
$$

Template-like components:

$$
m_{\nu p} = A_\nu\,t_p,
$$

and monopole special case $t_p=1$.

Code anchors: `src/dang_powlaw_comp_mod.f90`, `src/dang_mbb_comp_mod.f90`, `src/dang_template_comp_mod.f90`, `src/dang_monopole_comp_mod.f90`.

## Consistency checks

Use this as a quick intended-vs-actual checklist during doc/code updates.

1. **Data subtraction before CG RHS**
   - Expected: $d_\nu^* = d_\nu - \sum_{c\notin g}\hat d_{c,\nu}$ before building $b$.
   - Check: `compute_rhs` subtracts non-group and non-sampled components.

2. **Noise weighting location**
   - Expected: one factor of $\sigma^{-2}$ in `compute_rhs`, and one in `compute_Ax` after $T_\nu x$.
   - Check: divisions by `rms_map**2` appear in those exact stages.

3. **Template-like parameterization**
   - Expected: `template`/`monopole`/`hi_fit` use per-band coefficients, not per-pixel amplitude vectors.
   - Check: vector offsets advance by `c%nfit` for these component types.

4. **Bandpass consistency across components**
   - Expected: single-$\nu$ and integrated-band evaluation paths map to the same physical SED definition.
   - Check: component `eval` methods and mixmat precompute/integration routines agree.

5. **Pol-flag branch reachability (important code-level check)**
   - Intended behavior includes modes like `T`, `Q`, `U`, `Q+U`, `T+Q+U`.
   - Check: conditions using `iand(flag,0)` are mathematically always zero; if used as a branch test, that branch is unreachable and should be reviewed against intent.