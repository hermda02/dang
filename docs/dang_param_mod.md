# dang_param_mod

Defines the `dang_params` configuration structure and the parameter-file parsing pipeline.
It handles include expansion and keyed lookups, then maps global/data/component/CG settings into structured fields used by constructors and samplers.

In practice, this module is the run configuration loader: it turns text parameter files into strongly structured runtime settings.
Most other modules depend on these parsed fields to know what data to read, what components to build, and how to sample.

Operationally, it performs a mapping

$$
\texttt{key=value text} \longrightarrow \texttt{dang\_params fields},
$$

including include expansion and component-wise arrays such as \(\{\texttt{fg\_type}_i,\texttt{fg\_cg\_group}_i,\ldots\}_{i=1}^{n_{\mathrm{comp}}}\).
