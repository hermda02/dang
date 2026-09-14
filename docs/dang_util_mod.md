# dang_util_mod

Provides shared constants, global state, and general utilities used across the project.
It includes physical constants, RNG/priors helpers, string and flag parsers, masking/map I/O helpers, and infrastructure for MPI/OpenMP-aware workflows.

In practice, this is the common toolbox every other module leans on.
It centralizes low-level helpers and global run context so core science modules can stay focused on modeling and inference logic.

It also hosts common probability utilities used across modules, e.g. Gaussian priors of the form

$$
\pi(x\mid\mu,\sigma)\propto \exp\!\left[-\frac{(x-\mu)^2}{2\sigma^2}\right],
$$

plus RNG and masking helpers that keep sampling and map operations consistent.
