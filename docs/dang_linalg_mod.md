# dang_linalg_mod

Contains core linear algebra utilities used by solvers and preconditioners.
It includes Cholesky/LU-style routines, substitutions, sparse-like format conversions, matrix-vector products, and helper products such as `A^T A`.

In practice, this is the numerical backbone underneath CG and likelihood calculations.
It packages the low-level matrix operations needed for speed and memory efficiency in large map-based systems.

Typical usage is factoring and solving systems like

$$
A x = b,\qquad A = L L^{\mathsf T}\;\text{(Cholesky)}
$$

then applying triangular solves, sparse transforms, and matrix products needed by the solver stack.
