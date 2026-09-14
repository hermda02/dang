# dang_cg_mod

Implements conjugate-gradient machinery for solving or sampling grouped component amplitudes.
It builds right-hand-side vectors, applies the implicit normal-equation operator (`compute_Ax`), runs CG iterations, and maps solved values back into component amplitude/template structures.
