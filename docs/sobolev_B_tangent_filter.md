# Sobolev B matvec

`B` and `B0` share the residual in Appendix B.3.4 Eq. 31 of Yu, Schumacher, Crane, *Repulsive Curves* (`reference/crane2021-repulsive-curves.pdf`):

```text
r_I = Σ_J w_IJ (φ_I − φ_J)
```

`B` lifts each coordinate by the unit tangent. Section 5.2.1 sets \(D_I u = ((u_\mathrm{end} - u_\mathrm{start})/\ell)\,T\), and Equation 21 builds the block from the inner product of those edge vectors, so a cross term carries \(\langle T_I, T_J\rangle\). The matvec applies \(r_I = \sum w_{IJ}(\varphi_I - \varphi_J)\) to \(\varphi = (\Delta u/\ell)\,T\) and scatters \((T\cdot r)/\ell\). `B0` lifts with the edge average. Eq. 19 copies the same scalar matrix onto x, y, and z. Kernel lengths are geometric `ℓ`. The metric is `A = α_B B + α_0 B_0`. Translations are a kernel of both blocks; the barycenter constraint removes them.

A leaf partner is that edge's `φ_J`. A far node in Appendix B.1 contributes `Σ_J ℓ_J φ_J`, one factor of length. The covariance pyramid stores sums, and `w` already contains the cluster length `L = Σ ℓ`, so the partner is `(w / L) Σ (ℓ φ)`. Using the raw sum of `φ` multiplies the partner by the cluster's edge count.

Test: `sobolev_far_constant_field_is_kernel` in [`include/gaudi/test/tp_sobolev_jacobi_tests.hpp`](../include/gaudi/test/tp_sobolev_jacobi_tests.hpp). A constant vertex field has a constant edge average, so `(φ_I − φ_J)` is zero and `B0 u` is zero. A summed-`φ` partner is not.

The scalar kernel weight is not multiplied by \(\langle T_I, T_J\rangle\). That product is already inside \(\varphi = (\Delta u/\ell)\,T\). A perpendicular partner is not dropped from the kernel.

Test: `sobolev_B_eq31_residual` and `sobolev_matvec_matches_eq31` in [`include/gaudi/test/tp_sobolev_jacobi_tests.hpp`](../include/gaudi/test/tp_sobolev_jacobi_tests.hpp).
