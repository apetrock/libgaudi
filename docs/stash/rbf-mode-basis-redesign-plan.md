# STASHED: Construct/Init Split + RBF Mode-Basis Redesign

> Stashed 2026-08-31. **Superseded by** [`../rbf-mode-basis-redesign.md`](../rbf-mode-basis-redesign.md).

---

## User correction (supersedes prior RBF fixes)

The **variational bulk energy is unchanged** (`|∇ΔD|² + λ_gap ΔD²`). The original method was right:

- **Same geometry access** for rest embed (`Q_r`) and live embed (`Q_k`)
- **`ΔD = D(Q_r,x) − D(Q_k,x)`** with both fields evaluated on the **same foot-local samples**

Recent RBF patches (fixed rest centers passed to live refit, inhomogeneous LS as primary gauge) **did not fix the real problem** and should be **backed out** once the replacement is in place.

### Root cause reframing (not "pick perfect interpolant")

| Issue | Wrong framing | Correct framing |
|-------|---------------|-----------------|
| Identity at rest fails | SVD null space / need tolerance skip | **Representation mismatch** between `Q_r` and `Q_k` |
| Wild response to small bend | Too many RBF centers / bad subsampling | **Over-parameterized MLS** |
| Singular / unstable fits | Need homogeneous SVD | **Regularized low-order basis** for probable deformation modes |

**Design principle:** RBF field is a **smooth low-order prior** (cubic/quartic), not a perfect interpolant.

---

## What to revert (RBF band-aids)

In `include/gaudi/hepworth/block/local_slab_field_fit.hpp`:

- `fit_rbf_mls_on_stencil(..., fixed_centers*, fixed_sigma)` and `rest_basis` on `fit_field_coeffs_on_stencil`
- Live refit passing `&_Q_rest` only to freeze centers

**Keep:** ring-only hub normals, stack MLS matrices, no energy-tolerance skip_gd.

---

## New RBF / MLS design (sketch)

```text
D(Q, x) = p(x) + Σ_k w_k φ_k(x)
```

- **Shared basis** at Q_r: poly order, RBF centers, σ frozen
- **Q_k refit**: coefficients only on same basis
- Target: cubic (10) + 3–4 wide RBF modes, Tikhonov regularization
- Same `pos_at` from blocks (requires construct/init split)

---

## Implementation sequence (when un-stashed)

1. Revert RBF band-aids
2. Phase 0 + 1 construct/init split
3. New regularized low-order RBF fit
4. Identity + small/large bend tests
5. Demo default Darboux until RBF passes

See full detail in git history of `.cursor/plans/Construct init split migration-4b5aaf09.plan.md` from this date.
