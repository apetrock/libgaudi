#!/usr/bin/env python3
"""Convolution for harmonic / mollified TP: Φ = κ ∗ T (not composition).

Leaf kernels (separable):
  T(dp; N, l0, p)  — regularized TP density (tangent_point_grad.py)
  κ_ℓ              — Gaussian mollifier (ambient)

Discrete convolution at sample i:
  Φ_i = ∑_j κ(π_i − π_j) T_j

Derivative of the convolution (product rule on each κ_ij T_j):
  ∇Φ ≈ κ ∗ ∇T  +  (∇κ) ∗ T     =  Gs + Ks

Hessian in direction v (compose leaf Hessians):
  H_Φ v = κ ∗ (H_T v)
        + (Hess_κ · δdp) ∗ T
        + (δκ) ∗ ∇T
        + (∇κ) ∗ (G·v)          # δT_j ≈ G_j·v_j

where δdp = v_j − v_i,  δκ = (∂κ/∂dp)·δdp,  ∇_i κ = −∂κ/∂dp.

Gaussian leaf (matches weight_functions.hpp):
  κ = (C/ℓ) exp(−|dp|²/ℓ²),  C = 1/√(2π)
  ∂κ/∂dp = −(2/ℓ²) κ dp
  Hess_dp κ = −(2/ℓ²) κ I + (4/ℓ⁴) κ (dp ⊗ dp)

Soft TP is composition (softfloor ∘ R_tp); harmonic is this convolution.

Run:
  python3 scripts/sympy/tp_harmonic_compose.py
  python3 scripts/sympy/tp_harmonic_compose.py --dump-hess-kappa
"""

from __future__ import annotations

import argparse

import numpy as np
import sympy as sp


def gaussian_leaf():
    x, y, z = sp.symbols("x y z", real=True)
    ell = sp.symbols("ell", positive=True)
    C = 1 / sp.sqrt(2 * sp.pi)
    dp = sp.Matrix([x, y, z])
    r2 = dp.dot(dp)
    kappa = (C / ell) * sp.exp(-r2 / ell**2)
    dkappa = sp.Matrix([sp.diff(kappa, v) for v in (x, y, z)])
    H = sp.hessian(kappa, [x, y, z])
    return dp, ell, kappa, dkappa, H


def closed_forms():
    dp, ell, kappa, dkappa, H = gaussian_leaf()
    # Expected: dκ = -(2/ℓ²) κ dp
    expect_d = sp.simplify(-(2 / ell**2) * kappa * dp)
    # Expected: H = -(2/ℓ²)κ I + (4/ℓ⁴)κ dp⊗dp
    I = sp.eye(3)
    expect_H = sp.simplify(
        -(2 / ell**2) * kappa * I + (4 / ell**4) * kappa * (dp * dp.T)
    )
    print("=== Gaussian mollifier leaf ===")
    print("κ = (C/ℓ) exp(−|dp|²/ℓ²)")
    print()
    print("∂κ/∂dp − (−2/ℓ² κ dp)  (should be 0):")
    sp.pprint(sp.simplify(dkappa - expect_d))
    print()
    print("Hess κ − closed form (should be 0):")
    sp.pprint(sp.simplify(H - expect_H))
    print()
    print("Convolution (Φ = κ ∗ T):")
    print("  ∇Φ  = κ∗∇T + (∇κ)∗T")
    print("  Hv  = κ∗(H_T v) + (Hess_κ·δdp)∗T + (δκ)∗∇T + (∇κ)∗(G·v)")
    print("Leaf T / H_T: regularized TP (tangent_point_grad.py --dump-hess)")
    print("Soft TP is composition; harmonic is convolution.")


def dump_hess_kappa_cpp() -> str:
    """Pasteable closed form (already simple — emit the closed form, not CSE soup)."""
    return """inline mat3 calc_d2_gaussian(const vec3 &dp, real l) {
  // Hess_dp of κ = (C/l) exp(−|dp|²/l²),  C = 1/√(2π)
  // generated/verified: python3 scripts/sympy/tp_harmonic_compose.py
  const real lp = l * l;
  const real kappa = calc_gaussian(dp, l);
  const real c0 = real(-2.0) * kappa / lp;
  const real c1 = real(4.0) * kappa / (lp * lp);
  mat3 H = c0 * mat3::Identity();
  H += c1 * (dp * dp.transpose());
  return H;
}"""


def numeric_hess_check(n: int = 30, seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    dp_s, ell_s, kappa_s, dkappa_s, H_s = gaussian_leaf()
    x, y, z = dp_s
    fH = sp.lambdify((x, y, z, ell_s), H_s, "numpy")
    fk = sp.lambdify((x, y, z, ell_s), kappa_s, "numpy")
    errs = []
    for _ in range(n):
        dp = rng.normal(size=3)
        ell = float(rng.uniform(0.05, 0.5))
        H = np.asarray(fH(*dp, ell), dtype=float).reshape(3, 3)
        # FD Hessian
        eps = 1e-5
        Hfd = np.zeros((3, 3))
        for a in range(3):
            for b in range(3):
                def k_at(d):
                    return float(fk(*d, ell))
                dpp, dpm = dp.copy(), dp.copy()
                dpp[a] += eps
                dpp[b] += eps
                dmm = dp.copy()
                dmm[a] -= eps
                dmm[b] -= eps
                dpm = dp.copy()
                dpm[a] += eps
                dpm[b] -= eps
                dmp = dp.copy()
                dmp[a] -= eps
                dmp[b] += eps
                Hfd[a, b] = (k_at(dpp) - k_at(dpm) - k_at(dmp) + k_at(dmm)) / (
                    4 * eps * eps
                )
        denom = np.linalg.norm(H) + 1e-14
        errs.append(np.linalg.norm(H - Hfd) / denom)
    print("=== Hess_κ vs FD ===")
    print(f"  mean rel err = {np.mean(errs):.3e}  max = {np.max(errs):.3e}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-hess-kappa", action="store_true")
    args = ap.parse_args()
    if args.dump_hess_kappa:
        print(dump_hess_kappa_cpp())
        return
    closed_forms()
    print()
    numeric_hess_check()


if __name__ == "__main__":
    main()
