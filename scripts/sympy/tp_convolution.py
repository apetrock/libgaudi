#!/usr/bin/env python3
"""Can we compute K_ε = (K_TP * ζ_ε) with SymPy / numeric tools?

Singular TP kernel (N = e_z fixed):
  K(x) = |z|^p / |x|^{2p}

Gaussian blob:
  ζ_ε(y) = (2π ε²)^{-3/2} exp(-|y|²/(2ε²))

  K_ε(x) = ∫ K(x-y) ζ_ε(y) dy

Verdict from practice:
  - SymPy: can *write* the integral; full 3D Abs/power closed form usually hangs.
  - SymPy: sometimes reduces the ρ piece in cylindrical coords for small p.
  - Reliable path: numeric quadrature (or Fourier/Hankel) for K_ε(r), then fit m(ρ).

Run:
  python3 scripts/sympy/tp_convolution.py
"""

from __future__ import annotations

import numpy as np


def Ke_monte_carlo(r: float, eps: float, p: float, n: int = 200_000, seed: int = 0):
    """Monte Carlo: y ~ ζ_ε, estimate E[K(x-y)] with x=(0,0,r)."""
    rng = np.random.default_rng(seed)
    y = rng.normal(scale=eps, size=(n, 3))  # N(0,ε² I) has density ζ_ε
    v = np.array([0.0, 0.0, r]) - y
    r2 = np.einsum("ij,ij->i", v, v)
    # K=0 at coincidence; elsewhere |vz|^p / r2^p
    mask = r2 > 1e-32
    K = np.zeros(n)
    K[mask] = (np.abs(v[mask, 2]) ** p) / (r2[mask] ** p)
    return float(np.mean(K)), float(np.std(K) / np.sqrt(n))


def main() -> None:
    eps = 0.1
    p = 6.0

    print("=== TP mollification by convolution K_ε = K * ζ_ε ===")
    print("N = e_z, x = r e_z, ζ = Gaussian(0, ε^2 I)")
    print()
    print("SymPy status:")
    print("  - Can express the triple integral for K*ζ")
    print("  - Full symbolic integrate typically hangs (Abs + high powers).")
    print("  - So: use numeric (K*ζ), optionally fit radial factor m(r/ε).")
    print()

    print(f"{'r/ε':>8s} {'K_ε (MC)':>14s} {'±stderr':>10s} "
          f"{'k_reg along-N':>14s} {'1/r^p':>12s}")

    for scale in [0.05, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]:
        r = scale * eps
        Ke, se = Ke_monte_carlo(r, eps, p)
        k_reg = (r**p) / (r ** (2 * p) + eps**p)
        k_raw = r ** (-p)
        print(f"{scale:8.2f} {Ke:14.4e} {se:10.2e} {k_reg:14.4e} {k_raw:12.4e}")

    print()
    print("Interpretation:")
    print("  K_ε stays finite as r/ε → 0 (true blob mollification).")
    print("  That IS the integral (f*g). SymPy sets it up; numeric evaluates it.")
    print("  Next step if useful: tabulate K_ε(ρ)/(|N·û|^p/r^p) → m(ρ) for C++.")


if __name__ == "__main__":
    main()
