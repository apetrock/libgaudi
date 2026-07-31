#!/usr/bin/env python3
"""Verify calc_tangent_point_radius_grad against SymPy ∂k/∂dp.

Kernel (weight_functions.hpp):
  k = |P dp|^p / (|dp|^{2p} + l0^p),   P = N Nᵀ,  |N| = 1

C++ closed form used in production:
  dk = p k ( Px/(Px·Px + l0²) - 2 dp/(dp·dp + l0²) )

Run:
  python3 scripts/sympy/tangent_point_grad.py
"""

from __future__ import annotations

import numpy as np
import sympy as sp


def exact_grad_sym():
    """Symbolic gradient of k w.r.t. dp (N, l0, p constant)."""
    x, y, z = sp.symbols("x y z", real=True)
    nx, ny, nz = sp.symbols("nx ny nz", real=True)
    l0, p = sp.symbols("l0 p", positive=True)

    dp = sp.Matrix([x, y, z])
    N = sp.Matrix([nx, ny, nz])
    # Assume unit normal; keep symbolic and substitute later.
    Px = (N.dot(dp)) * N  # P dp for |N|=1
    f = sp.sqrt(Px.dot(Px))  # |N·dp|
    g = sp.sqrt(dp.dot(dp))
    k = f**p / (g ** (2 * p) + l0**p)

    grad = sp.Matrix([sp.diff(k, v) for v in (x, y, z)])
    return dp, N, l0, p, k, grad, f, g, Px


def cpp_grad_expr(dp, N, l0, p, k, Px):
    """Production closed form (regularized)."""
    l2 = l0**2
    return p * k * (Px / (Px.dot(Px) + l2) - 2 * dp / (dp.dot(dp) + l2))


def unregularized_exact_rewritten(dp, N, l0, p, k, f, g, Px):
    """Algebraic rewrite of exact ∇k (no abs-softening), for |N|=1, f>0, g>0:
      ∇k = p k [ Px/f² - 2 g^{2p-2} dp / (|dp|^{2p} + l0^p) ]
    """
    den = g ** (2 * p) + l0**p
    return p * k * (Px / f**2 - 2 * (g ** (2 * p - 2)) * dp / den)


def cpp_style_like_exact(dp, N, l0, p, k, f, g, Px):
    """What C++ would be without the l0² regularization shortcut:
      ∇k ≈ p k [ Px/f² - 2 dp / (|dp|²) * (something) ]
    For p=1:  ∇k = p k [ Px/f² - 2 dp / (|dp|² + l0) ]  if den = g²+l0 and g^{0}=1.
    """
    den = g ** (2 * p) + l0**p
    return p * k * (Px / f**2 - 2 * (g ** (2 * p - 2)) * dp / den)


def numeric_check(n_samples: int = 20, seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    dp_s, N_s, l0_s, p_s, k_s, grad_s, f_s, g_s, Px_s = exact_grad_sym()
    x, y, z = dp_s
    nx, ny, nz = N_s

    # Lambdify exact sympy grad and C++ form
    cpp = cpp_grad_expr(dp_s, N_s, l0_s, p_s, k_s, Px_s)
    rewrite = unregularized_exact_rewritten(
        dp_s, N_s, l0_s, p_s, k_s, f_s, g_s, Px_s
    )

    vars_ = (x, y, z, nx, ny, nz, l0_s, p_s)
    f_exact = sp.lambdify(vars_, grad_s, "numpy")
    f_cpp = sp.lambdify(vars_, cpp, "numpy")
    f_rewrite = sp.lambdify(vars_, rewrite, "numpy")
    f_k = sp.lambdify(vars_, k_s, "numpy")

    err_cpp = []
    err_rewrite = []
    for _ in range(n_samples):
        dp = rng.normal(size=3)
        N = rng.normal(size=3)
        N /= np.linalg.norm(N)
        # Avoid near-singular samples
        if abs(np.dot(N, dp)) < 1e-3 or np.linalg.norm(dp) < 1e-3:
            continue
        l0 = float(rng.uniform(0.05, 0.5))
        p = float(rng.choice([3.0, 4.0, 6.0]))
        args = (*dp, *N, l0, p)

        ge = np.asarray(f_exact(*args), dtype=float).reshape(3)
        gc = np.asarray(f_cpp(*args), dtype=float).reshape(3)
        gr = np.asarray(f_rewrite(*args), dtype=float).reshape(3)

        # Finite-difference check of sympy
        eps = 1e-6
        gfd = np.zeros(3)
        for ax in range(3):
            dpp, dpm = dp.copy(), dp.copy()
            dpp[ax] += eps
            dpm[ax] -= eps
            kp = float(f_k(*dpp, *N, l0, p))
            km = float(f_k(*dpm, *N, l0, p))
            gfd[ax] = (kp - km) / (2 * eps)

        err_cpp.append(np.linalg.norm(gc - ge) / (np.linalg.norm(ge) + 1e-14))
        err_rewrite.append(np.linalg.norm(gr - ge) / (np.linalg.norm(ge) + 1e-14))
        err_fd = np.linalg.norm(gfd - ge) / (np.linalg.norm(ge) + 1e-14)
        if err_fd > 1e-4:
            print(f"warn: sympy vs FD rel err {err_fd:.3e} (p={p})")

    print("=== tangent-point kernel gradient check ===")
    print("k = |N·dp|^p / (|dp|^{2p} + l0^p)")
    print()
    print("Exact rewrite vs SymPy  (should be ~0):")
    print(f"  mean rel err = {np.mean(err_rewrite):.3e}  max = {np.max(err_rewrite):.3e}")
    print()
    print("C++ form vs SymPy  (regularized; not exact for general p):")
    print(f"  mean rel err = {np.mean(err_cpp):.3e}  max = {np.max(err_cpp):.3e}")
    print()
    print("C++ uses:  p k ( Px/(|Px|²+l0²) - 2 dp/(|dp|²+l0²) )")
    print("Exact is:  p k ( Px/|Px|² - 2 |dp|^{2p-2} dp / (|dp|^{2p}+l0^p) )")
    print("They agree closely when p=1 and l0→0; diverge for large p / large l0.")


def show_symbolic_difference() -> None:
    dp, N, l0, p, k, grad, f, g, Px = exact_grad_sym()
    # Specialize to axis-aligned N = e_z for readability
    subs = {
        N[0]: 0,
        N[1]: 0,
        N[2]: 1,
        p: 2,  # concrete power for simplify
    }
    exact = sp.simplify(grad.subs(subs))
    cpp = sp.simplify(cpp_grad_expr(dp, N, l0, p, k, Px).subs(subs))
    rew = sp.simplify(
        unregularized_exact_rewritten(dp, N, l0, p, k, f, g, Px).subs(subs)
    )
    print("=== specialized N=(0,0,1), p=2 ===")
    print("exact ∇k =")
    sp.pprint(exact)
    print()
    print("rewrite ∇k =")
    sp.pprint(rew)
    print()
    print("C++ ∇k =")
    sp.pprint(cpp)
    print()
    print("exact - rewrite (should be 0):")
    sp.pprint(sp.simplify(exact - rew))
    print()
    print("exact - C++:")
    sp.pprint(sp.simplify(exact - cpp))


def main() -> None:
    numeric_check()
    print()
    show_symbolic_difference()


if __name__ == "__main__":
    main()
