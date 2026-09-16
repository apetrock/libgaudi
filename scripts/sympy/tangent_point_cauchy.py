#!/usr/bin/env python3
"""Exact ∇ / Hess of Cauchy TP density vs C++ closed forms.

  K = |N·dp|^p / (|dp|² + ε²)^p

On the N·dp>0 chart (|N|=1):
  K = (N·dp)^p / (|dp|² + ε²)^p
  ∇K = p K ( Px / |Px|² - 2 dp / (|dp|²+ε²) ),  Px = (N·dp) N

ε already regularizes the kernel — no extra leaf-grad smoothing.

Run:
  python3 scripts/sympy/tangent_point_cauchy.py
  python3 scripts/sympy/tangent_point_cauchy.py --dump-hess
"""

from __future__ import annotations

import argparse

import numpy as np
import sympy as sp


def cauchy_sym():
    x, y, z = sp.symbols("x y z", real=True)
    nx, ny, nz = sp.symbols("nx ny nz", real=True)
    eps, p = sp.symbols("eps p", positive=True)
    dp = sp.Matrix([x, y, z])
    N = sp.Matrix([nx, ny, nz])
    f = N.dot(dp)
    R2 = dp.dot(dp) + eps**2
    K = f**p / R2**p
    grad = sp.Matrix([sp.diff(K, v) for v in (x, y, z)])
    Px = f * N
    cpp = p * K * (Px / f**2 - 2 * dp / R2)
    return dp, N, eps, p, K, grad, cpp, R2


def numeric_check(n_samples: int = 30, seed: int = 0) -> None:
    dp_s, N_s, eps_s, p_s, K_s, grad_s, cpp_s, _ = cauchy_sym()
    x, y, z = dp_s
    nx, ny, nz = N_s
    vars_ = (x, y, z, nx, ny, nz, eps_s, p_s)
    f_exact = sp.lambdify(vars_, grad_s, "numpy")
    f_cpp = sp.lambdify(vars_, cpp_s, "numpy")
    f_K = sp.lambdify(vars_, K_s, "numpy")
    f_H = sp.lambdify(vars_, sp.hessian(K_s, [x, y, z]), "numpy")

    rng = np.random.default_rng(seed)
    max_g = 0.0
    max_h = 0.0
    for _ in range(n_samples):
        dpv = rng.normal(size=3)
        Nv = rng.normal(size=3)
        Nv /= np.linalg.norm(Nv)
        if Nv @ dpv < 0:
            Nv = -Nv
        ev, pv = float(rng.uniform(0.05, 0.5)), float(rng.choice([2, 4, 6]))
        args = (*dpv, *Nv, ev, pv)
        ge = np.asarray(f_exact(*args), dtype=float).ravel()
        gc = np.asarray(f_cpp(*args), dtype=float).ravel()
        max_g = max(max_g, float(np.max(np.abs(ge - gc))))

        Hn = np.asarray(f_H(*args), dtype=float)
        h = 1e-6
        Hfd = np.zeros((3, 3))
        for i in range(3):
            dpp = dpv.copy()
            dpm = dpv.copy()
            dpp[i] += h
            dpm[i] -= h
            gp = np.asarray(f_exact(*dpp, *Nv, ev, pv), dtype=float).ravel()
            gm = np.asarray(f_exact(*dpm, *Nv, ev, pv), dtype=float).ravel()
            Hfd[:, i] = (gp - gm) / (2 * h)
        max_h = max(max_h, float(np.max(np.abs(Hn - Hfd))))

    print(f"max |∇_exact - ∇_cpp| over {n_samples}: {max_g:.3e}")
    print(f"max |H_exact - H_FD|   over {n_samples}: {max_h:.3e}")
    print(f"sample K(1,0,0; ez): {f_K(1, 0, 0, 0, 0, 1, 0.1, 6):.6e}")


def dump_hessian_cpp() -> str:
    dp, N, eps, p, K, _, _, _ = cauchy_sym()
    x, y, z = dp
    nx, ny, nz = N
    H = sp.hessian(K, [x, y, z])
    flat = [H[i, j] for i in range(3) for j in range(3)]
    repl, red = sp.cse(flat, symbols=sp.numbered_symbols("x"))
    lines = [
        "inline mat3 calc_tangent_point_radius_hessian_cauchy_chart(",
        "    const vec3 &dp, const vec3 &N, const real &eps, const real &p) {",
        "  // generated: python3 scripts/sympy/tangent_point_cauchy.py --dump-hess",
        "  const real x = dp[0], y = dp[1], z = dp[2];",
        "  const real nx = N[0], ny = N[1], nz = N[2];",
    ]
    for s, e in repl:
        lines.append(f"  const real {s} = {sp.cxxcode(e)};")
    lines.append("  mat3 H;")
    idx = 0
    for i in range(3):
        for j in range(3):
            lines.append(f"  H({i}, {j}) = {sp.cxxcode(red[idx])};")
            idx += 1
    lines.append(
        "  if (!std::isfinite(H(0, 0)) || !std::isfinite(H(1, 1)) || "
        "!std::isfinite(H(2, 2)))"
    )
    lines.append("    return mat3::Zero();")
    lines.append("  return H;")
    lines.append("}")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--dump-hess",
        action="store_true",
        help="print pasteable Cauchy Hessian C++",
    )
    args = ap.parse_args()
    if args.dump_hess:
        print(dump_hessian_cpp())
        return
    _, _, _, _, _, grad, cpp, _ = cauchy_sym()
    print("symbolic ∇K − cpp (should be 0):")
    sp.pprint(sp.simplify(grad - cpp))
    print()
    numeric_check()


if __name__ == "__main__":
    main()
