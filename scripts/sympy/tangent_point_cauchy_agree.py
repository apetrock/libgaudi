#!/usr/bin/env python3
"""Circle-agreement term EΔ for Cauchy TP pairs.

Geometric TP circles (N held fixed in ∇; caller orients chart):
  Ri = |dp|² / (2 Ni·dp)     with Ni·dp > 0,  ci = pi + Ri Ni
  Rj = |dp|² / (2 Nj·(-dp))  with Nj·(-dp) > 0, cj = pj + Rj Nj
  δc = ci - cj = -dp + Ri Ni - Rj Nj
  EΔ = (1/2) |δc|² / (|dp|² + eps²)

EΔ = 0 when the two triad circles coincide.

Closed-form ∇_dp EΔ (Ni, Nj fixed):
  ∇Ri = dp/ai - Ri Ni/ai,   ai = Ni·dp
  ∇Rj = dp/aj + Rj Nj/aj,   aj = Nj·(-dp) = -Nj·dp
  J^T δc = -δc + (Ni·δc) ∇Ri - (Nj·δc) ∇Rj
  ∇EΔ = J^T δc / D - |δc|² dp / D²

Run:
  python3 scripts/sympy/tangent_point_cauchy_agree.py
  python3 scripts/sympy/tangent_point_cauchy_agree.py --dump-hess
"""

from __future__ import annotations

import argparse

import numpy as np
import sympy as sp


def agree_sym():
    x, y, z = sp.symbols("x y z", real=True)
    nix, niy, niz = sp.symbols("nix niy niz", real=True)
    njx, njy, njz = sp.symbols("njx njy njz", real=True)
    eps = sp.symbols("eps", positive=True)
    dp = sp.Matrix([x, y, z])
    Ni = sp.Matrix([nix, niy, niz])
    Nj = sp.Matrix([njx, njy, njz])
    ai = Ni.dot(dp)          # > 0 on chart
    aj = -Nj.dot(dp)         # Nj·(-dp) > 0 on chart
    r2 = dp.dot(dp)
    Ri = r2 / (2 * ai)
    Rj = r2 / (2 * aj)
    dc = -dp + Ri * Ni - Rj * Nj
    D = r2 + eps**2
    E = sp.Rational(1, 2) * dc.dot(dc) / D
    grad = sp.Matrix([sp.diff(E, v) for v in (x, y, z)])

    # Closed form: ∂ai/∂dp = Ni, ∂aj/∂dp = -Nj
    # ∇Ri = (2 dp ai - r2 Ni)/(2 ai²) = dp/ai - Ri Ni/ai
    # ∇Rj = (2 dp aj - r2 (-Nj))/(2 aj²) = dp/aj + Rj Nj/aj
    gRi = dp / ai - Ri * Ni / ai
    gRj = dp / aj + Rj * Nj / aj
    JTdc = -dc + (Ni.dot(dc)) * gRi - (Nj.dot(dc)) * gRj
    cpp = JTdc / D - (dc.dot(dc)) * dp / D**2
    return dp, Ni, Nj, eps, E, grad, cpp, D, dc


def _orient_pair(dpv, Ni, Nj):
    """Flip normals onto the geometric chart."""
    if Ni @ dpv < 0:
        Ni = -Ni
    if Nj @ (-dpv) < 0:
        Nj = -Nj
    return Ni, Nj


def numeric_check(n_samples: int = 40, seed: int = 0) -> None:
    dp_s, Ni_s, Nj_s, eps_s, E_s, grad_s, cpp_s, _, _ = agree_sym()
    x, y, z = dp_s
    nix, niy, niz = Ni_s
    njx, njy, njz = Nj_s
    vars_ = (x, y, z, nix, niy, niz, njx, njy, njz, eps_s)
    f_exact = sp.lambdify(vars_, grad_s, "numpy")
    f_cpp = sp.lambdify(vars_, cpp_s, "numpy")
    f_E = sp.lambdify(vars_, E_s, "numpy")
    f_H = sp.lambdify(vars_, sp.hessian(E_s, [x, y, z]), "numpy")

    rng = np.random.default_rng(seed)
    max_g = 0.0
    max_h = 0.0
    n_ok = 0
    for _ in range(n_samples):
        dpv = rng.normal(size=3)
        Ni = rng.normal(size=3)
        Nj = rng.normal(size=3)
        Ni /= np.linalg.norm(Ni)
        Nj /= np.linalg.norm(Nj)
        Ni, Nj = _orient_pair(dpv, Ni, Nj)
        if Ni @ dpv < 0.05 or Nj @ (-dpv) < 0.05:
            continue
        ev = float(rng.uniform(0.05, 0.5))
        args = (*dpv, *Ni, *Nj, ev)
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
            gp = np.asarray(f_exact(*dpp, *Ni, *Nj, ev), dtype=float).ravel()
            gm = np.asarray(f_exact(*dpm, *Ni, *Nj, ev), dtype=float).ravel()
            Hfd[:, i] = (gp - gm) / (2 * h)
        max_h = max(max_h, float(np.max(np.abs(Hn - Hfd))))
        n_ok += 1

    print(f"max |∇_exact - ∇_cpp| over {n_ok}: {max_g:.3e}")
    print(f"max |H_exact - H_FD|   over {n_ok}: {max_h:.3e}")

    # Coincident circle: pi=0, pj=ex, c=(0.5,0.5,0)
    dpv = np.array([1.0, 0.0, 0.0])
    c = np.array([0.5, 0.5, 0.0])
    R = np.linalg.norm(c)
    Ni = c / R
    Nj = (c - dpv) / R
    Ni, Nj = _orient_pair(dpv, Ni, Nj)
    e = float(f_E(*dpv, *Ni, *Nj, 0.1))
    print(f"sample EΔ coincident circle: {e:.6e}")


def dump_hessian_cpp() -> str:
    dp, Ni, Nj, eps, E, _, _, _, _ = agree_sym()
    x, y, z = dp
    nix, niy, niz = Ni
    njx, njy, njz = Nj
    H = sp.hessian(E, [x, y, z])
    flat = [H[i, j] for i in range(3) for j in range(3)]
    repl, red = sp.cse(flat, symbols=sp.numbered_symbols("x"))
    lines = [
        "inline mat3 calc_circle_agree_hessian_chart(",
        "    const vec3 &dp, const vec3 &Ni, const vec3 &Nj, const real &eps) {",
        "  // generated: python3 scripts/sympy/tangent_point_cauchy_agree.py --dump-hess",
        "  // Chart: Ni·dp > 0, Nj·(-dp) > 0.",
        "  const real x = dp[0], y = dp[1], z = dp[2];",
        "  const real nix = Ni[0], niy = Ni[1], niz = Ni[2];",
        "  const real njx = Nj[0], njy = Nj[1], njz = Nj[2];",
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
        help="print pasteable circle-agree Hessian C++",
    )
    args = ap.parse_args()
    if args.dump_hess:
        print(dump_hessian_cpp())
        return
    _, _, _, _, _, grad, cpp, _, _ = agree_sym()
    print("symbolic ∇EΔ − cpp (should be 0):")
    sp.pprint(sp.simplify(grad - cpp))
    print()
    numeric_check()


if __name__ == "__main__":
    main()
