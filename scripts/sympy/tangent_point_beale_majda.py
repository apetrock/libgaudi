#!/usr/bin/env python3
"""Beale–Majda mollification for the tangent-point kernel (SymPy + numeric).

Brochu/Keeler/Bridson SCA 2012 (vortex sheet smoke) mollify Biot–Savart as:
  grad K_moll(r) = -(x-y) * (1 - exp(-(r/β)^3)) / r^3
so the factor (1 - exp(-(r/β)^m)) cancels the 1/r^m singularity as r→0.

Tangent-point energy density (weight_functions.hpp):
  k_reg = |P dp|^p / (|dp|^{2p} + l0^p)          # current (denom regularized)
  k_raw = |P dp|^p /  |dp|^{2p}                   # singular radial factor 1/r^p

Beale–Majda-style on the raw radial singularity:
  k_bm  = |P dp|^p * (1 - exp(-(|dp|/β)^{2p})) / |dp|^{2p}
        = |N·û|^p  * (1 - exp(-(r/β)^{2p})) / r^p

As r→0: 1-exp(-(r/β)^{2p}) ~ (r/β)^{2p}, so k_bm ~ |N·û|^p * r^p / β^{2p} → 0.

Also compare soft-norm composition:
  f_eps(v) = sqrt(|v|^2 + eps^2)
  k_soft = f_eps(P dp)^p / (f_eps(dp)^{2p} + l0^p)

Run:
  python3 scripts/sympy/tangent_point_beale_majda.py
"""

from __future__ import annotations

import numpy as np
import sympy as sp


def symbols_setup():
    x, y, z = sp.symbols("x y z", real=True)
    nx, ny, nz = sp.symbols("nx ny nz", real=True)
    l0, beta, eps, p = sp.symbols("l0 beta eps p", positive=True)
    dp = sp.Matrix([x, y, z])
    N = sp.Matrix([nx, ny, nz])
    return dp, N, l0, beta, eps, p, (x, y, z)


def kernels(dp, N, l0, beta, eps, p):
    """Return (name -> scalar sympy expr) for several mollifications."""
    # |N|=1 ⇒ P dp = (N·dp) N, |P dp| = |N·dp|
    ndp = N.dot(dp)
    f = sp.Abs(ndp)
    g = sp.sqrt(dp.dot(dp))

    k_raw = f**p / g ** (2 * p)
    k_reg = f**p / (g ** (2 * p) + l0**p)

    # Beale–Majda on the 1/r^{2p} factor (same spirit as (1-e^{-(r/β)^3})/r^3).
    # Power 2p matches the singular radial order of the TP denominator.
    moll = 1 - sp.exp(-((g / beta) ** (2 * p)))
    k_bm = f**p * moll / g ** (2 * p)

    # Optional: BM *and* keep a tiny floor so far-field matches k_reg asymptotics.
    # (Not required; included for comparison.)
    k_bm_floor = f**p * moll / (g ** (2 * p) + l0**p)

    # Soft-norm composition (variational, chain-rule safe).
    f_s = sp.sqrt(ndp**2 + eps**2)
    g_s = sp.sqrt(dp.dot(dp) + eps**2)
    k_soft = f_s**p / (g_s ** (2 * p) + l0**p)

    # Existing calc_mollified-style on |dp|: (1-exp(-(r/l0)^p))/r^p times |N·û|^p
    # ⇒ (1-exp(-(r/l0)^p)) |N·û|^p / r^p
    moll_lib = 1 - sp.exp(-((g / l0) ** p))
    k_lib = (f / g) ** p * moll_lib / g**p  # = f^p * moll_lib / g^{2p}

    return {
        "raw": k_raw,
        "reg_l0": k_reg,
        "beale_majda": k_bm,
        "bm+l0_floor": k_bm_floor,
        "soft_norm": k_soft,
        "lib_mollified": k_lib,
    }, f, g, ndp


def grad_expr(k, vars_xyz):
    return sp.Matrix([sp.diff(k, v) for v in vars_xyz])


def lambdify_grad(k, dp, N, params, vars_xyz):
    """params: ordered symbols after (x,y,z,nx,ny,nz)."""
    x, y, z = vars_xyz
    nx, ny, nz = N[0], N[1], N[2]
    g = grad_expr(k, vars_xyz)
    args = (x, y, z, nx, ny, nz, *params)
    return sp.lambdify(args, (k, g), "numpy")


def sample_near_singularities(rng, n_each=40):
    """Cases: (a) r→0, (b) glancing f→0 with fixed r, (c) generic."""
    samples = []

    # (a) approach origin along / across N
    for _ in range(n_each):
        N = rng.normal(size=3)
        N /= np.linalg.norm(N)
        r = float(10 ** rng.uniform(-4, -1))
        # mixed direction
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        samples.append(("r_small", r * d, N))

    # (b) glancing: N·dp ≈ 0, |dp| ~ O(1)
    for _ in range(n_each):
        N = rng.normal(size=3)
        N /= np.linalg.norm(N)
        t = rng.normal(size=3)
        t -= t.dot(N) * N
        t /= np.linalg.norm(t)
        r = float(rng.uniform(0.2, 1.5))
        # tiny normal component
        f = float(10 ** rng.uniform(-4, -2))
        samples.append(("glancing", r * t + f * N, N))

    # (c) generic bulk
    for _ in range(n_each):
        N = rng.normal(size=3)
        N /= np.linalg.norm(N)
        dp = rng.normal(size=3)
        if np.linalg.norm(dp) < 1e-3:
            continue
        samples.append(("bulk", dp, N))

    return samples


def finite_diff_grad(fn, dp, N, extras, eps=1e-6):
    g = np.zeros(3)
    for ax in range(3):
        dpp, dpm = dp.copy(), dp.copy()
        dpp[ax] += eps
        dpm[ax] -= eps
        kp = float(np.asarray(fn(*dpp, *N, *extras)[0]).reshape(()))
        km = float(np.asarray(fn(*dpm, *N, *extras)[0]).reshape(()))
        g[ax] = (kp - km) / (2 * eps)
    return g


def main() -> None:
    dp, N, l0, beta, eps, p, vars_xyz = symbols_setup()
    ks, f_s, g_s, ndp = kernels(dp, N, l0, beta, eps, p)

    print("=== kernels (symbolic, N unit) ===")
    for name, k in ks.items():
        print(f"\n{name}:")
        # Specialize N=e_z for readability
        kk = sp.simplify(k.subs({N[0]: 0, N[1]: 0, N[2]: 1, p: 2}))
        sp.pprint(kk)

    # Limit r→0 along N for p=2, N=ez, dp=(0,0,r)
    print("\n=== limit r→0 along N (p=2, N=ez, dp=(0,0,r)) ===")
    r = sp.symbols("r", positive=True)
    subs_line = {
        N[0]: 0,
        N[1]: 0,
        N[2]: 1,
        dp[0]: 0,
        dp[1]: 0,
        dp[2]: r,
        p: 2,
    }
    for name, k in ks.items():
        expr = sp.simplify(k.subs(subs_line))
        lim = sp.limit(expr, r, 0)
        print(f"  {name:14s}  k→ {lim}   (expr = {expr})")

    # Numeric stress test near singularities
    rng = np.random.default_rng(0)
    samples = sample_near_singularities(rng)

    # Fixed scales (typical rod-demo-ish)
    l0_v, beta_v, eps_v, p_v = 0.1, 0.1, 0.1, 6.0

    # Build lambdas. Param lists differ slightly.
    fns = {}
    fns["reg_l0"] = lambdify_grad(ks["reg_l0"], dp, N, (l0, p), vars_xyz)
    fns["beale_majda"] = lambdify_grad(ks["beale_majda"], dp, N, (beta, p), vars_xyz)
    fns["bm+l0_floor"] = lambdify_grad(
        ks["bm+l0_floor"], dp, N, (beta, l0, p), vars_xyz
    )
    fns["soft_norm"] = lambdify_grad(ks["soft_norm"], dp, N, (eps, l0, p), vars_xyz)
    fns["lib_mollified"] = lambdify_grad(
        ks["lib_mollified"], dp, N, (l0, p), vars_xyz
    )

    extras = {
        "reg_l0": (l0_v, p_v),
        "beale_majda": (beta_v, p_v),
        "bm+l0_floor": (beta_v, l0_v, p_v),
        "soft_norm": (eps_v, l0_v, p_v),
        "lib_mollified": (l0_v, p_v),
    }

    # Also C++ production gradient of k_reg (ad-hoc +l2 form)
    Px = ndp * N
    k_reg = ks["reg_l0"]
    l2 = l0**2
    cpp = p * k_reg * (
        Px / (Px.dot(Px) + l2) - 2 * dp / (dp.dot(dp) + l2)
    )
    f_cpp = sp.lambdify(
        (vars_xyz[0], vars_xyz[1], vars_xyz[2], N[0], N[1], N[2], l0, p),
        cpp,
        "numpy",
    )

    print("\n=== near-singularity numeric probe ===")
    print(f"l0=beta=eps={l0_v}, p={p_v}")
    print(
        f"{'case':10s} {'kernel':14s} {'max|k|':>10s} {'max|grad|':>12s} "
        f"{'mean|grad|':>12s} {'FD relerr':>10s}"
    )

    cases = ("r_small", "glancing", "bulk")
    for case in cases:
        bucket = [(dp_, N_) for c, dp_, N_ in samples if c == case]
        for name, fn in fns.items():
            ks_abs = []
            gs_abs = []
            fd_errs = []
            for dpi, Ni in bucket:
                ex = extras[name]
                kv, gv = fn(*dpi, *Ni, *ex)
                kv = float(np.asarray(kv).reshape(()))
                gv = np.asarray(gv, dtype=float).reshape(3)
                if not np.all(np.isfinite(gv)) or not np.isfinite(kv):
                    ks_abs.append(np.inf)
                    gs_abs.append(np.inf)
                    continue
                ks_abs.append(abs(kv))
                gs_abs.append(np.linalg.norm(gv))
                # FD check a few
                if len(fd_errs) < 8:
                    gfd = finite_diff_grad(fn, dpi, Ni, ex)
                    denom = np.linalg.norm(gv) + 1e-14
                    fd_errs.append(np.linalg.norm(gfd - gv) / denom)

            print(
                f"{case:10s} {name:14s} {np.max(ks_abs):10.3e} "
                f"{np.max(gs_abs):12.3e} {np.mean(gs_abs):12.3e} "
                f"{np.mean(fd_errs) if fd_errs else float('nan'):10.2e}"
            )

        # C++ ad-hoc grad magnitude on same samples (of k_reg)
        gs_cpp = []
        for dpi, Ni in bucket:
            gv = np.asarray(f_cpp(*dpi, *Ni, l0_v, p_v), dtype=float).reshape(3)
            gs_cpp.append(np.linalg.norm(gv) if np.all(np.isfinite(gv)) else np.inf)
        print(
            f"{case:10s} {'cpp_l2_hack':14s} {'(n/a)':>10s} "
            f"{np.max(gs_cpp):12.3e} {np.mean(gs_cpp):12.3e} {'(n/a)':>10s}"
        )
        print()

    # Exact ∇ of beale_majda vs limit structure
    print("=== verdict notes ===")
    print(
        "- beale_majda: kills 1/r^{2p} like Biot–Savart mollify; k→0 and ∇k stays\n"
        "  finite as r→0 (variational: differentiate the mollified energy).\n"
        "- soft_norm: also variational; handles glancing f→0 via hypot(|N·dp|,eps).\n"
        "- reg_l0: energy OK; current C++ grad is NOT exact ∇ of reg_l0.\n"
        "- For TP contact (glancing), soft_norm / BM address different axes:\n"
        "  BM ↔ radial self-distance, soft_norm ↔ |N·dp| vanishing.\n"
        "- Best combo for TP forces: soft_norm (or BM) energy + exact autodiff grad,\n"
        "  optionally psi(k)=k/(1+a k) saturation on top."
    )


if __name__ == "__main__":
    main()
