#!/usr/bin/env python3
"""Composable soft-floor TP density: r → soft_floor → (1/rs)^p (+ full Hess).

Stages (N fixed, N·dp > 0 for FD):
  r  = tan_point(dp, N)           # |dp|² / (2 N·dp)
  rs = soft_floor(r, R_min, τ)    # softmax_τ(r, R_min)
  f  = density(rs, p)             # (1/rs)^p

Chain:
  ∇f = f'(r) ∇r
  H  = f''(r) (∇r)(∇r)ᵀ + f'(r) ∇²r     # full Hess — free from composition

Run:
  python3 scripts/sympy/tangent_point_soft.py
  python3 scripts/sympy/tangent_point_soft.py --dump
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import sympy as sp


# --- stage definitions (symbolic / numeric) ---------------------------------


def soft_floor(a: sp.Expr, b: sp.Expr, tau: sp.Expr) -> sp.Expr:
    """softmax_τ(a,b) ≈ max(a,b)."""
    return tau * sp.log(sp.exp(a / tau) + sp.exp(b / tau))


def density(rs: sp.Expr, p: sp.Expr) -> sp.Expr:
    return rs ** (-p)


def tan_point_sym(dp: sp.Matrix, N: sp.Matrix) -> sp.Expr:
    """Geometric TP radius (N·dp > 0 branch)."""
    return dp.dot(dp) / (2 * N.dot(dp))


def soft_floor_n(a: float, b: float, tau: float) -> float:
    M = max(a, b)
    return M + tau * np.log(np.exp((a - M) / tau) + np.exp((b - M) / tau))


def soft_floor_dr_n(a: float, b: float, tau: float) -> float:
    """∂softmax(a,b)/∂a."""
    M = max(a, b)
    ea = np.exp((a - M) / tau)
    eb = np.exp((b - M) / tau)
    return ea / (ea + eb)


def soft_floor_d2r_n(a: float, b: float, tau: float) -> float:
    """∂²softmax(a,b)/∂a² = w(1-w)/τ."""
    w = soft_floor_dr_n(a, b, tau)
    return w * (1.0 - w) / tau


def tan_point_n(dp: np.ndarray, N: np.ndarray) -> float:
    return float(np.dot(dp, dp) / (2.0 * np.dot(N, dp)))


def density_n(rs: float, p: float) -> float:
    return (1.0 / rs) ** p


def density_drs_n(rs: float, p: float) -> float:
    return -p * density_n(rs, p) / rs


def density_d2rs_n(rs: float, p: float) -> float:
    return p * (p + 1.0) * density_n(rs, p) / (rs * rs)


def grad_tan_point_n(dp: np.ndarray, N: np.ndarray) -> np.ndarray:
    fs = float(np.dot(N, dp))
    g2 = float(np.dot(dp, dp))
    return (2.0 * fs * dp - g2 * N) / (2.0 * fs * fs)


def hess_tan_point_n(dp: np.ndarray, N: np.ndarray) -> np.ndarray:
    """Full ∇²r for r = |dp|²/(2 N·dp), N·dp > 0."""
    fs = float(np.dot(N, dp))
    g2 = float(np.dot(dp, dp))
    # From sympy hessian of g2/(2*fs):
    # H = I/fs - (dp⊗N + N⊗dp)/fs² + (g2 N⊗N)/fs³
    I = np.eye(3)
    dpN = np.outer(dp, N)
    Ndp = np.outer(N, dp)
    NN = np.outer(N, N)
    return I / fs - (dpN + Ndp) / (fs * fs) + (g2 * NN) / (fs ** 3)


def compose_grad_hess_n(dp: np.ndarray, N: np.ndarray, R_min: float, tau: float,
                        p: float):
    """f = density(soft_floor(tan_point)). Returns f, ∇f, H."""
    r = tan_point_n(dp, N)
    rs = soft_floor_n(r, R_min, tau)
    f = density_n(rs, p)

    w = soft_floor_dr_n(r, R_min, tau)
    w2 = soft_floor_d2r_n(r, R_min, tau)
    fp = density_drs_n(rs, p)   # df/drs
    fpp = density_d2rs_n(rs, p)  # d²f/drs²

    # df/dr, d²f/dr² through soft_floor
    d1 = fp * w
    d2 = fpp * w * w + fp * w2

    gr = grad_tan_point_n(dp, N)
    Hr = hess_tan_point_n(dp, N)
    gk = d1 * gr
    Hk = d2 * np.outer(gr, gr) + d1 * Hr
    return f, gk, Hk


# --- sympy: confirm Hess_r formula ------------------------------------------


def symbolic_hess_r() -> None:
    x, y, z = sp.symbols("x y z")
    nx, ny, nz = sp.symbols("nx ny nz")
    dp = sp.Matrix([x, y, z])
    N = sp.Matrix([nx, ny, nz])
    r = tan_point_sym(dp, N)
    H = sp.simplify(sp.hessian(r, [x, y, z]))
    fs = N.dot(dp)
    g2 = dp.dot(dp)
    H_closed = sp.simplify(
        sp.eye(3) / fs
        - (dp * N.T + N * dp.T) / fs ** 2
        + (g2 * (N * N.T)) / fs ** 3
    )
    print("=== Hess_r closed form vs sympy.hessian ===")
    print("residual:", sp.simplify(H - H_closed))


def symbolic_chain() -> None:
    R, Rmin, tau, p = sp.symbols("R R_min tau p", positive=True)
    rs_sym = sp.symbols("rs", positive=True)
    rs = soft_floor(R, Rmin, tau)
    f = density(rs_sym, p)
    print("=== composable stages ===")
    print("rs = soft_floor(R) = softmax(R, R_min)")
    print("f  = density(rs)   = rs^{-p}")
    print("∂rs/∂R =")
    sp.pprint(sp.simplify(sp.diff(rs, R)))
    print("∂f/∂rs =")
    sp.pprint(sp.simplify(sp.diff(f, rs_sym)))
    print("\nFull Hess_dp(f) = f''(∇r)(∇r)ᵀ + f' ∇²r  (free from composition)")


# --- FD --------------------------------------------------------------------


def numeric_fd_check(n_samples: int = 40, seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    R_min = 0.05
    tau = 0.005
    p = 4.0
    g_errs = []
    H_errs = []

    for _ in range(n_samples * 10):
        if len(g_errs) >= n_samples:
            break
        dp = rng.normal(size=3)
        N = rng.normal(size=3)
        N /= np.linalg.norm(N)
        if float(np.dot(N, dp)) < 1e-2:
            continue
        if float(np.dot(dp, dp)) < 1e-4:
            continue
        r = tan_point_n(dp, N)
        if not (R_min * 1.5 < r < R_min * 20.0):
            continue

        f0, gk, Hk = compose_grad_hess_n(dp, N, R_min, tau, p)

        eps = 1e-6

        def f_of(d):
            if float(np.dot(N, d)) <= 0.0:
                return np.nan
            return compose_grad_hess_n(d, N, R_min, tau, p)[0]

        def g_of(d):
            if float(np.dot(N, d)) <= 0.0:
                return None
            return compose_grad_hess_n(d, N, R_min, tau, p)[1]

        gfd = np.zeros(3)
        Hfd = np.zeros((3, 3))
        ok = True
        for a in range(3):
            dpp, dpm = dp.copy(), dp.copy()
            dpp[a] += eps
            dpm[a] -= eps
            fp, fm = f_of(dpp), f_of(dpm)
            gp, gm = g_of(dpp), g_of(dpm)
            if not (np.isfinite(fp) and np.isfinite(fm)) or gp is None or gm is None:
                ok = False
                break
            gfd[a] = (fp - fm) / (2 * eps)
            Hfd[:, a] = (gp - gm) / (2 * eps)
        if not ok:
            continue

        g_errs.append(np.linalg.norm(gk - gfd) / (np.linalg.norm(gfd) + 1e-14))
        H_errs.append(np.linalg.norm(Hk - Hfd) / (np.linalg.norm(Hfd) + 1e-14))

    print("=== composable soft TP FD ===")
    print("r=tan_point → rs=soft_floor → f=(1/rs)^p")
    print(
        f"grad: samples={len(g_errs)}  rel err mean={np.mean(g_errs):.3e}  max={np.max(g_errs):.3e}"
    )
    print(
        f"Hess: samples={len(H_errs)}  rel err mean={np.mean(H_errs):.3e}  max={np.max(H_errs):.3e}"
    )
    if not g_errs or np.max(g_errs) > 1e-4 or np.max(H_errs) > 1e-3:
        print("WARNING: FD residual high", file=sys.stderr)


# --- C++ dump --------------------------------------------------------------


def dump_cpp() -> str:
    """Pasteable C++: stage helpers + composed grad/Hess."""
    return r'''
    // Stages: r = tan_point(dp,N); rs = soft_floor(r); f = density(rs)
    // (from scripts/sympy/tangent_point_soft.py)

    /// Geometric TP radius R = |dp|² / (2 |N·dp|).
    inline real tan_point_radius(const vec3 &dp, const vec3 &N) {
      const real f = std::abs(N.dot(dp));
      const real g2 = dp.squaredNorm();
      if (!(f > real(1.0e-16)) || !(g2 > real(1.0e-32)))
        return real(0.0);
      return real(0.5) * g2 / f;
    }

    /// Soft contact floor: rs = softmax_τ(r, R_min) ≈ max(r, R_min).
    inline real soft_floor_radius(real r, real R_min, real tau) {
      return tp_softmax(r, R_min, tau);
    }

    /// Classic density f = (1/rs)^p.
    inline real soft_tp_density(real rs, real p) {
      return std::pow(real(1.0) / rs, p);
    }

    /// ∂softmax(a,b)/∂a (stable).
    inline real soft_floor_dr(real a, real b, real tau) {
      const real inv_tau = real(1.0) / tau;
      const real M = std::max(a, b);
      const real ea = std::exp((a - M) * inv_tau);
      const real eb = std::exp((b - M) * inv_tau);
      return ea / (ea + eb);
    }

    /// ∂²softmax(a,b)/∂a² = w(1-w)/τ.
    inline real soft_floor_d2r(real a, real b, real tau) {
      const real w = soft_floor_dr(a, b, tau);
      return w * (real(1.0) - w) / tau;
    }

    /// ∇_dp R for R = |dp|²/(2|N·dp|) (N fixed).
    inline vec3 tan_point_radius_grad(const vec3 &dp, const vec3 &N) {
      const real f_s = N.dot(dp);
      const real f = std::abs(f_s);
      const real g2 = dp.squaredNorm();
      if (!(f > real(1.0e-16)) || !(g2 > real(1.0e-32)))
        return vec3::Zero();
      const real sgn = (f_s >= real(0.0)) ? real(1.0) : real(-1.0);
      return (real(2.0) * f * dp - g2 * (sgn * N)) / (real(2.0) * f * f);
    }

    /// ∇²_dp R for R = |dp|²/(2 N·dp) on the N·dp > 0 chart:
    ///   H = I/f_s − (dp⊗N + N⊗dp)/f_s² + (|dp|² N⊗N)/f_s³
    /// For N·dp < 0, flip N → −N (same as abs-branch chart).
    inline mat3 tan_point_radius_hess(const vec3 &dp, const vec3 &N) {
      const real f_s = N.dot(dp);
      const real g2 = dp.squaredNorm();
      if (!(std::abs(f_s) > real(1.0e-16)) || !(g2 > real(1.0e-32)))
        return mat3::Zero();
      const vec3 Ns = (f_s >= real(0.0)) ? N : vec3(-N);
      const real fs = std::abs(f_s);
      mat3 H = mat3::Identity() / fs;
      H.noalias() -= (dp * Ns.transpose() + Ns * dp.transpose()) / (fs * fs);
      H.noalias() += (g2 * (Ns * Ns.transpose())) / (fs * fs * fs);
      return H;
    }

    /// Softmax-floor TP density f = density(soft_floor(tan_point)).
    inline real calc_tangent_point_radius_soft(const vec3 &dp, const vec3 &N,
                                              real R_min, real tau, real p) {
      if (!(R_min > 0.0) || !(tau > 0.0))
        return real(0.0);
      const real r = tan_point_radius(dp, N);
      if (!(r > 0.0))
        return real(0.0);
      const real rs = soft_floor_radius(r, R_min, tau);
      return soft_tp_density(rs, p);
    }

    /// ∇_dp f (N fixed). ∇f = f'(r) ∇r.
    inline vec3 calc_tangent_point_radius_gradient_soft(const vec3 &dp,
                                                       const vec3 &N,
                                                       real R_min, real tau,
                                                       real p) {
      if (!(R_min > 0.0) || !(tau > 0.0))
        return vec3::Zero();
      const real r = tan_point_radius(dp, N);
      if (!(r > 0.0))
        return vec3::Zero();
      const real rs = soft_floor_radius(r, R_min, tau);
      const real f = soft_tp_density(rs, p);
      const real w = soft_floor_dr(r, R_min, tau);
      const real df_dr = (-p * f / rs) * w;
      const vec3 g = df_dr * tan_point_radius_grad(dp, N);
      if (!std::isfinite(g[0]) || !std::isfinite(g[1]) || !std::isfinite(g[2]))
        return vec3::Zero();
      return g;
    }

    /// Full Hess_dp f (N fixed):
    ///   H = f''(r) (∇r)(∇r)ᵀ + f'(r) ∇²r
    inline mat3 calc_tangent_point_radius_hessian_soft(const vec3 &dp,
                                                      const vec3 &N,
                                                      real R_min, real tau,
                                                      real p) {
      if (!(R_min > 0.0) || !(tau > 0.0))
        return mat3::Zero();
      const real r = tan_point_radius(dp, N);
      if (!(r > 0.0))
        return mat3::Zero();
      const real rs = soft_floor_radius(r, R_min, tau);
      const real f = soft_tp_density(rs, p);
      const real w = soft_floor_dr(r, R_min, tau);
      const real w2 = soft_floor_d2r(r, R_min, tau);
      const real df_drs = -p * f / rs;
      const real d2f_drs2 = p * (p + real(1.0)) * f / (rs * rs);
      const real df_dr = df_drs * w;
      const real d2f_dr2 = d2f_drs2 * w * w + df_drs * w2;
      const vec3 gr = tan_point_radius_grad(dp, N);
      const mat3 Hr = tan_point_radius_hess(dp, N);
      mat3 H = d2f_dr2 * (gr * gr.transpose()) + df_dr * Hr;
      if (!std::isfinite(H(0, 0)) || !std::isfinite(H(1, 1)) ||
          !std::isfinite(H(2, 2)))
        return mat3::Zero();
      return H;
    }
'''


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", action="store_true", help="print pasteable C++ only")
    args = ap.parse_args()

    if args.dump:
        print(dump_cpp())
        return

    symbolic_hess_r()
    print()
    symbolic_chain()
    print()
    numeric_fd_check()
    print()
    print("=== pasteable C++ (also: --dump) ===")
    print(dump_cpp())


if __name__ == "__main__":
    main()
