#!/usr/bin/env python3
"""Filament-style TP mollification for rods — how close can SymPy get?

Rod kernel (T fixed; N = dp_⊥/|dp_⊥| from your integrator):
  K(dp; T) = |dp_⊥|^p / |dp|^{2p}
           = ρ^p / (ρ² + s²)^p
where s = dp·T (along-tangent separation), ρ = |dp_⊥|.

Vortex-tube style mollify: hold T (and s) fixed, convolve in the normal plane
with a 2D blob ζ_ε(q), q ∈ R²:

  K_ε(ρ, s) = ∫_{R²}  |q|^p / (|q|^2 + s²)^p   ζ_ε(ρ_vec - q) dq
            or evaluate at offset ρ with smear on the source cross-section.

Often collapses to a radial function of (ρ, s, ε).

Also compare closed algebraic surrogates (no integral):
  K_alg = ρ^p / (ρ² + s² + ε²)^p          # Rosenhead / Cauchy-style
  K_reg = ρ^p / ((ρ² + s²)^p + l0^p)     # current-style floor

Run:
  python3 scripts/sympy/tp_rod_mollify.py
"""

from __future__ import annotations

import numpy as np
import sympy as sp


def algebraic_closed_forms():
    rho, s, eps, l0, p = sp.symbols("rho s eps l0 p", positive=True)
    K_raw = rho**p / (rho**2 + s**2) ** p
    K_alg = rho**p / (rho**2 + s**2 + eps**2) ** p
    K_reg = rho**p / ((rho**2 + s**2) ** p + l0**p)

    print("=== algebraic surrogates (closed form, no integral) ===")
    for name, K in [("raw", K_raw), ("Rosenhead", K_alg), ("l0-floor", K_reg)]:
        print(f"\n{name}:")
        sp.pprint(K)

    # Gradient w.r.t. dp in basis: e_ρ, e_s (T)
    # K=K(ρ,s) ⇒ ∇_⊥ K = (∂K/∂ρ) e_ρ,  ∂_s K along T
    print("\n=== exact ∇ of Rosenhead K_alg (p symbolic) ===")
    d_rho = sp.simplify(sp.diff(K_alg, rho))
    d_s = sp.simplify(sp.diff(K_alg, s))
    print("∂K/∂ρ =")
    sp.pprint(d_rho)
    print("∂K/∂s =")
    sp.pprint(d_s)

    # Specialize p=6
    print("\n=== Rosenhead p=6, limits ===")
    K6 = K_alg.subs(p, 6)
    print("ρ→0, s fixed:", sp.limit(K6, rho, 0))
    print("s→0, ρ fixed:", sp.limit(K6, s, 0))
    print("ρ=s=0 along ρ=s:", sp.limit(K6.subs(s, rho), rho, 0))
    print()
    return K_alg, d_rho, d_s, rho, s, eps, p


def attempt_2d_gaussian_convolution():
    """K_ε(ρ,s) = ∫ K(|q|,s) ζ_ε(ρ_vec - q) d²q with ζ 2D Gaussian.

    At ρ=0 (on-axis): most symmetric — try SymPy polar integral.
    """
    print("=== 2D Gaussian convolution on-axis (ρ=0) ===")
    q, s, eps, p = sp.symbols("q s eps p", positive=True)
    # On-axis: by rotation symmetry
    # K_ε(0,s) = ∫_0^∞ [q^p / (q²+s²)^p] * (1/(2πε²)) exp(-q²/(2ε²)) * 2π q dq
    #          = ∫_0^∞ q^{p+1} / (q²+s²)^p * (1/ε²) exp(-q²/(2ε²)) dq
    # wait: 2D gaussian ζ = 1/(2π ε²) exp(-r²/(2ε²)), ∫ ζ 2π q dq = 1
    # K_ε(0,s) = ∫_0^∞ K(q,s) ζ(q) 2π q dq
    zeta_ring = (1 / eps**2) * sp.exp(-(q**2) / (2 * eps**2)) * q  # 2π*(1/(2π ε²))=1/ε²
    # Actually: ∫_0^∞ (1/(2π ε²)) exp * 2π q dq = ∫ (1/ε²) q exp dq = 1. Good.
    integrand = (q**p / (q**2 + s**2) ** p) * (1 / eps**2) * sp.exp(
        -(q**2) / (2 * eps**2)
    ) * q

    print("integrand for K_ε(0,s):")
    sp.pprint(integrand)
    print()

    # Try p=2 first (friendliest)
    print("SymPy integrate p=2, q=0..∞ ...")
    I2 = integrand.subs(p, 2)
    try:
        # Use a time-box via manual heuristics: expand or meijerg
        ans = sp.integrate(I2, (q, 0, sp.oo))
        ans = sp.simplify(ans)
        print("closed form p=2:")
        sp.pprint(ans)
    except Exception as e:
        print(f"failed: {e}")

    # Try with concrete eps,s via meijerg / numerical check path
    print()
    print("Trying p=2 with meijerg=True ...")
    try:
        ans = sp.integrate(I2, (q, 0, sp.oo), meijerg=True)
        print(sp.simplify(ans))
    except Exception as e:
        print(f"meijerg failed: {e}")
    print()


def numeric_on_axis_and_fit():
    """Numeric K_ε(0,s) vs algebraic surrogates — are we close?"""
    from numpy.polynomial.hermite import hermgauss

    eps = 0.1
    p = 6.0
    # Gauss-Laguerre / plain quadrature on q∈[0,∞)
    # Use q = eps * sqrt(2) * u with Hermite-like: ∫_0^∞ f(q) e^{-q²/(2ε²)} dq
    # Let t = q/(ε√2), dq = ε√2 dt, integrand has e^{-t²}
    # K_ε(0,s) = ∫_0^∞ [q^p/(q²+s²)^p] (1/ε²) e^{-q²/(2ε²)} q dq

    def Ke0(s, n=80):
        # Gauss-Hermite on (-∞,∞) but integrand even in disguise — use t≥0
        # Simple: fixed quad with scipy if available, else trapezoid log-grid
        qs = np.logspace(-6, np.log10(12 * eps), 400)
        # add dense near 0
        w = (1 / eps**2) * np.exp(-(qs**2) / (2 * eps**2)) * qs
        K = (qs**p) / ((qs**2 + s**2) ** p)
        # trapz with d(log q) carefully — use linear trapz on qs
        return float(np.trapz(K * w, qs))

    print("=== numeric on-axis K_ε(0,s) vs surrogates (p=6, ε=0.1) ===")
    print(f"{'s/ε':>8s} {'K_ε':>12s} {'Rosenhead':>12s} {'l0-floor':>12s} "
          f"{'raw':>12s}")
    for scale in [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]:
        s = scale * eps
        # at ρ=0, raw = 0 for p>0; Rosenhead ρ=0 ⇒ 0; on-axis convolution smears ρ
        # Wait: K(q,s)=q^p/(q²+s²)^p, at evaluation ρ=0 we integrate over q≠0
        # so K_ε(0,s) > 0 even though K(0,s)=0 for p>0.
        Ke = Ke0(s) if True else 0.0
        # For comparison at ρ=0: Rosenhead=0, l0-floor=0, raw=0
        # Better compare at ρ=ε (typical core radius offset)
        rho = eps
        K_ros = (rho**p) / ((rho**2 + s**2 + eps**2) ** p)
        K_flr = (rho**p) / (((rho**2 + s**2) ** p) + eps**p)
        K_raw = (rho**p) / ((rho**2 + s**2) ** p) if (rho + s) > 0 else np.inf
        # Also report Ke at ρ=0 (smear)
        print(
            f"{scale:8.2f} {Ke:12.4e} {K_ros:12.4e} {K_flr:12.4e} {K_raw:12.4e}"
        )
    print("(K_ε column = on-axis ρ=0 convolution; surrogates at ρ=ε)")
    print()

    # Direct compare: evaluate numeric convolution at ρ=ε vs Rosenhead(ρ=ε)
    print("=== numeric K_ε(ρ=ε, s) vs Rosenhead(ρ=ε,s) ===")

    def Ke_rho(rho, s, nphi=32, nq=200):
        # Monte-Carlo / polar quad: q in R², source at 0, eval at (ρ,0)
        # K_ε(ρ) = ∫ K(|q|,s) ζ(ρ_vec - q) d²q
        # ζ(u) = 1/(2π ε²) exp(-|u|²/(2ε²))
        qs = np.linspace(0.0, 10 * eps, nq)
        phis = np.linspace(0.0, 2 * np.pi, nphi, endpoint=False)
        dphi = 2 * np.pi / nphi
        # dq area ~ q dq dphi; use trap on q
        acc = 0.0
        for phi in phis:
            # eval point e = (ρ,0); source q = q(cos,sin); u = e - q
            qx = qs * np.cos(phi)
            qy = qs * np.sin(phi)
            ux = rho - qx
            uy = 0.0 - qy
            u2 = ux**2 + uy**2
            zeta = (1.0 / (2 * np.pi * eps**2)) * np.exp(-u2 / (2 * eps**2))
            Kq = (qs**p) / ((qs**2 + s**2) ** p)
            acc += np.trapz(Kq * zeta * qs, qs) * dphi
        return float(acc)

    print(f"{'s/ε':>8s} {'K_ε(ρ=ε)':>12s} {'Rosenhead':>12s} {'ratio':>8s}")
    for scale in [0.0, 0.5, 1.0, 2.0, 4.0]:
        s = scale * eps
        rho = eps
        Ke = Ke_rho(rho, s)
        Kr = (rho**p) / ((rho**2 + s**2 + eps**2) ** p)
        ratio = Ke / Kr if Kr > 0 else np.nan
        print(f"{scale:8.2f} {Ke:12.4e} {Kr:12.4e} {ratio:8.3f}")
    print()


def cpp_ready_gradient():
    """Emit a compact C++-ish formula for Rosenhead rod TP grad."""
    print("=== practical close form for C++ (Rosenhead filament TP) ===")
    print(
        """
With T unit, dp, s = dp·T, dp_⊥ = dp - s T, ρ² = |dp_⊥|²,
  K = ρ^p / (ρ² + s² + ε²)^p

Let R2ε = ρ² + s² + ε² = |dp|² + ε².
  K = |dp_⊥|^p / R2ε^p

∇_dp K (exact) can be coded from ∂/∂ρ, ∂/∂s:
  ∂K/∂ρ = p ρ^{p-1} / R2ε^p  -  p ρ^p * 2ρ / R2ε^{p+1}
        = p ρ^{p-1} ( R2ε - 2 ρ² ) / R2ε^{p+1}
  ∂K/∂s = - p ρ^p * 2s / R2ε^{p+1}

  ∇K = (∂K/∂ρ) e_ρ + (∂K/∂s) T
     with e_ρ = dp_⊥/ρ  (ρ>0; else 0)

This is "close" to true (K*ζ_2D): same symmetry, same ε as core radius,
closed form, variationally exact for this energy. SymPy derives ∇ easily;
the integral K*ζ is what we compare numerically above.
"""
    )


def main() -> None:
    algebraic_closed_forms()
    # Keep symbolic integrate short — may be slow; skip if hangs by using timeout externally
    try:
        attempt_2d_gaussian_convolution()
    except Exception as e:
        print("symbolic convolution aborted:", e)
    numeric_on_axis_and_fit()
    cpp_ready_gradient()


if __name__ == "__main__":
    main()
