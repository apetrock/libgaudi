#!/usr/bin/env python3
"""Where Beale–Majda mollification comes from, and how to build a TP kernel.

Biot–Savart / Newtonian singular kernel in 3D (up to constants):
  K(x) = x / |x|^3          # ~ grad(1/r)

Vortex-blob idea (Beale–Majda, Chorin, …):
  replace the singular measure by a smooth approximate identity ζ_ε,
  and DEFINE the mollified kernel by CONVOLUTION:

    K_ε := K * ζ_ε
         = ∫ K(x - y) ζ_ε(y) dy

For radial blobs ζ_ε(x) = ε^{-3} ζ(|x|/ε) with ∫ ζ = 1, symmetry gives
  K_ε(x) = K(x) * m(|x|/ε)
i.e. the singular direction x/|x|^3 times a radial cutoff m(ρ).

Brochu SCA2012 quote one closed form (Beale–Majda style):
  K_ε(x) = x * (1 - exp(-(r/β)^3)) / r^3
so
  m(ρ) = 1 - exp(-ρ^3)     with ρ = r/β

That m is NOT an arbitrary f(g) guess — it is the radial factor of a
particular blob convolution. As ρ→0, m(ρ)~ρ^3 ⇒ K_ε ~ x/β^3 (bounded).
As ρ→∞, m→1 ⇒ K_ε → K.

Tangent-point energy density is ALSO a singular radial kernel:
  k_raw(dp, N) = |N·dp|^p / |dp|^{2p}
               = |N·û|^p  /  r^p

So the same two constructions apply:

  (A) CONVOLUTION  h = k_sing * ζ_ε
      mollify the singular radial profile 1/r^p (or the full kernel in R^3)
      against a blob → h(dp) = |N·û|^p / r^p * m(r/ε)
      with m from the blob (BM is one choice of m).

  (B) COMPOSITION  h = f(g)
      pick a geometric scalar g (soft TP radius, soft |N·dp|, …)
      and a saturating f (log, k/(1+αk), 1-e^{-k}, …).
      Variationally clean, but NOT automatically a blob convolution
      unless f∘g is chosen to match some K*ζ.

This script:
  1) shows the radial-factor structure of BM
  2) builds TP candidates of type (A) and (B)
  3) checks r→0 limits and compares peakiness

Run:
  python3 scripts/sympy/beale_majda_derivation.py
"""

from __future__ import annotations

import numpy as np
import sympy as sp


def radial_bm_factor():
    """K_ε = x/r^3 * m(r/β) with m(ρ)=1-exp(-ρ^3)."""
    r, beta = sp.symbols("r beta", positive=True)
    m = 1 - sp.exp(-((r / beta) ** 3))
    K_sing_mag = 1 / r**2  # |x/r^3| = 1/r^2
    K_moll_mag = m / r**2
    print("=== Beale–Majda as radial factor of a convolution ===")
    print("singular |K|     = 1/r^2")
    print("mollified |K_ε|  = (1 - exp(-(r/β)^3)) / r^2")
    print("limit r→0 |K_ε|  =", sp.limit(K_moll_mag, r, 0))
    print("limit r→∞ |K_ε|/|K| =", sp.limit(m, r, sp.oo))
    print()
    # series: m/r^2 = (r/β)^3/r^2 + O(r^4) = r/β^3 + …
    print("series r→0 of |K_ε|:")
    sp.pprint(K_moll_mag.series(r, 0, 4).removeO())
    print()


def tp_convolution_family():
    """Type (A): singular TP radial profile times blob factor m(r/ε)."""
    r, eps, p = sp.symbols("r eps p", positive=True)
    # angular factor |N·û|^p treated as O(1) multiplier A in [0,1]
    A = sp.symbols("A", nonnegative=True)

    k_sing = A / r**p

    # Several standard radial mollifier factors (all m(0)=0 behavior or bounded product)
    mollifiers = {
        # Beale–Majda 3D-ish: m = 1-exp(-ρ^3); product ~ A r^{3-p}/ε^3 near 0
        "bm3": (1 - sp.exp(-((r / eps) ** 3))) / r**p,
        # Match singular order: m = 1-exp(-ρ^p)  (your calc_mollified)
        "bm_p": (1 - sp.exp(-((r / eps) ** p))) / r**p,
        # Algebraic blob (Cauchy / Rosenhead): 1/(r^2+ε^2)^{p/2}
        "cauchy": 1 / (r**2 + eps**2) ** (p / 2),
        # Gaussian approximate-identity convolution surrogate (not exact K*ζ)
        "gauss_cutoff": sp.exp(-((r / eps) ** 2)) * 0 + (1 - sp.exp(-((r / eps) ** 2))) / r**p,
    }

    print("=== (A) convolution-style TP: A * (singular radial * blob) ===")
    print("k_sing = A / r^p   (A = |N·û|^p)")
    print()
    for name, radial in mollifiers.items():
        k = sp.simplify(A * radial)
        lim = sp.limit(k.subs(p, 6), r, 0)
        print(f"{name:14s}  lim_r→0 (p=6) = {lim}")
        # peak scale rough: evaluate on a grid numerically below
    print()
    return mollifiers, A, r, eps, p


def tp_composition_family():
    """Type (B): h = f(g) with geometric g."""
    # g1 = classical TP inverse radius ρ = |N·dp| / |dp|^2 = A^{1/p} / r
    # use scalars f = |N·dp|, g = |dp|
    f, g, eps, l0, p, alpha = sp.symbols("f g eps l0 p alpha", positive=True)

    rho_raw = f / g**2
    rho_soft = f / (g**2 + eps**2)  # Rosenhead-style on denominator
    rho_hypot = sp.sqrt(f**2 + eps**2) / (g**2 + eps**2)

    compositions = {
        "power_reg": (f**p) / (g ** (2 * p) + l0**p),  # current
        "f_of_rho_soft": (rho_soft) ** p,  # f(g)=g^p with g=soft ρ
        "f_of_rho_hypot": (rho_hypot) ** p,
        "saturate_reg": (f**p) / (g ** (2 * p) + l0**p)
        / (1 + alpha * (f**p) / (g ** (2 * p) + l0**p)),  # k/(1+αk)
        "log_reg": sp.log(1 + (f**p) / (g ** (2 * p) + l0**p)),
    }

    print("=== (B) composition-style TP: h = f(g) ===")
    for name, expr in compositions.items():
        # approach along normal: f=g=r → 0
        r = sp.symbols("r", positive=True)
        lim = sp.limit(
            expr.subs({f: r, g: r, p: 6, l0: 1, eps: 1, alpha: 1}), r, 0
        )
        print(f"{name:16s}  along-N lim r→0 (p=6) = {lim}")
    print()
    return compositions


def numeric_peakiness():
    """Where does each radial profile peak? (A = 1)"""
    rs = np.logspace(-4, 1, 400)
    eps = 0.1
    p = 6.0
    A = 1.0

    def peak(name, vals):
        i = int(np.nanargmax(vals))
        print(
            f"  {name:16s}  peak at r={rs[i]:.4g}  value={vals[i]:.4g}  "
            f"at r=eps: {vals[np.argmin(np.abs(rs - eps))]:.4g}"
        )

    print("=== numeric peakiness of radial profiles (A=1, eps=0.1, p=6) ===")
    peak("raw 1/r^p", A / rs**p)
    peak("reg l0=eps", A / (rs ** (2 * p) + eps**p) * rs**p)  # wait: |N·dp|^p/(r^{2p}+l0^p) along N is r^p/(r^{2p}+l0^p)
    # along N: A=1 means |N·û|=1, k_reg = r^p / (r^{2p} + l0^p)
    peak("reg along-N", rs**p / (rs ** (2 * p) + eps**p))
    peak("bm3 / r^p", (1 - np.exp(-((rs / eps) ** 3))) / rs**p)
    peak("bm_p / r^p", (1 - np.exp(-((rs / eps) ** p))) / rs**p)
    peak("cauchy", 1.0 / (rs**2 + eps**2) ** (p / 2))
    # convolution-faithful TP: A/r^p * m_bm3, along N ⇒
    peak("TP*bm3 along-N", (1 - np.exp(-((rs / eps) ** 3))) / rs**p)
    peak(
        "TP*bm3 * l0floor",
        rs**p
        * (1 - np.exp(-((rs / eps) ** 3)))
        / (rs ** (2 * p) + eps**p),
    )
    print()


def design_recipe():
    print("=== design recipe for a new TP kernel ===")
    print(
        """
Convolution route (h = K * ζ) — closest to BI literature:
  1. Identify singular kernel K(dp,N) = |N·û|^p / r^p
  2. Pick radial blob ζ_ε (Gaussian, BM, algebraic, …)
  3. Set K_ε = K * ζ_ε  ⇒  for radial blobs, K_ε = |N·û|^p / r^p * m(r/ε)
  4. Energy density = K_ε (optionally + angular gates)
  5. Force = exact ∇_dp of that scalar (chain rule on m and û)

  BM choice m(ρ)=1-exp(-ρ^3) is one blob; your calc_mollified is
  m(ρ)=1-exp(-ρ^p) matched to order p. Algebraic m from Cauchy blobs
  is often smoother for large p (less shell peaking).

Composition route (h = f(g)) — energy shaping:
  1. Pick geometric g: soft TP radius ρ_ε = |N·dp| / (|dp|^2+ε^2)
                     or hypot forms, or current k_reg
  2. Pick f that is C^∞, f(0)=0, f'≥0, f'→0 at ∞ if you want saturation
  3. h = f(g); force = f'(g) ∇g   ← singularity dies if f' decays fast enough

  These coincide when f(g) equals some K*ζ closed form.
  Otherwise composition is a different (still valid) regularizer.

Practical compose for libgaudi:
  g = |N·dp|^p / (|dp|^{2p} + l0^p)           # current, or soft-norm variant
  h = g / (1 + α g)                           # saturate
  or
  h = |N·û|^p * (1 - exp(-(r/β)^q)) / (r^p + l0^p)   # BM radial * floor
"""
    )


def main() -> None:
    radial_bm_factor()
    tp_convolution_family()
    tp_composition_family()
    numeric_peakiness()
    design_recipe()


if __name__ == "__main__":
    main()
