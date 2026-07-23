"""Implicit surface shape operators — shared by quadric and Darboux cyclide."""

from __future__ import annotations

from typing import Dict, Tuple

import sympy as sp
from sympy import Matrix, Symbol, symbols
from sympy.core.expr import Expr

x, y, z = symbols("x y z", real=True)
t = sp.Symbol("t", real=True)


def gradient(D: Expr) -> Matrix:
    return Matrix([sp.diff(D, x), sp.diff(D, y), sp.diff(D, z)])


def hessian(D: Expr) -> Matrix:
    g = gradient(D)
    return Matrix([[sp.diff(g[i], v) for v in (x, y, z)] for i in range(3)])


def unit_normal(D: Expr) -> Matrix:
    g = gradient(D)
    return g / sp.sqrt(g.dot(g))


def shape_operator(D: Expr) -> Matrix:
    """Weingarten map W = (1/|∇D|) P H P."""
    g = gradient(D)
    H = hessian(D)
    n = unit_normal(D)
    I3 = sp.eye(3)
    P = I3 - n * n.T
    return (P * H * P) / sp.sqrt(g.dot(g))


def principal_curvatures_at(D: Expr, point: Dict[Symbol, float]) -> Tuple[float, float]:
    """Numeric κ₁, κ₂ at a point (eigenvalues of W; drop the ~0 eigenvalue)."""
    import numpy as np

    W = shape_operator(D).subs(point)
    Wn = np.array([[complex(sp.N(W[i, j])).real for j in range(3)] for i in range(3)])
    vals = np.linalg.eigvalsh(Wn)
    vals = sorted(vals, key=lambda v: abs(v), reverse=True)
    return float(vals[0]), float(vals[1])


def principal_curvatures_symbolic(D: Expr) -> Tuple[Expr, Expr]:
    """Symbolic κ₁, κ₂ — only use for simple D (e.g. sphere). Slow on general quadrics."""
    W = shape_operator(D)
    lam = sp.Symbol("lam")
    poly = sp.Poly(sp.expand(W.charpoly(lam).as_expr()), lam)
    roots = list(sp.roots(poly, lam).keys())
    roots_sorted = sorted(roots, key=lambda r: abs(complex(sp.N(r.subs({}), 6))), reverse=True)
    if len(roots_sorted) < 2:
        return sp.S.Zero, sp.S.Zero
    return sp.simplify(roots_sorted[0]), sp.simplify(roots_sorted[1])


def march_position(x0: Matrix, D: Expr) -> Matrix:
    """x(t) = x₀ − t · ∇D/|∇D| at x₀ (inward along −g)."""
    g0 = gradient(D).subs({x: x0[0], y: x0[1], z: x0[2]})
    n0 = g0 / sp.sqrt(g0.dot(g0))
    return x0 - t * n0


def curvature_along_march(k1: Expr, x_march: Matrix) -> Tuple[Expr, Expr, Expr]:
    """κ₁(t), dκ₁/dt, d²κ₁/dt² along march x(t)."""
    k1_t = sp.simplify(k1.subs({x: x_march[0], y: x_march[1], z: x_march[2]}))
    dk = sp.diff(k1_t, t)
    d2k = sp.diff(dk, t)
    return sp.simplify(k1_t), sp.simplify(dk), sp.simplify(d2k)
