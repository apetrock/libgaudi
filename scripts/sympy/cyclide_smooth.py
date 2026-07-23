"""Jet-matching smooth energy for Darboux cyclide fits (Phase 1).

Per-vertex local energy with one neighbor slot (C++ sums over mesh neighbors):

  E_anchor   = wi * ||Qi - Q0||^2
  E_neighbor = (1-wi) * w_ij * [ alpha_G * frob(G_i - G_j)^2
                               + alpha_H * frob(H_i - H_j)^2 ]

G/H are Darboux gradient/Hessian, linear in Q coefficients → quadratic energy,
gradient linear in Qi → normal equations are a 14x14 linear solve.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import sympy as sp
from sympy import Matrix

from cache_utils import CACHE_DIR, disk_cache
from darboux import C, G_expr, H_expr
from implicit import x, y, z

N_Q = 14
CoeffSyms = Tuple[sp.Symbol, ...]

_BUNDLE_CACHE: Optional["CyclideSmoothBundle"] = None


def cyclide_smooth_cache_key() -> str:
    root = Path(__file__).resolve().parent
    h = hashlib.sha256()
    for name in ("cyclide_smooth.py", "darboux.py", "implicit.py"):
        h.update((root / name).read_bytes())
    h.update(sp.__version__.encode())
    return h.hexdigest()[:16]


def _coeff_symbols(tag: str) -> CoeffSyms:
    return sp.symbols(f"{tag}0:14", real=True)


def subs_Q(expr: sp.Basic, Q_syms: Sequence[sp.Symbol]) -> sp.Basic:
    return expr.subs({C[k]: Q_syms[k] for k in range(N_Q)})


def G_at(Q_syms: Sequence[sp.Symbol], pos: Matrix) -> Matrix:
    G = subs_Q(G_expr(), Q_syms)
    return G.subs({x: pos[0], y: pos[1], z: pos[2]})


def H_at(Q_syms: Sequence[sp.Symbol], pos: Matrix) -> Matrix:
    H = subs_Q(H_expr(), Q_syms)
    return H.subs({x: pos[0], y: pos[1], z: pos[2]})


def frobenius_sq(M: Matrix) -> sp.Expr:
    return sum(M[i, j] ** 2 for i in range(3) for j in range(3))


def linear_map_from_vector(expr_vec: Matrix, Q_syms: Sequence[sp.Symbol]) -> Matrix:
    """expr_vec = M @ Q when expr is linear in Q."""
    n = len(Q_syms)
    return Matrix(
        [
            [sp.expand(expr_vec[i].diff(Q_syms[k])) for k in range(n)]
            for i in range(expr_vec.rows)
        ]
    )


def flatten_H(H: Matrix) -> Matrix:
    return Matrix([H[i, j] for i in range(3) for j in range(3)])


@dataclass(frozen=True)
class CyclideSmoothBundle:
    Qi: CoeffSyms
    Q0: CoeffSyms
    Qj: CoeffSyms
    x_j_in_i: Matrix
    wi: sp.Symbol
    w_ij: sp.Symbol
    alpha_G: sp.Symbol
    alpha_H: sp.Symbol
    G_i: Matrix
    G_j: Matrix
    H_i: Matrix
    H_j: Matrix
    M_G_i: Matrix
    M_H_i: Matrix
    E_anchor: sp.Expr
    E_neighbor: sp.Expr
    grad_anchor: Matrix
    grad_neighbor: Matrix
    hess_anchor: Matrix
    hess_neighbor: Matrix


def _derive_cyclide_smooth() -> CyclideSmoothBundle:
    Qi = _coeff_symbols("Qi")
    Q0 = _coeff_symbols("Qg")
    Qj = _coeff_symbols("Qj")
    x_j_in_i = Matrix(sp.symbols("xj0:3", real=True))
    wi, w_ij, alpha_G, alpha_H = sp.symbols(
        "wi w_ij alpha_G alpha_H", real=True
    )

    G_i = G_at(Qi, x_j_in_i)
    G_j = G_at(Qj, Matrix([0, 0, 0]))
    H_i = H_at(Qi, x_j_in_i)
    H_j = H_at(Qj, Matrix([0, 0, 0]))

    M_G_i = linear_map_from_vector(G_i, Qi)
    M_H_i = linear_map_from_vector(flatten_H(H_i), Qi)

    dG = G_i - G_j
    dH = H_i - H_j

    E_anchor = wi * sum((Qi[k] - Q0[k]) ** 2 for k in range(N_Q))
    E_neighbor = (1 - wi) * w_ij * (
        alpha_G * dG.dot(dG) + alpha_H * frobenius_sq(dH)
    )

    grad_anchor = Matrix([sp.expand(E_anchor.diff(Qi[k])) for k in range(N_Q)])
    grad_neighbor = Matrix(
        [sp.expand(E_neighbor.diff(Qi[k])) for k in range(N_Q)]
    )

    hess_anchor = Matrix(
        [
            [sp.expand(grad_anchor[i].diff(Qi[j])) for j in range(N_Q)]
            for i in range(N_Q)
        ]
    )
    hess_neighbor = Matrix(
        [
            [sp.expand(grad_neighbor[i].diff(Qi[j])) for j in range(N_Q)]
            for i in range(N_Q)
        ]
    )

    return CyclideSmoothBundle(
        Qi=Qi,
        Q0=Q0,
        Qj=Qj,
        x_j_in_i=x_j_in_i,
        wi=wi,
        w_ij=w_ij,
        alpha_G=alpha_G,
        alpha_H=alpha_H,
        G_i=G_i,
        G_j=G_j,
        H_i=H_i,
        H_j=H_j,
        M_G_i=M_G_i,
        M_H_i=M_H_i,
        E_anchor=E_anchor,
        E_neighbor=E_neighbor,
        grad_anchor=grad_anchor,
        grad_neighbor=grad_neighbor,
        hess_anchor=hess_anchor,
        hess_neighbor=hess_neighbor,
    )


def compose_cyclide_smooth() -> CyclideSmoothBundle:
    global _BUNDLE_CACHE
    if _BUNDLE_CACHE is not None:
        return _BUNDLE_CACHE
    cache_path = CACHE_DIR / f"cyclide_smooth_{cyclide_smooth_cache_key()}.pkl"
    _BUNDLE_CACHE = disk_cache(cache_path, _derive_cyclide_smooth)
    return _BUNDLE_CACHE


INSPECT_KEYS: Dict[str, str] = {
    "E_anchor": "E_anchor",
    "E_neighbor": "E_neighbor",
    "grad_anchor": "grad_anchor",
    "grad_neighbor": "grad_neighbor",
    "hess_anchor": "hess_anchor",
    "hess_neighbor": "hess_neighbor",
    "M_G_i": "M_G_i",
    "M_H_i": "M_H_i",
    "G_i": "G_i",
    "H_i": "H_i",
}


def get_inspect_expr(bundle: CyclideSmoothBundle, key: str) -> sp.Basic:
    if key not in INSPECT_KEYS:
        raise KeyError(f"unknown key {key!r} (choose from {list(INSPECT_KEYS)})")
    return getattr(bundle, INSPECT_KEYS[key])
