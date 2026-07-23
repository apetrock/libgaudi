"""Composed medial-energy graph and abstract kinematic jet bundle."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import sympy as sp
from sympy import Matrix

from cache_utils import CACHE_DIR, disk_cache, energy_cache_key
from implicit import x, y, z
from line import f_expr, t
from shape_fcn import abstract_gh_symbols, shape_operator

DEFAULT_EPS = 1e-12

_GRAPH_CACHE: Dict[str, "MedialGraph"] = {}
_ENERGY_CACHE: Optional[Dict[str, sp.Basic]] = None


@dataclass(frozen=True)
class MedialGraph:
    """Composed march: E = 1 / (trace(W²) + ε)."""

    f: Matrix
    g: Matrix
    H: Matrix
    W: Matrix
    s2: sp.Expr
    E: sp.Expr
    E_prime: sp.Expr
    E_pprime: sp.Expr
    eps: sp.Symbol


def subs_GH_at_f(G_expr: Matrix, H_expr: Matrix, f: Matrix) -> Tuple[Matrix, Matrix]:
    pos = {x: f[0], y: f[1], z: f[2]}
    return G_expr.subs(pos), H_expr.subs(pos)


def frobenius_sq(W: Matrix) -> sp.Expr:
    return sp.trace(W * W)


def energy_from_W(W: Matrix, eps: sp.Basic, mode: str = "inv_frob") -> sp.Expr:
    s2 = frobenius_sq(W)
    if mode == "frob":
        return s2
    if mode == "inv_frob":
        return 1 / (s2 + eps)
    raise ValueError(f"unknown energy mode: {mode}")


def abstract_kinematic_symbols() -> Tuple[Matrix, Matrix, Matrix, Matrix]:
    g_dot = Matrix(sp.symbols("gd0:3", real=True))
    H_dot = Matrix(3, 3, lambda i, j: sp.symbols(f"Hd{i}{j}", real=True))
    g_ddot = Matrix(sp.symbols("gdd0:3", real=True))
    H_ddot = Matrix(3, 3, lambda i, j: sp.symbols(f"Hdd{i}{j}", real=True))
    return g_dot, H_dot, g_ddot, H_ddot


def shape_operator_jet(
    g: Matrix,
    H: Matrix,
    g_dot: Matrix,
    H_dot: Matrix,
    g_ddot: Matrix,
    H_ddot: Matrix,
) -> Tuple[Matrix, Matrix, Matrix]:
    eps_p = sp.Symbol("eps", real=True)
    W = shape_operator(g, H)
    W_eps = shape_operator(g + eps_p * g_dot, H + eps_p * H_dot)
    W_dot = sp.diff(W_eps, eps_p).subs(eps_p, 0)
    W_eps2 = shape_operator(
        g + eps_p * g_dot + eps_p**2 * g_ddot / 2,
        H + eps_p * H_dot + eps_p**2 * H_ddot / 2,
    )
    W_ddot = sp.diff(W_eps2, eps_p, 2).subs(eps_p, 0)
    return W, W_dot, W_ddot


def energy_jet(
    W: Matrix,
    W_dot: Matrix,
    W_ddot: Matrix,
    eps: sp.Basic,
    mode: str = "inv_frob",
) -> Tuple[sp.Expr, sp.Expr, sp.Expr]:
    """E, dE/dt, d²E/dt² via ε-perturbation on W(t)."""
    eps_p = sp.Symbol("eps", real=True)
    E0 = energy_from_W(W, eps, mode=mode)
    E1 = energy_from_W(W + eps_p * W_dot, eps, mode=mode)
    E_prime = sp.diff(E1, eps_p).subs(eps_p, 0)
    E2 = energy_from_W(W + eps_p * W_dot + eps_p**2 * W_ddot / 2, eps, mode=mode)
    E_pprime = sp.diff(E2, eps_p, 2).subs(eps_p, 0)
    return E0, E_prime, E_pprime


def H_directional(H: Matrix, v: Matrix) -> Matrix:
    pos = (x, y, z)
    rows = []
    for i in range(3):
        row = []
        for j in range(3):
            val = sum(sp.diff(H[i, j], pos[k]) * v[k] for k in range(3))
            row.append(val)
        rows.append(row)
    return Matrix(rows)


def H_directional_rate(H: Matrix, fp: Matrix, fpp: Matrix) -> Matrix:
    H_dot = H_directional(H, fp)
    return H_directional(H_dot, fp) + H_directional(H, fpp)


def medial_energy_graph(G_expr: Matrix, H_expr: Matrix) -> MedialGraph:
    """Compose f(t) → G,H → W → E = 1/(trace(W²)+ε); diff w.r.t. t."""
    eps = sp.Symbol("eps", real=True, positive=True)

    f = f_expr()
    g, H = subs_GH_at_f(G_expr, H_expr, f)
    W = shape_operator(g, H)
    s2 = frobenius_sq(W)
    E = energy_from_W(W, eps, mode="inv_frob")
    E_prime = sp.diff(E, t)
    E_pprime = sp.diff(E_prime, t)

    return MedialGraph(
        f=f,
        g=g,
        H=H,
        W=W,
        s2=s2,
        E=E,
        E_prime=E_prime,
        E_pprime=E_pprime,
        eps=eps,
    )


def compose_quadric_graph() -> MedialGraph:
    from quadric import G_expr, H_expr

    return medial_energy_graph(G_expr(), H_expr())


def compose_darboux_graph() -> MedialGraph:
    from darboux import G_expr, H_expr

    return medial_energy_graph(G_expr(), H_expr())


def get_graph(name: str) -> MedialGraph:
    if name in _GRAPH_CACHE:
        return _GRAPH_CACHE[name]

    builders = {
        "quadric": compose_quadric_graph,
        "darboux": compose_darboux_graph,
    }
    if name not in builders:
        raise ValueError(f"unknown graph: {name}")

    cache_path = CACHE_DIR / f"march_graph_{name}_{energy_cache_key()}.pkl"
    graph = disk_cache(cache_path, builders[name])
    _GRAPH_CACHE[name] = graph
    return graph


def _derive_energy_abstract() -> Dict[str, sp.Basic]:
    g, H = abstract_gh_symbols()
    g_dot, H_dot, g_ddot, H_ddot = abstract_kinematic_symbols()
    eps = sp.Symbol("eps", real=True, positive=True)

    W, W_dot, W_ddot = shape_operator_jet(g, H, g_dot, H_dot, g_ddot, H_ddot)
    E, E_prime, E_pprime = energy_jet(W, W_dot, W_ddot, eps, mode="inv_frob")

    return {
        "g": g,
        "H": H,
        "g_dot": g_dot,
        "H_dot": H_dot,
        "g_ddot": g_ddot,
        "H_ddot": H_ddot,
        "eps": eps,
        "W": W,
        "W_dot": W_dot,
        "W_ddot": W_ddot,
        "s2": frobenius_sq(W),
        "E": E,
        "E_prime": E_prime,
        "E_pprime": E_pprime,
    }


def compose_energy_abstract() -> Dict[str, sp.Basic]:
    """Abstract jet bundle for lambdify and C++ emit."""
    global _ENERGY_CACHE
    if _ENERGY_CACHE is not None:
        return _ENERGY_CACHE

    cache_path = CACHE_DIR / f"energy_abstract_{energy_cache_key()}.pkl"
    _ENERGY_CACHE = disk_cache(cache_path, _derive_energy_abstract)
    return _ENERGY_CACHE
