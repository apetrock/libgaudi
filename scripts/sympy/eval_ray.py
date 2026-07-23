"""Numeric inv_frob energy evaluation along a fixed ray."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np
import sympy as sp
from sympy import Matrix

from compose import DEFAULT_EPS, H_directional, H_directional_rate, compose_energy_abstract


def shape_operator_np(g: np.ndarray, H: np.ndarray) -> np.ndarray:
    s2 = float(np.dot(g, g))
    if s2 < 1e-28:
        raise ValueError("zero gradient scale")
    s = np.sqrt(s2)
    n = g / s
    P = np.eye(3) - np.outer(n, n)
    return P @ H @ P / s


def frobenius_sq_np(W: np.ndarray) -> float:
    return float(np.sum(W * W))


def inv_frob_np(W: np.ndarray, eps: float = DEFAULT_EPS) -> float:
    return 1.0 / (frobenius_sq_np(W) + eps)


def principal_curvatures_np(g: np.ndarray, H: np.ndarray) -> Tuple[float, float]:
    """Debug/profile only — not used in sympy derive."""
    W = shape_operator_np(g, H)
    vals = np.linalg.eigvalsh(W)
    vals = sorted(vals, key=lambda v: abs(v), reverse=True)
    return float(vals[0]), float(vals[1])


_ENERGY_LAMBDA = None


def _build_energy_lambda():
    e = compose_energy_abstract()
    g = e["g"]
    H = e["H"]
    g_dot, H_dot, g_ddot, H_ddot = e["g_dot"], e["H_dot"], e["g_ddot"], e["H_ddot"]

    syms: List[sp.Symbol] = []
    for i in range(3):
        syms.append(g[i])
    for i in range(3):
        for j in range(3):
            syms.append(H[i, j])
    for i in range(3):
        syms.append(g_dot[i])
    for i in range(3):
        for j in range(3):
            syms.append(H_dot[i, j])
    for i in range(3):
        syms.append(g_ddot[i])
    for i in range(3):
        for j in range(3):
            syms.append(H_ddot[i, j])
    syms.append(e["eps"])

    return sp.lambdify(
        syms,
        [e["E"], e["E_prime"], e["E_pprime"]],
        modules="numpy",
    )


def _energy_lambda():
    global _ENERGY_LAMBDA
    if _ENERGY_LAMBDA is not None:
        return _ENERGY_LAMBDA

    _ENERGY_LAMBDA = _build_energy_lambda()
    return _ENERGY_LAMBDA


def _pack_kinematics(
    g: np.ndarray,
    H: np.ndarray,
    g_dot: np.ndarray,
    H_dot: np.ndarray,
    g_ddot: np.ndarray,
    H_ddot: np.ndarray,
) -> Tuple[float, ...]:
    args: List[float] = []
    args.extend(float(x) for x in g)
    args.extend(float(H[i, j]) for i in range(3) for j in range(3))
    args.extend(float(x) for x in g_dot)
    args.extend(float(H_dot[i, j]) for i in range(3) for j in range(3))
    args.extend(float(x) for x in g_ddot)
    args.extend(float(H_ddot[i, j]) for i in range(3) for j in range(3))
    return tuple(args)


def eval_energy_kinematics(
    g: np.ndarray,
    H: np.ndarray,
    g_dot: np.ndarray,
    H_dot: np.ndarray,
    g_ddot: np.ndarray,
    H_ddot: np.ndarray,
    *,
    eps: float = DEFAULT_EPS,
) -> Tuple[float, float, float]:
    fn = _energy_lambda()
    out = fn(*_pack_kinematics(g, H, g_dot, H_dot, g_ddot, H_ddot), float(eps))
    return float(out[0]), float(out[1]), float(out[2])


@dataclass(frozen=True)
class LineRay:
    p0: np.ndarray
    n0: np.ndarray

    def f(self, t: float) -> np.ndarray:
        return self.p0 + t * self.n0

    def fp(self, t: float) -> np.ndarray:
        del t
        return self.n0.copy()

    def fpp(self, t: float) -> np.ndarray:
        del t
        return np.zeros(3)


class QuadricGeom:
    def __init__(self, coeffs: Sequence[float]):
        self.coeffs = np.asarray(coeffs, dtype=float)

    def G(self, x: np.ndarray) -> np.ndarray:
        Q = self.coeffs
        return np.array(
            [
                2 * Q[0] * x[0] + Q[3] * x[1] + Q[4] * x[2] + 2 * Q[6],
                Q[3] * x[0] + 2 * Q[1] * x[1] + Q[5] * x[2] + 2 * Q[7],
                Q[4] * x[0] + Q[5] * x[1] + 2 * Q[2] * x[2] + 2 * Q[8],
            ]
        )

    def H(self, x: np.ndarray) -> np.ndarray:
        del x
        Q = self.coeffs
        return np.array(
            [
                [2 * Q[0], Q[3], Q[4]],
                [Q[3], 2 * Q[1], Q[5]],
                [Q[4], Q[5], 2 * Q[2]],
            ]
        )

    def H_directional(self, x: np.ndarray, fp: np.ndarray) -> np.ndarray:
        del x, fp
        return np.zeros((3, 3))

    def H_directional_rate(self, x: np.ndarray, fp: np.ndarray, fpp: np.ndarray) -> np.ndarray:
        del x, fp, fpp
        return np.zeros((3, 3))


class DarbouxGeom:
    def __init__(self, coeffs: Sequence[float]):
        self.coeffs = list(coeffs)
        from darboux import H_expr

        self._H = H_expr()
        self._G = __import__("darboux", fromlist=["G_expr"]).G_expr()
        from implicit import x, y, z

        self._pos = (x, y, z)

    def _subs_point(self, x: np.ndarray) -> dict:
        from darboux import C

        subs = {C[i]: self.coeffs[i] for i in range(14)}
        subs.update({self._pos[0]: float(x[0]), self._pos[1]: float(x[1]), self._pos[2]: float(x[2])})
        return subs

    def D(self, x: np.ndarray) -> float:
        from darboux import D_expr

        return float(D_expr().subs(self._subs_point(x)))

    def G(self, x: np.ndarray) -> np.ndarray:
        g = self._G.subs(self._subs_point(x))
        return np.array([float(g[i]) for i in range(3)], dtype=float)

    def H(self, x: np.ndarray) -> np.ndarray:
        Hm = self._H.subs(self._subs_point(x))
        return np.array([[float(Hm[i, j]) for j in range(3)] for i in range(3)])

    def H_directional(self, x: np.ndarray, fp: np.ndarray) -> np.ndarray:
        fp_sym = Matrix(sp.symbols("fp0:3", real=True))
        H_dir = H_directional(self._H, fp_sym)
        subs = self._subs_point(x)
        for i in range(3):
            subs[fp_sym[i]] = float(fp[i])
        Hd = H_dir.subs(subs)
        return np.array([[float(Hd[i, j]) for j in range(3)] for i in range(3)])

    def H_directional_rate(self, x: np.ndarray, fp: np.ndarray, fpp: np.ndarray) -> np.ndarray:
        fp_sym = Matrix(sp.symbols("fp0:3", real=True))
        fpp_sym = Matrix(sp.symbols("fpp0:3", real=True))
        H_rate = H_directional_rate(self._H, fp_sym, fpp_sym)
        subs = self._subs_point(x)
        for i in range(3):
            subs[fp_sym[i]] = float(fp[i])
            subs[fpp_sym[i]] = float(fpp[i])
        Hd = H_rate.subs(subs)
        return np.array([[float(Hd[i, j]) for j in range(3)] for i in range(3)])


def march_kinematics_at_t(geom, line: LineRay, t: float):
    x = line.f(t)
    fp = line.fp(t)
    fpp = line.fpp(t)
    g = geom.G(x)
    H = geom.H(x)
    H_dot = geom.H_directional(x, fp)
    H_ddot = geom.H_directional_rate(x, fp, fpp)
    g_dot = H @ fp
    g_ddot = H_dot @ fp + H @ fpp
    return g, H, g_dot, H_dot, g_ddot, H_ddot


def eval_energy_at_t(
    geom,
    line: LineRay,
    t: float,
    *,
    eps: float = DEFAULT_EPS,
) -> Tuple[float, float, float]:
    g, H, g_dot, H_dot, g_ddot, H_ddot = march_kinematics_at_t(geom, line, t)
    try:
        return eval_energy_kinematics(
            g, H, g_dot, H_dot, g_ddot, H_ddot, eps=eps
        )
    except (ValueError, FloatingPointError):
        def E_at(dt: float) -> float:
            xt = line.f(t + dt)
            fp = line.fp(t)
            fpp = line.fpp(t)
            gs = geom.G(xt)
            Hs = geom.H(xt)
            Ws = shape_operator_np(gs, Hs)
            return inv_frob_np(Ws, eps=eps)

        h = 1e-7
        e0 = E_at(0.0)
        ep = (E_at(h) - E_at(-h)) / (2 * h)
        epp = (E_at(h) - 2 * e0 + E_at(-h)) / (h * h)
        return e0, ep, epp


def ray_from_foot(geom, p0: np.ndarray) -> LineRay:
    """Ray along ∇D/|∇D|; medial may lie at t* < 0 or t* > 0."""
    g = geom.G(p0)
    return LineRay(np.asarray(p0, dtype=float), g / np.linalg.norm(g))


def inward_ray_from_foot(geom, p0: np.ndarray, **kwargs) -> LineRay:
    del kwargs
    return ray_from_foot(geom, p0)


def make_eval_energy(
    geom,
    line: LineRay,
    *,
    eps: float = DEFAULT_EPS,
) -> Callable[[float], Tuple[float, float, float]]:
    def eval_all(t: float) -> Tuple[float, float, float]:
        return eval_energy_at_t(geom, line, t, eps=eps)

    return eval_all
