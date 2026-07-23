"""Dump scalar profiles along a bidirectional march ray."""

from __future__ import annotations

from typing import List, Optional, TextIO

import numpy as np

from compose import DEFAULT_EPS
from eval_ray import (
    DarbouxGeom,
    QuadricGeom,
    eval_energy_at_t,
    frobenius_sq_np,
    principal_curvatures_np,
    ray_from_foot,
    shape_operator_np,
)
from quadric import sphere_coefficients


def _profile_ray(
    geom,
    p0: np.ndarray,
    t_max: float,
    steps: int,
    eps: float,
    outfile: Optional[TextIO],
) -> List[dict]:
    line = ray_from_foot(geom, p0)
    header = "t\tx\ty\tz\t|g|\tk1\tk2\ts2\tE_inv_frob\tE'\tE''"
    print(header, file=outfile)

    rows: List[dict] = []
    for i in range(steps + 1):
        frac = 2.0 * i / steps - 1.0
        t = t_max * frac
        x = line.f(t)
        try:
            g = geom.G(x)
            H = geom.H(x)
            gnorm = float(np.linalg.norm(g))
            k1, k2 = principal_curvatures_np(g, H)
            W = shape_operator_np(g, H)
            s2 = frobenius_sq_np(W)
            e, e_p, e_pp = eval_energy_at_t(geom, line, t, eps=eps)
        except (ValueError, FloatingPointError):
            continue

        row = {
            "t": t,
            "x": x,
            "gnorm": gnorm,
            "k1": k1,
            "k2": k2,
            "s2": s2,
            "e": e,
            "e_p": e_p,
            "e_pp": e_pp,
        }
        rows.append(row)
        print(
            f"{t:.4f}\t{x[0]:.4f}\t{x[1]:.4f}\t{x[2]:.4f}\t"
            f"{gnorm:.6f}\t{k1:.6f}\t{k2:.6f}\t{s2:.6f}\t"
            f"{e:.6f}\t{e_p:.6f}\t{e_pp:.6f}",
            file=outfile,
        )
    return rows


def profile_sphere_ray(
    R: float = 2.5,
    t_max: Optional[float] = None,
    steps: int = 25,
    *,
    eps: float = DEFAULT_EPS,
    outfile: Optional[TextIO] = None,
) -> List[dict]:
    t_max = float(R if t_max is None else t_max)
    geom = QuadricGeom(sphere_coefficients(R))
    p0 = np.array([R, 0.0, 0.0])
    return _profile_ray(geom, p0, t_max, steps, eps, outfile)


def profile_torus_ray(
    major_R: float = 1.25,
    minor_r: float = 0.35,
    t_max: Optional[float] = None,
    steps: int = 25,
    *,
    eps: float = DEFAULT_EPS,
    outfile: Optional[TextIO] = None,
) -> List[dict]:
    from darboux import canonical_torus_coefficients, torus_point

    t_max = float(minor_r if t_max is None else t_max)
    geom = DarbouxGeom(canonical_torus_coefficients())
    p0 = torus_point(0.0, 0.0, major_R, minor_r)
    return _profile_ray(geom, p0, t_max, steps, eps, outfile)
