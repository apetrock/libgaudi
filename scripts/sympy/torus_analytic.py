"""Parametric torus ground truth for curvature and inward-ray validation."""

from __future__ import annotations

from typing import List, Tuple

import numpy as np

from darboux import torus_normal, torus_point


def torus_principal_curvatures(u: float, v: float, major_R: float, minor_r: float) -> Tuple[float, float]:
    """κ_v (tube), κ_u (major); |κ_v| ≥ |κ_u| at outer equator."""
    cv = np.cos(v)
    denom = major_R + minor_r * cv
    k_v = cv / minor_r
    k_u = cv / denom if abs(denom) > 1e-14 else 0.0
    if abs(k_v) >= abs(k_u):
        return float(k_v), float(k_u)
    return float(k_u), float(k_v)


def radii_sum_at_uv(u: float, v: float, major_R: float, minor_r: float, eps_k: float = 1e-6) -> float:
    k1, k2 = torus_principal_curvatures(u, v, major_R, minor_r)
    return 1.0 / (abs(k1) + eps_k) + 1.0 / (abs(k2) + eps_k)


def torus_inward_ray_outer_equator(major_R: float, minor_r: float) -> Tuple[np.ndarray, np.ndarray]:
    """Foot at u=v=0; inward = -outward parametric normal."""
    p0 = torus_point(0.0, 0.0, major_R, minor_r)
    outward = torus_normal(p0, major_R)
    return p0, -outward


def arc_profile_inward(
    major_R: float = 1.25,
    minor_r: float = 0.35,
    steps: int = 25,
    eps_k: float = 1e-6,
) -> List[dict]:
    """On-surface arc v: 0 → π/2 (outer equator toward waist)."""
    rows: List[dict] = []
    for i in range(steps + 1):
        v = (np.pi / 2) * i / steps
        p = torus_point(0.0, v, major_R, minor_r)
        k1, k2 = torus_principal_curvatures(0.0, v, major_R, minor_r)
        rows.append(
            {
                "v": v,
                "x": p,
                "k1": k1,
                "k2": k2,
                "inv_k1": 1.0 / max(abs(k1), 1e-14),
                "inv_k2": 1.0 / max(abs(k2), 1e-14),
                "E_rad": radii_sum_at_uv(0.0, v, major_R, minor_r, eps_k),
            }
        )
    return rows
