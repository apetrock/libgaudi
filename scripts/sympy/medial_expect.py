"""Expected medial signatures for validation at x(t*)."""

from __future__ import annotations

from typing import Any, Dict

import numpy as np

from eval_ray import principal_curvatures_np


def inspect_medial(geom, x_medial: np.ndarray) -> Dict[str, Any]:
    g = geom.G(x_medial)
    H = geom.H(x_medial)
    gnorm = float(np.linalg.norm(g))
    if gnorm < 1e-12:
        return {
            "k1": float("inf"),
            "k2": float("inf"),
            "inv_k1": 0.0,
            "inv_k2": 0.0,
            "gnorm": gnorm,
        }
    k1, k2 = principal_curvatures_np(g, H)
    return {
        "k1": k1,
        "k2": k2,
        "inv_k1": 1.0 / max(abs(k1), 1e-14),
        "inv_k2": 1.0 / max(abs(k2), 1e-14),
        "gnorm": gnorm,
    }


def inspect_medial_near(
    geom,
    line,
    t_star: float,
    *,
    gnorm_min: float = 1e-4,
) -> Dict[str, Any]:
    """Inspect κ at medial; step slightly back if |g|≈0 at t*."""
    for scale in (1.0, 0.98, 0.95, 0.92, 0.88):
        x = line.f(t_star * scale)
        out = inspect_medial(geom, x)
        if out["gnorm"] >= gnorm_min:
            return out
    return inspect_medial(geom, line.f(t_star))
