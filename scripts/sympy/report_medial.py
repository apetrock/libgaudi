#!/usr/bin/env python3
"""Report expected vs found medial search results."""

from __future__ import annotations

import numpy as np

from compose import DEFAULT_EPS
from darboux import canonical_torus_coefficients, torus_point
from eval_ray import (
    DarbouxGeom,
    QuadricGeom,
    frobenius_sq_np,
    inv_frob_np,
    make_eval_energy,
    principal_curvatures_np,
    ray_from_foot,
    shape_operator_np,
)
from medial_expect import inspect_medial, inspect_medial_near
from quadric import ellipsoid_coefficients, sphere_coefficients
from search import newton_min_energy


def _row(label: str, expected: str, found: str) -> str:
    return f"{label:<28} {expected:>18} {found:>18}"


def _fmt(x: float) -> str:
    if abs(x) > 1e6:
        return f"{x:.3e}"
    return f"{x:.6f}"


def report_sphere(R: float = 2.5) -> None:
    geom = QuadricGeom(sphere_coefficients(R))
    p0 = np.array([R, 0.0, 0.0])
    g0, H0 = geom.G(p0), geom.H(p0)
    k1f, k2f = principal_curvatures_np(g0, H0)
    W0 = shape_operator_np(g0, H0)
    e0 = inv_frob_np(W0, eps=DEFAULT_EPS)

    line = ray_from_foot(geom, p0)
    eval_all = make_eval_energy(geom, line)
    t_star, e_star, ok = newton_min_energy(eval_all, t0=-0.1 * R, t_max=R)
    x_medial = line.f(t_star)
    medial = inspect_medial(geom, x_medial)

    print(f"\n=== Sphere R={R} (QuadricGeom) ===")
    print(_row("metric", "expected", "found"))
    print(_row("foot |1/k1|", _fmt(R), _fmt(1 / abs(k1f))))
    print(_row("foot |1/k2|", _fmt(R), _fmt(1 / abs(k2f))))
    print(_row("foot E=1/(k1²+k2²)", _fmt(R**2 / 2), _fmt(e0)))
    print(_row("|t*|", _fmt(R), _fmt(abs(t_star))))
    print(_row("|x(t*)|", "0", _fmt(float(np.linalg.norm(x_medial)))))
    print(_row("medial |g|", "~0", _fmt(medial["gnorm"])))
    print(_row("Newton converged", "True", str(ok)))


def report_ellipsoid(a: float = 3.0, b: float = 1.5, c: float = 1.0) -> None:
    geom = QuadricGeom(ellipsoid_coefficients(a, b, c))
    for axis, p0, t_exp in [
        ("x", np.array([a, 0.0, 0.0]), a),
        ("y", np.array([0.0, b, 0.0]), b),
    ]:
        line = ray_from_foot(geom, p0)
        eval_all = make_eval_energy(geom, line)
        t_star, _, ok = newton_min_energy(eval_all, t0=-0.1 * t_exp, t_max=t_exp)
        x_medial = line.f(t_star)
        medial = inspect_medial(geom, x_medial)
        g0, H0 = geom.G(p0), geom.H(p0)
        k1f, k2f = principal_curvatures_np(g0, H0)

        print(f"\n=== Ellipsoid axis={axis} a={a} b={b} c={c} ===")
        print(_row("metric", "expected", "found"))
        print(_row("foot |1/k1|", f"~{axis} radius", _fmt(1 / abs(k1f))))
        print(_row("foot |1/k2|", f"~{axis} radius", _fmt(1 / abs(k2f))))
        print(_row("|t*|", _fmt(t_exp), _fmt(abs(t_star))))
        print(_row("|x(t*)|", "0", _fmt(float(np.linalg.norm(x_medial)))))
        print(_row("medial |g|", "~0", _fmt(medial["gnorm"])))
        print(_row("Newton converged", "True", str(ok)))


def report_torus(major_R: float = 1.25, minor_r: float = 0.35) -> None:
    geom = DarbouxGeom(canonical_torus_coefficients())
    p0 = torus_point(0.0, 0.0, major_R, minor_r)
    g0, H0 = geom.G(p0), geom.H(p0)
    k1f, k2f = principal_curvatures_np(g0, H0)

    line = ray_from_foot(geom, p0)
    eval_all = make_eval_energy(geom, line)
    t_star, e_star, ok = newton_min_energy(eval_all, t0=0.05 * minor_r, t_max=minor_r)
    x_medial = line.f(t_star)
    waist = np.array([major_R, 0.0, 0.0])
    medial = inspect_medial_near(geom, line, t_star)
    inv_ks = sorted([medial["inv_k1"], medial["inv_k2"]], reverse=True)

    print(f"\n=== Torus outer equator R={major_R} r={minor_r} (DarbouxGeom) ===")
    print(_row("metric", "expected", "found"))
    print(_row("foot |1/k1|", _fmt(minor_r), _fmt(1 / abs(k1f))))
    print(_row("foot |1/k2|", _fmt(major_R + minor_r), _fmt(1 / abs(k2f))))
    print(_row("|t*|", _fmt(minor_r), _fmt(abs(t_star))))
    print(_row("|x-waist|", "0", _fmt(float(np.linalg.norm(x_medial - waist)))))
    print(_row("medial inv_k (spine)", f"~{major_R}", _fmt(inv_ks[0])))
    print(_row("medial inv_k (tube)", "~0", _fmt(inv_ks[1])))
    print(_row("medial k1", "large", _fmt(medial["k1"])))
    print(_row("medial k2", f"~-1/{major_R}", _fmt(medial["k2"])))
    print(_row("Newton converged", "True", str(ok)))


def main() -> None:
    print("Medial search: expected vs found (pure W, hybrid Newton)")
    report_sphere()
    report_ellipsoid()
    report_torus()


if __name__ == "__main__":
    main()
