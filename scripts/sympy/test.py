"""Tests for flat sympy medial pipeline."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent


def _quadric_grad_cpp(Q, x):
    return np.array(
        [
            2 * Q[0] * x[0] + 2 * Q[3] * x[1] + 2 * Q[4] * x[2] + 2 * Q[6],
            2 * Q[3] * x[0] + 2 * Q[1] * x[1] + 2 * Q[5] * x[2] + 2 * Q[7],
            2 * Q[4] * x[0] + 2 * Q[5] * x[1] + 2 * Q[2] * x[2] + 2 * Q[8],
        ]
    )


def test_compose_energy_abstract():
    from compose import compose_energy_abstract

    e = compose_energy_abstract()
    assert "E" in e
    assert "s2" in e
    assert "E_prime" in e
    assert "E_pprime" in e


def test_compose_matches_jet_sphere():
    from compose import DEFAULT_EPS, get_graph
    from eval_ray import QuadricGeom, eval_energy_at_t, inward_ray_from_foot
    from line import n0, p0, t
    from quadric import Q, sphere_coefficients

    R = 2.5
    t_val = 0.3 * R
    geom = QuadricGeom(sphere_coefficients(R))
    p0_np = np.array([R, 0.0, 0.0])
    line = inward_ray_from_foot(geom, p0_np)

    graph = get_graph("quadric")
    subs = {Q[i]: sphere_coefficients(R)[i] for i in range(10)}
    subs[graph.eps] = DEFAULT_EPS
    for i in range(3):
        subs[p0[i]] = float(p0_np[i])
        subs[n0[i]] = float(line.n0[i])
    subs[t] = t_val

    E_prime_sym = float(graph.E_prime.subs(subs))
    _, E_prime_num, _ = eval_energy_at_t(geom, line, t_val)
    assert abs(E_prime_sym - E_prime_num) < 1e-3


def test_quadric_grad_matches_cpp():
    from implicit import x, y, z
    from quadric import G_expr, Q, sphere_coefficients

    G = G_expr()
    coeffs = sphere_coefficients(2.5)
    subs = {Q[i]: coeffs[i] for i in range(10)}
    subs.update({x: 2.5, y: 0.0, z: 0.0})
    g_sym = np.array([float(G[i].subs(subs)) for i in range(3)])
    g_cpp = _quadric_grad_cpp(coeffs, np.array([2.5, 0.0, 0.0]))
    assert np.allclose(g_sym, g_cpp)


def test_inv_frob_at_sphere_foot():
    from compose import DEFAULT_EPS
    from eval_ray import QuadricGeom, frobenius_sq_np, inv_frob_np, shape_operator_np
    from quadric import sphere_coefficients

    R = 2.5
    geom = QuadricGeom(sphere_coefficients(R))
    p0 = np.array([R, 0.0, 0.0])
    g = geom.G(p0)
    H = geom.H(p0)
    W = shape_operator_np(g, H)
    s2 = frobenius_sq_np(W)
    assert abs(s2 - 2.0 / R**2) < 1e-6
    e = inv_frob_np(W, eps=DEFAULT_EPS)
    assert abs(e - R**2 / 2.0) < 1e-4


def test_trace_s2_equals_kappa_sum_sq():
    from darboux import canonical_torus_coefficients, torus_point
    from eval_ray import DarbouxGeom, frobenius_sq_np, principal_curvatures_np, shape_operator_np

    geom = DarbouxGeom(canonical_torus_coefficients())
    p0 = torus_point(0.0, 0.0, 1.25, 0.35)
    g = geom.G(p0)
    H = geom.H(p0)
    k1, k2 = principal_curvatures_np(g, H)
    W = shape_operator_np(g, H)
    assert abs(frobenius_sq_np(W) - (k1 * k1 + k2 * k2)) < 1e-4


def test_newton_step_matches_formula():
    from search import damped_newton_step, newton_step

    assert abs(newton_step(0.3, 0.1, 2.0) - 0.25) < 1e-12
    assert abs(damped_newton_step(0.3, 0.1, 2.0) - 0.25) < 1e-12
    assert abs(damped_newton_step(0.3, 0.1, 0.0) - 0.29) < 1e-12


def test_search_damped_newton_flat_hessian():
    from search import damped_newton_step

    t_next = damped_newton_step(1.0, 0.5, 0.0, min_denom=1e-14, lm_lambda=1e-8)
    assert t_next != 1.0
    assert t_next < 1.0


def test_search_bracket_fallback_near_singularity():
    from search import bracket_min_energy, newton_min_energy

    def flat_bowl(t: float):
        e = t * t + 0.01
        e_p = 2.0 * t
        e_pp = 2.0 if abs(t) > 1e-3 else 0.0
        return e, e_p, e_pp

    t_br, _ = bracket_min_energy(flat_bowl, -1.0, 1.0, n_samples=100)
    assert abs(t_br) < 0.1
    t_star, _, ok = newton_min_energy(flat_bowl, t0=0.5, t_max=1.0, max_iters=8)
    assert ok
    assert abs(t_star) < 0.2


def test_generated_has_inv_frob():
    path = ROOT / "generated" / "medial_kernels.h"
    if not path.exists():
        pytest.skip("run generate.py first")
    text = path.read_text()
    assert "inv_frob_derivs_from_kinematics" in text
    assert "eval_medial_energy_at_t" in text


def test_sphere_medial_newton():
    from eval_ray import QuadricGeom, make_eval_energy, ray_from_foot
    from medial_expect import inspect_medial_near
    from quadric import sphere_coefficients
    from search import newton_min_energy

    R = 2.5
    geom = QuadricGeom(sphere_coefficients(R))
    p0 = np.array([R, 0.0, 0.0])
    line = ray_from_foot(geom, p0)
    eval_all = make_eval_energy(geom, line)
    t_star, e_star, ok = newton_min_energy(eval_all, t0=-0.1 * R, t_max=R)
    assert ok
    assert abs(abs(t_star) - R) < 0.05
    x_medial = line.f(t_star)
    assert np.linalg.norm(x_medial) < 0.05
    e_foot, _, _ = eval_all(0.0)
    assert e_star < e_foot
    medial = inspect_medial_near(geom, line, t_star)
    assert medial["gnorm"] < 1e-3 or np.linalg.norm(x_medial) < 0.05


def test_sphere_inv_frob_decreases_inward():
    from eval_ray import QuadricGeom, make_eval_energy, ray_from_foot
    from quadric import sphere_coefficients

    R = 2.5
    geom = QuadricGeom(sphere_coefficients(R))
    line = ray_from_foot(geom, np.array([R, 0.0, 0.0]))
    eval_all = make_eval_energy(geom, line)
    e0, _, _ = eval_all(0.0)
    e_in, _, _ = eval_all(-0.5 * R)
    assert e_in < e0


def test_ellipsoid_medial_newton():
    from eval_ray import QuadricGeom, make_eval_energy, ray_from_foot
    from medial_expect import inspect_medial_near
    from quadric import ellipsoid_coefficients
    from search import newton_min_energy

    a, b, c = 3.0, 1.5, 1.0
    geom = QuadricGeom(ellipsoid_coefficients(a, b, c))
    p0 = np.array([a, 0.0, 0.0])
    line = ray_from_foot(geom, p0)
    eval_all = make_eval_energy(geom, line)
    t_star, _, ok = newton_min_energy(eval_all, t0=-0.1 * a, t_max=a)
    assert ok
    assert abs(abs(t_star) - a) < 0.08
    x_medial = line.f(t_star)
    assert np.linalg.norm(x_medial) < 0.08
    medial = inspect_medial_near(geom, line, t_star)
    assert medial["gnorm"] < 1e-3 or np.linalg.norm(x_medial) < 0.05


def test_ellipsoid_medial_newton_y_axis():
    from eval_ray import QuadricGeom, make_eval_energy, ray_from_foot
    from medial_expect import inspect_medial_near
    from quadric import ellipsoid_coefficients
    from search import newton_min_energy

    a, b, c = 3.0, 1.5, 1.0
    geom = QuadricGeom(ellipsoid_coefficients(a, b, c))
    p0 = np.array([0.0, b, 0.0])
    line = ray_from_foot(geom, p0)
    eval_all = make_eval_energy(geom, line)
    t_star, _, ok = newton_min_energy(eval_all, t0=-0.1 * b, t_max=b)
    assert ok
    assert abs(abs(t_star) - b) < 0.08
    x_medial = line.f(t_star)
    assert np.linalg.norm(x_medial) < 0.08
    medial = inspect_medial_near(geom, line, t_star)
    assert medial["gnorm"] < 1e-3 or np.linalg.norm(x_medial) < 0.05


def test_darboux_sphere_matches_quadric_foot():
    from compose import DEFAULT_EPS
    from eval_ray import DarbouxGeom, QuadricGeom, inv_frob_np, shape_operator_np
    from darboux import sphere_coefficients as darboux_sphere
    from quadric import sphere_coefficients

    R = 2.5
    p0 = np.array([R, 0.0, 0.0])
    Wq = shape_operator_np(QuadricGeom(sphere_coefficients(R)).G(p0), QuadricGeom(sphere_coefficients(R)).H(p0))
    Wd = shape_operator_np(DarbouxGeom(darboux_sphere(R)).G(p0), DarbouxGeom(darboux_sphere(R)).H(p0))
    assert abs(inv_frob_np(Wq, DEFAULT_EPS) - inv_frob_np(Wd, DEFAULT_EPS)) < 1e-4


def test_darboux_torus_surface_fit():
    from darboux import canonical_torus_coefficients, torus_point
    from eval_ray import DarbouxGeom

    geom = DarbouxGeom(canonical_torus_coefficients())
    errs = []
    for i in range(24):
        u = 2.0 * np.pi * i / 24.0
        for j in range(12):
            v = 2.0 * np.pi * j / 12.0
            p = torus_point(u, v, 1.25, 0.35)
            errs.append(abs(geom.D(p)))
    assert max(errs) < 1e-10


def test_torus_k_radii_at_foot():
    from darboux import canonical_torus_coefficients, torus_point
    from eval_ray import DarbouxGeom, principal_curvatures_np

    major_R, minor_r = 1.25, 0.35
    geom = DarbouxGeom(canonical_torus_coefficients())
    p0 = torus_point(0.0, 0.0, major_R, minor_r)
    g = geom.G(p0)
    H = geom.H(p0)
    k1, k2 = principal_curvatures_np(g, H)
    assert abs(1 / abs(k1) - minor_r) < 0.02
    assert abs(1 / abs(k2) - (major_R + minor_r)) < 0.05


def test_torus_gradient_parallel_to_parametric_normal():
    from darboux import canonical_torus_coefficients, torus_normal, torus_point
    from eval_ray import DarbouxGeom

    major_R = 1.25
    geom = DarbouxGeom(canonical_torus_coefficients())
    p0 = torus_point(0.0, 0.0, major_R, 0.35)
    g = geom.G(p0)
    outward = torus_normal(p0, major_R)
    dot = float(np.dot(g / np.linalg.norm(g), outward))
    assert abs(abs(dot) - 1.0) < 0.05


def test_torus_medial_newton():
    from darboux import canonical_torus_coefficients, torus_point
    from eval_ray import DarbouxGeom, make_eval_energy, ray_from_foot
    from medial_expect import inspect_medial_near
    from search import newton_min_energy

    major_R, minor_r = 1.25, 0.35
    geom = DarbouxGeom(canonical_torus_coefficients())
    p0 = torus_point(0.0, 0.0, major_R, minor_r)
    line = ray_from_foot(geom, p0)
    eval_all = make_eval_energy(geom, line)
    t_star, e_star, ok = newton_min_energy(eval_all, t0=0.05 * minor_r, t_max=minor_r)
    assert ok
    assert abs(abs(t_star) - minor_r) < 0.08
    x_medial = line.f(t_star)
    waist = np.array([major_R, 0.0, 0.0])
    assert np.linalg.norm(x_medial - waist) < 0.08
    e_foot, _, _ = eval_all(0.0)
    assert e_star < e_foot
    medial = inspect_medial_near(geom, line, t_star)
    inv_ks = sorted([medial["inv_k1"], medial["inv_k2"]], reverse=True)
    assert inv_ks[0] > 0.5 * major_R
    assert inv_ks[1] < 0.1 * minor_r


def test_generate_inspect_runs(capsys):
    from generate import inspect_symbols

    inspect_symbols("E", use_cse=False)
    captured = capsys.readouterr()
    assert "E" in captured.out


def test_emit_kernels_inv_frob():
    from albers_printer import emit_kernels

    text = emit_kernels()
    assert "inv_frob_derivs_from_kinematics" in text
    assert "eval_medial_energy_at_t" in text
    assert "inv_frob_derivs_from_W_jets" not in text
    assert "sigma" not in text
    assert "objective_log_inv_from_e2" not in text


def test_emit_has_cse_blob():
    from albers_printer import emit_kernels

    text = emit_kernels()
    assert "const real x0 =" in text
    assert "inv_frob_derivs_from_kinematics" in text


def test_write_generated_header(tmp_path):
    from albers_printer import emit_kernels, write_header

    path = tmp_path / "medial_kernels.h"
    write_header(path, emit_kernels())
    text = path.read_text()
    assert "AUTO-GENERATED" in text
    assert "shape_operator_from_GH" in text


def test_cyclide_smooth_compose():
    from cyclide_smooth import compose_cyclide_smooth

    b = compose_cyclide_smooth()
    assert b.E_anchor is not None
    assert b.grad_neighbor.rows == 14
    assert b.M_G_i.rows == 3 and b.M_G_i.cols == 14


def test_cyclide_smooth_G_H_linear_in_Q():
    import sympy as sp
    from cyclide_smooth import compose_cyclide_smooth

    b = compose_cyclide_smooth()
    Qi = b.Qi
    for k in range(14):
        for expr in list(b.G_i) + list(b.H_i.reshape(9, 1)):
            assert sp.diff(expr, Qi[k], 2) == 0


def test_cyclide_smooth_anchor_grad():
    import sympy as sp
    from cyclide_smooth import compose_cyclide_smooth

    b = compose_cyclide_smooth()
    Qi = b.Qi
    Q0 = b.Q0
    wi = sp.Rational(3, 4)
    E = b.E_anchor.subs(b.wi, wi)
    grad_sym = b.grad_anchor.subs(b.wi, wi)
    subs = {Qi[i]: float(i) * 0.1 + 0.5 for i in range(14)}
    subs.update({Q0[i]: float(i) * 0.07 for i in range(14)})
    eps = 1e-6
    for k in range(14):
        up = dict(subs)
        dn = dict(subs)
        up[Qi[k]] += eps
        dn[Qi[k]] -= eps
        num = float(E.subs(up) - E.subs(dn)) / (2 * eps)
        sym = float(grad_sym.subs(subs)[k])
        assert abs(num - sym) < 1e-4


def test_cyclide_smooth_neighbor_grad_matches_numeric():
    from cyclide_smooth import compose_cyclide_smooth

    b = compose_cyclide_smooth()
    Qi = b.Qi
    Qj = b.Qj
    x = b.x_j_in_i
    fixed = {
        b.wi: 0.85,
        b.w_ij: 0.5,
        b.alpha_G: 1.0,
        b.alpha_H: 1.0,
        x[0]: 0.12,
        x[1]: -0.05,
        x[2]: 0.08,
    }
    fixed.update({Qj[i]: 0.25 + 0.015 * i for i in range(14)})
    qi_vals = {Qi[i]: 0.3 + 0.02 * i for i in range(14)}
    E = b.E_neighbor.subs(fixed)
    grad = b.grad_neighbor.subs({**fixed, **qi_vals})
    eps = 1e-6
    for k in range(14):
        up = dict(qi_vals)
        dn = dict(qi_vals)
        up[Qi[k]] += eps
        dn[Qi[k]] -= eps
        num = float(E.subs(up) - E.subs(dn)) / (2 * eps)
        sym = float(grad[k])
        assert abs(num - sym) < 2e-3


def test_cyclide_smooth_hess_neighbor_constant_in_Qi():
    """Quadratic in Qi => Hessian w.r.t. Qi is independent of Qi."""
    from cyclide_smooth import compose_cyclide_smooth

    b = compose_cyclide_smooth()
    Qi = b.Qi
    Qj = b.Qj
    x = b.x_j_in_i
    common = {
        b.wi: 0.85,
        b.w_ij: 0.5,
        b.alpha_G: 1.0,
        b.alpha_H: 1.0,
        x[0]: 0.12,
        x[1]: -0.05,
        x[2]: 0.08,
    }
    common.update({Qj[i]: 0.25 + 0.01 * i for i in range(14)})
    subs_a = {**common, **{Qi[i]: 0.1 * i for i in range(14)}}
    subs_b = {**common, **{Qi[i]: 0.2 * i + 0.3 for i in range(14)}}
    Ha = b.hess_neighbor.subs(subs_a)
    Hb = b.hess_neighbor.subs(subs_b)
    for i in range(14):
        for j in range(14):
            assert float(Ha[i, j]) == pytest.approx(float(Hb[i, j]), rel=0, abs=1e-9)


def test_emit_cyclide_smooth_header():
    from albers_printer import emit_cyclide_smooth

    text = emit_cyclide_smooth()
    assert "cyclide_smooth_anchor_energy" in text
    assert "cyclide_smooth_neighbor_grad" in text
    assert "Q0[0]" in text
    assert "Q[0]" not in text or "Q0[0]" in text
