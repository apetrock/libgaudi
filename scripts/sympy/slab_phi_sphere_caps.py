#!/usr/bin/env python3
"""Discrete incompressible phi on synthetic sphere caps — mirrors C++ wedge basis.

Builds fan-wedge div-nullspace B, fixed quadric Q_r, and runs GD on q for cap-scale
deformations (radial stretch + cap separation). Uses FD gradient for GD; compares
analytic vs FD at the end.

Usage:
  python3 scripts/sympy/slab_phi_sphere_caps.py
  python3 scripts/sympy/slab_phi_sphere_caps.py --eps 0.05 --gd-iters 20
  python3 scripts/sympy/slab_phi_sphere_caps.py --separation 0.08
  python3 scripts/sympy/slab_phi_sphere_caps.py --plot
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field

import numpy as np

from quadric import sphere_coefficients


@dataclass
class Stencil:
    foot: int
    pierce_foot: int
    top_ids: list[int] = field(default_factory=list)
    bot_ids: list[int] = field(default_factory=list)
    top_weights: list[float] = field(default_factory=list)
    bot_weights: list[float] = field(default_factory=list)


def ring_points(lat_deg: float, n: int, hub_z: float) -> np.ndarray:
    """n points on a latitude circle."""
    lat = np.deg2rad(lat_deg)
    z = np.sin(lat) if hub_z > 0 else -np.sin(lat)
    r = np.cos(lat)
    pts = []
    for k in range(n):
        th = 2.0 * np.pi * k / n
        pts.append([r * np.cos(th), r * np.sin(th), z])
    return np.array(pts, dtype=float)


def make_unit_sphere_stencil(n_ring: int = 6, lat_deg: float = 55.0) -> tuple[Stencil, np.ndarray]:
    """Synthetic pierce stencil: north foot + ring, south pierce + ring."""
    foot = np.array([0.0, 0.0, 1.0])
    pierce = np.array([0.0, 0.0, -1.0])
    top_ring = ring_points(lat_deg, n_ring, hub_z=1.0)
    bot_ring = ring_points(lat_deg, n_ring, hub_z=-1.0)

    n = n_ring
    pos = np.zeros((2 * n + 2, 3))
    pos[0] = foot
    pos[1 : 1 + n] = top_ring
    pos[1 + n] = pierce
    pos[2 + n :] = bot_ring

    top_ids = [0] + list(range(1, 1 + n))
    bot_ids = [1 + n] + list(range(2 + n, 2 + n + n))
    w = [1.0 / n] * n
    stencil = Stencil(
        foot=0,
        pierce_foot=1 + n,
        top_ids=top_ids,
        bot_ids=bot_ids,
        top_weights=w,
        bot_weights=w,
    )
    return stencil, pos


def build_fan_wedges(stencil: Stencil) -> list[tuple[int, int, int]]:
    """Return (hub, i0, i1) for each fan triangle (top + bottom)."""
    wedges = []
    n_top = len(stencil.top_ids) - 1
    for k in range(n_top):
        i0 = stencil.top_ids[k + 1]
        i1 = stencil.top_ids[(k + 1) % n_top + 1]
        wedges.append((stencil.foot, i0, i1))
    n_bot = len(stencil.bot_ids) - 1
    for k in range(n_bot):
        i0 = stencil.bot_ids[k + 1]
        i1 = stencil.bot_ids[(k + 1) % n_bot + 1]
        wedges.append((stencil.pierce_foot, i0, i1))
    return wedges


def build_div_free_basis(
    stencil: Stencil, pos_rest: np.ndarray
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """Build B (n_dof x dim_q) and foot_rest; mirrors phi_div_free_basis.hpp."""
    foot_rest = pos_rest[stencil.foot].copy()
    vert_ids = list(stencil.top_ids) + [v for v in stencil.bot_ids if v not in stencil.top_ids]
    vert_idx = {vi: i for i, vi in enumerate(vert_ids)}
    n_dof = 3 * len(vert_ids)

    wedges = build_fan_wedges(stencil)
    rows = []
    for hub, i0, i1 in wedges:
        p0 = pos_rest[hub] - foot_rest
        p1 = pos_rest[i0] - foot_rest
        p2 = pos_rest[i1] - foot_rest
        e01 = p1 - p0
        e02 = p2 - p0
        n = np.cross(e01, e02)
        area2 = np.linalg.norm(n)
        if area2 < 1e-24:
            continue
        grad_l0 = np.cross(e02, n) / area2
        grad_l1 = np.cross(n, e01) / area2
        grad_l2 = -(grad_l0 + grad_l1)
        row = np.zeros(n_dof)
        for vi, gl in ((hub, grad_l0), (i0, grad_l1), (i1, grad_l2)):
            j = vert_idx[vi]
            row[3 * j : 3 * j + 3] = gl
        rows.append(row)

    if not rows:
        B = np.eye(n_dof, 3)
        return B, foot_rest, vert_ids

    D = np.vstack(rows)
    _, s, Vt = np.linalg.svd(D, full_matrices=True)
    tol = max(D.shape) * s[0] * 1e-10 if s.size else 1e-10
    rank = int(np.sum(s > tol))
    null_cols = [Vt[k].T for k in range(rank, Vt.shape[0])]
    if not null_cols:
        B = np.eye(n_dof, 3)
        return B, foot_rest, vert_ids
    B = np.column_stack(null_cols)
    return B, foot_rest, vert_ids


@dataclass
class WedgePatch:
    v0: int
    v1: int
    v2: int
    grad_l0: np.ndarray
    grad_l1: np.ndarray
    grad_l2: np.ndarray


@dataclass
class DivFreeBasis:
    B: np.ndarray
    foot_rest: np.ndarray
    vert_ids: list[int]
    wedges: list[WedgePatch] = field(default_factory=list)

    @property
    def dim_q(self) -> int:
        return self.B.shape[1]

    def u_at(self, vi: int, q: np.ndarray) -> np.ndarray:
        return u_at_vertex(self.B, self.vert_ids, vi, q)

    def max_div(self, q: np.ndarray) -> float:
        m = 0.0
        for wp in self.wedges:
            div = 0.0
            for vi, gl in ((wp.v0, wp.grad_l0), (wp.v1, wp.grad_l1), (wp.v2, wp.grad_l2)):
                div += np.dot(gl, self.u_at(vi, q))
            m = max(m, abs(div))
        return m


def u_at_vertex(B: np.ndarray, vert_ids: list[int], vi: int, q: np.ndarray) -> np.ndarray:
    j = vert_ids.index(vi)
    base = 3 * j
    u = np.zeros(3)
    for m in range(B.shape[1]):
        u += q[m] * B[base : base + 3, m]
    return u


def total_fan_wedge_volume(stencil: Stencil, pos: np.ndarray) -> float:
    """Sum |scalar triple| over top fan wedges extruded along spine."""
    c_top = pos[stencil.foot]
    c_bot = pos[stencil.pierce_foot]
    d = c_bot - c_top
    vol = 0.0
    n_top = len(stencil.top_ids) - 1
    for k in range(n_top):
        i0 = stencil.top_ids[k + 1]
        i1 = stencil.top_ids[(k + 1) % n_top + 1]
        e0 = pos[i0] - c_top
        e1 = pos[i1] - c_top
        vol += abs(np.dot(np.cross(e0, e1), d)) / 6.0
    n_bot = len(stencil.bot_ids) - 1
    for k in range(n_bot):
        i0 = stencil.bot_ids[k + 1]
        i1 = stencil.bot_ids[(k + 1) % n_bot + 1]
        e0 = pos[i0] - c_bot
        e1 = pos[i1] - c_bot
        vol += abs(np.dot(np.cross(e0, e1), -d)) / 6.0
    return vol


def displacement_hoop_fraction(pos: np.ndarray, d: np.ndarray, vi: int) -> float:
    """Fraction of displacement toward pierce axis (hoop inward in xy)."""
    p = pos[vi]
    hoop = np.array([-p[0], -p[1], 0.0])
    nh = np.linalg.norm(hoop)
    nd = np.linalg.norm(d)
    if nh < 1e-12 or nd < 1e-12:
        return 0.0
    return max(0.0, np.dot(d, hoop / nh) / nd)


def max_div_on_wedges(
    stencil: Stencil,
    pos_rest: np.ndarray,
    B: np.ndarray,
    vert_ids: list[int],
    q: np.ndarray,
) -> float:
    foot_rest = pos_rest[stencil.foot]
    vert_idx = {vi: i for i, vi in enumerate(vert_ids)}
    n_dof = 3 * len(vert_ids)
    u_dof = B @ q
    m = 0.0
    for hub, i0, i1 in build_fan_wedges(stencil):
        p0 = pos_rest[hub] - foot_rest
        p1 = pos_rest[i0] - foot_rest
        p2 = pos_rest[i1] - foot_rest
        e01 = p1 - p0
        e02 = p2 - p0
        n = np.cross(e01, e02)
        area2 = np.linalg.norm(n)
        if area2 < 1e-24:
            continue
        grad_l0 = np.cross(e02, n) / area2
        grad_l1 = np.cross(n, e01) / area2
        grad_l2 = -(grad_l0 + grad_l1)
        div = 0.0
        for vi, gl in ((hub, grad_l0), (i0, grad_l1), (i1, grad_l2)):
            j = vert_idx[vi]
            u_v = u_dof[3 * j : 3 * j + 3]
            div += np.dot(gl, u_v)
        m = max(m, abs(div))
    return m


def basis_field(B: np.ndarray, vert_ids: list[int], vi: int, mode: int) -> np.ndarray:
    j = vert_ids.index(vi)
    base = 3 * j
    return B[base : base + 3, mode].copy()


def quadric_design_row(local: np.ndarray) -> np.ndarray:
    x, y, z = local
    return np.array(
        [x * x, y * y, z * z, x * y, x * z, y * z, 2 * x, 2 * y, 2 * z, 1.0]
    )


def fit_quadric_on_stencil(
    stencil: Stencil, pos: np.ndarray, foot_normal_weight: float = 1.0
) -> np.ndarray:
    """Foot-local LS quadric (spoke samples f=0 + foot grad·n=1), mirrors C++ fit."""
    foot = pos[stencil.foot]
    rows: list[np.ndarray] = []
    wsum = 0.0
    for k in range(len(stencil.top_weights)):
        vi = stencil.top_ids[k + 1]
        w = stencil.top_weights[k]
        if w <= 0:
            continue
        rows.append(np.sqrt(w) * quadric_design_row(pos[vi] - foot))
        wsum += w
    for k in range(len(stencil.bot_weights)):
        vi = stencil.bot_ids[k + 1]
        w = stencil.bot_weights[k]
        if w <= 0:
            continue
        rows.append(np.sqrt(w) * quadric_design_row(pos[vi] - foot))
        wsum += w
    if foot_normal_weight > 0 and wsum > 0:
        # grad f at origin = (2Q6, 2Q7, 2Q8); pin ∇f·e_z ≈ 1 at north cap foot
        n_top = np.array([0.0, 0.0, 1.0])
        row = np.zeros(10)
        row[6:9] = 2.0 * n_top
        rows.append(np.sqrt(foot_normal_weight * wsum) * row)
    A = np.vstack(rows)
    Q, _, _, _ = np.linalg.lstsq(A, np.zeros(A.shape[0]), rcond=None)
    return Q


def eval_quadric(Q: np.ndarray, x: np.ndarray) -> float:
    return (
        Q[0] * x[0] ** 2
        + Q[1] * x[1] ** 2
        + Q[2] * x[2] ** 2
        + Q[3] * x[0] * x[1]
        + Q[4] * x[0] * x[2]
        + Q[5] * x[1] * x[2]
        + 2 * Q[6] * x[0]
        + 2 * Q[7] * x[1]
        + 2 * Q[8] * x[2]
        + Q[9]
    )


def eval_quadric_grad(Q: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.array(
        [
            2 * Q[0] * x[0] + 2 * Q[3] * x[1] + 2 * Q[4] * x[2] + 2 * Q[6],
            2 * Q[3] * x[0] + 2 * Q[1] * x[1] + 2 * Q[5] * x[2] + 2 * Q[7],
            2 * Q[4] * x[0] + 2 * Q[5] * x[1] + 2 * Q[2] * x[2] + 2 * Q[8],
        ]
    )


def quadric_hess(Q: np.ndarray) -> np.ndarray:
    return 2.0 * np.array(
        [[Q[0], Q[3], Q[4]], [Q[3], Q[1], Q[5]], [Q[4], Q[5], Q[2]]]
    )


def ring_vert_ids(stencil: Stencil) -> list[int]:
    ids = list(stencil.top_ids) + [v for v in stencil.bot_ids if v not in stencil.top_ids]
    return ids


def gauge_sample_vert_ids(stencil: Stencil) -> list[int]:
    """Spoke ring verts only (exclude foot / pierce hubs)."""
    out = [vi for vi in stencil.top_ids[1:]]
    for vi in stencil.bot_ids[1:]:
        if vi not in out:
            out.append(vi)
    return out


def gauge_energy(
    Q: np.ndarray,
    pos_bent: np.ndarray,
    stencil: Stencil,
    B: np.ndarray,
    vert_ids: list[int],
    foot_rest: np.ndarray,
    q: np.ndarray,
    lam_gap: float = 1.0,
    mu: float = 1.0,
) -> float:
    E = 0.0
    for vi in gauge_sample_vert_ids(stencil):
        local = pos_bent[vi] - foot_rest
        u = u_at_vertex(B, vert_ids, vi, q)
        x = local + u
        f = eval_quadric(Q, x)
        if not np.isfinite(f):
            continue
        E += mu * lam_gap * f * f
    return E


def gauge_grad_analytic(
    Q: np.ndarray,
    pos_bent: np.ndarray,
    stencil: Stencil,
    B: np.ndarray,
    vert_ids: list[int],
    foot_rest: np.ndarray,
    q: np.ndarray,
    lam_gap: float = 1.0,
    mu: float = 1.0,
) -> np.ndarray:
    dim = B.shape[1]
    grad = np.zeros(dim)
    for vi in gauge_sample_vert_ids(stencil):
        local = pos_bent[vi] - foot_rest
        u = u_at_vertex(B, vert_ids, vi, q)
        x = local + u
        f = eval_quadric(Q, x)
        g = eval_quadric_grad(Q, x)
        if not np.isfinite(f):
            continue
        for j in range(dim):
            bj = basis_field(B, vert_ids, vi, j)
            dF = np.dot(g, bj)
            grad[j] += mu * lam_gap * 2.0 * f * dF
    return grad


def gauge_grad_fd(
    Q: np.ndarray,
    pos_bent: np.ndarray,
    stencil: Stencil,
    B: np.ndarray,
    vert_ids: list[int],
    foot_rest: np.ndarray,
    q: np.ndarray,
    lam_gap: float = 1.0,
    mu: float = 1.0,
    eps: float = 1e-6,
) -> np.ndarray:
    dim = len(q)
    grad = np.zeros(dim)
    for i in range(dim):
        qp = q.copy()
        qm = q.copy()
        qp[i] += eps
        qm[i] -= eps
        ep = gauge_energy(Q, pos_bent, stencil, B, vert_ids, foot_rest, qp, lam_gap, mu)
        em = gauge_energy(Q, pos_bent, stencil, B, vert_ids, foot_rest, qm, lam_gap, mu)
        grad[i] = (ep - em) / (2.0 * eps)
    return grad


def fan_wedge_taus_and_weights(
    stencil: Stencil, pos: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Oriented fan-wedge scalars tau (matches trivector3 coeff in C++)."""
    c_top = pos[stencil.foot]
    c_bot = pos[stencil.pierce_foot]
    d = c_bot - c_top
    taus: list[float] = []
    weights: list[float] = []
    n_top = len(stencil.top_ids) - 1
    for k in range(n_top):
        i0 = stencil.top_ids[k + 1]
        i1 = stencil.top_ids[(k + 1) % n_top + 1]
        e0 = pos[i0] - c_top
        e1 = pos[i1] - c_top
        taus.append(float(np.dot(np.cross(e0, e1), d)))
        weights.append(stencil.top_weights[k])
    n_bot = len(stencil.bot_ids) - 1
    for k in range(n_bot):
        i0 = stencil.bot_ids[k + 1]
        i1 = stencil.bot_ids[(k + 1) % n_bot + 1]
        e0 = pos[i0] - c_bot
        e1 = pos[i1] - c_bot
        taus.append(float(np.dot(np.cross(e0, e1), -d)))
        weights.append(stencil.bot_weights[k])
    return np.array(taus), np.array(weights)


def wedge_strain_energy(
    stencil: Stencil,
    pos_bent: np.ndarray,
    B: np.ndarray,
    vert_ids: list[int],
    q: np.ndarray,
    tau_rest: np.ndarray,
    weights: np.ndarray,
) -> float:
    pos = apply_phi(pos_bent, stencil, B, vert_ids, q)
    taus, w = fan_wedge_taus_and_weights(stencil, pos)
    diff = taus - tau_rest
    return float(np.sum(w * diff * diff))


def combined_energy(
    Q: np.ndarray,
    pos_bent: np.ndarray,
    stencil: Stencil,
    B: np.ndarray,
    vert_ids: list[int],
    foot_rest: np.ndarray,
    q: np.ndarray,
    tau_rest: np.ndarray,
    wedge_weights: np.ndarray,
    lam_gap: float = 1.0,
    mu: float = 1.0,
    nu_wedge: float = 0.0,
) -> float:
    e = gauge_energy(Q, pos_bent, stencil, B, vert_ids, foot_rest, q, lam_gap, mu)
    if nu_wedge > 0.0:
        e += nu_wedge * wedge_strain_energy(
            stencil, pos_bent, B, vert_ids, q, tau_rest, wedge_weights
        )
    return e


def combined_grad_fd(
    Q: np.ndarray,
    pos_bent: np.ndarray,
    stencil: Stencil,
    B: np.ndarray,
    vert_ids: list[int],
    foot_rest: np.ndarray,
    q: np.ndarray,
    tau_rest: np.ndarray,
    wedge_weights: np.ndarray,
    lam_gap: float = 1.0,
    mu: float = 1.0,
    nu_wedge: float = 0.0,
    eps: float = 1e-6,
) -> np.ndarray:
    dim = len(q)
    grad = np.zeros(dim)
    for i in range(dim):
        qp = q.copy()
        qm = q.copy()
        qp[i] += eps
        qm[i] -= eps
        ep = combined_energy(
            Q, pos_bent, stencil, B, vert_ids, foot_rest, qp,
            tau_rest, wedge_weights, lam_gap, mu, nu_wedge,
        )
        em = combined_energy(
            Q, pos_bent, stencil, B, vert_ids, foot_rest, qm,
            tau_rest, wedge_weights, lam_gap, mu, nu_wedge,
        )
        grad[i] = (ep - em) / (2.0 * eps)
    return grad


def resolve_gd(
    Q: np.ndarray,
    pos_bent: np.ndarray,
    stencil: Stencil,
    B: np.ndarray,
    vert_ids: list[int],
    foot_rest: np.ndarray,
    q0: np.ndarray,
    lam_gap: float,
    iters: int,
    tau_rest: np.ndarray | None = None,
    wedge_weights: np.ndarray | None = None,
    step_scale: float = 0.25,
    mu: float = 1.0,
    nu_wedge: float = 0.0,
    use_fd: bool = True,
) -> np.ndarray:
    if tau_rest is None:
        tau_rest = np.zeros(0)
    if wedge_weights is None:
        wedge_weights = np.zeros(0)
    q = q0.copy()
    e_rest = combined_energy(
        Q, pos_bent, stencil, B, vert_ids, foot_rest, np.zeros_like(q0),
        tau_rest, wedge_weights, lam_gap, mu, nu_wedge,
    )
    if e_rest < 1e-20:
        return np.zeros_like(q0)
    for _ in range(iters):
        if use_fd or nu_wedge > 0.0:
            grad = combined_grad_fd(
                Q, pos_bent, stencil, B, vert_ids, foot_rest, q,
                tau_rest, wedge_weights, lam_gap, mu, nu_wedge,
            )
        else:
            grad = gauge_grad_analytic(
                Q, pos_bent, stencil, B, vert_ids, foot_rest, q, lam_gap, mu
            )
        g2 = np.dot(grad, grad)
        if g2 < 1e-32:
            break
        e0 = combined_energy(
            Q, pos_bent, stencil, B, vert_ids, foot_rest, q,
            tau_rest, wedge_weights, lam_gap, mu, nu_wedge,
        )
        alpha = step_scale / np.sqrt(g2)
        improved = False
        for _ in range(8):
            trial = q - alpha * grad
            e1 = combined_energy(
                Q, pos_bent, stencil, B, vert_ids, foot_rest, trial,
                tau_rest, wedge_weights, lam_gap, mu, nu_wedge,
            )
            if np.isfinite(e1) and e1 < e0:
                q = trial
                improved = True
                break
            alpha *= 0.5
        if not improved:
            break
    return q


def apply_phi(
    pos_bent: np.ndarray,
    stencil: Stencil,
    B: np.ndarray,
    vert_ids: list[int],
    q: np.ndarray,
) -> np.ndarray:
    out = pos_bent.copy()
    for vi in ring_vert_ids(stencil):
        out[vi] = pos_bent[vi] + u_at_vertex(B, vert_ids, vi, q)
    return out


def radial_stretch(pos: np.ndarray, stencil: Stencil, eps: float) -> np.ndarray:
    out = pos.copy()
    c_top = pos[stencil.foot]
    for vi in stencil.top_ids[1:]:
        out[vi] = c_top + (1.0 + eps) * (pos[vi] - c_top)
    c_bot = pos[stencil.pierce_foot]
    for vi in stencil.bot_ids[1:]:
        out[vi] = c_bot + (1.0 + eps) * (pos[vi] - c_bot)
    return out


def cap_separation(pos: np.ndarray, stencil: Stencil, delta: float) -> np.ndarray:
    """Move foot up and pierce down along z by delta/2 each; ring follows rigidly."""
    out = pos.copy()
    dz_top = 0.5 * delta
    dz_bot = -0.5 * delta
    c_top_rest = pos[stencil.foot]
    c_bot_rest = pos[stencil.pierce_foot]
    out[stencil.foot] = c_top_rest + np.array([0, 0, dz_top])
    for vi in stencil.top_ids[1:]:
        out[vi] = pos[vi] + np.array([0, 0, dz_top])
    out[stencil.pierce_foot] = c_bot_rest + np.array([0, 0, dz_bot])
    for vi in stencil.bot_ids[1:]:
        out[vi] = pos[vi] + np.array([0, 0, dz_bot])
    return out


def max_field_magnitude(
    Q: np.ndarray,
    pos: np.ndarray,
    stencil: Stencil,
    B: np.ndarray,
    vert_ids: list[int],
    foot_rest: np.ndarray,
    q: np.ndarray,
) -> float:
    m = 0.0
    for vi in gauge_sample_vert_ids(stencil):
        local = pos[vi] - foot_rest
        u = u_at_vertex(B, vert_ids, vi, q)
        f = abs(eval_quadric(Q, local + u))
        m = max(m, f)
    return m


def displacement_stats(
    pos_before: np.ndarray,
    pos_after: np.ndarray,
    stencil: Stencil,
    ring_only: list[int] | None = None,
) -> dict:
    ids = ring_only if ring_only is not None else ring_vert_ids(stencil)
    disp = [np.linalg.norm(pos_after[vi] - pos_before[vi]) for vi in ids]
    return {
        "mean": float(np.mean(disp)),
        "max": float(np.max(disp)),
        "per_vert": disp,
    }


def mean_hoop_inward(pos: np.ndarray, pos_before: np.ndarray, pos_after: np.ndarray,
                     ring_ids: list[int]) -> float:
    fracs = []
    for vi in ring_ids:
        d = pos_after[vi] - pos_before[vi]
        fracs.append(displacement_hoop_fraction(pos, d, vi))
    return float(np.mean(fracs)) if fracs else 0.0


def run_scenario(
    name: str,
    pos_rest: np.ndarray,
    pos_bent: np.ndarray,
    stencil: Stencil,
    Q: np.ndarray,
    B: np.ndarray,
    vert_ids: list[int],
    foot_rest: np.ndarray,
    tau_rest: np.ndarray,
    wedge_weights: np.ndarray,
    gd_iters: int,
    nu_wedge: float = 0.0,
    mu_gauge: float = 1.0,
) -> dict:
    q0 = np.zeros(B.shape[1])
    e0 = combined_energy(
        Q, pos_bent, stencil, B, vert_ids, foot_rest, q0,
        tau_rest, wedge_weights, 1.0, mu_gauge, nu_wedge,
    )
    q1 = resolve_gd(
        Q, pos_bent, stencil, B, vert_ids, foot_rest, q0, 1.0, gd_iters,
        tau_rest, wedge_weights, 0.25, mu_gauge, nu_wedge, use_fd=True,
    )
    e1 = combined_energy(
        Q, pos_bent, stencil, B, vert_ids, foot_rest, q1,
        tau_rest, wedge_weights, 1.0, mu_gauge, nu_wedge,
    )
    pos_star = apply_phi(pos_bent, stencil, B, vert_ids, q1)

    stats = displacement_stats(pos_bent, pos_star, stencil)
    f0 = max_field_magnitude(Q, pos_bent, stencil, B, vert_ids, foot_rest, q0)
    f1 = max_field_magnitude(Q, pos_bent, stencil, B, vert_ids, foot_rest, q1)
    max_div = max_div_on_wedges(stencil, pos_rest, B, vert_ids, q1)
    taus_bent, _ = fan_wedge_taus_and_weights(stencil, pos_bent)
    taus_star, _ = fan_wedge_taus_and_weights(stencil, pos_star)
    tau_err_bent = float(np.sum(wedge_weights * (taus_bent - tau_rest) ** 2))
    tau_err_star = float(np.sum(wedge_weights * (taus_star - tau_rest) ** 2))
    v_rest = total_fan_wedge_volume(stencil, pos_rest)
    v_bent = total_fan_wedge_volume(stencil, pos_bent)
    v_star = total_fan_wedge_volume(stencil, pos_star)

    top_ids = [v for v in stencil.top_ids[1:]]
    hoop_mean = mean_hoop_inward(pos_bent, pos_bent, pos_star, top_ids)

    print(f"\n=== {name} (mu={mu_gauge}, nu_wedge={nu_wedge}) ===")
    print(f"  dim_q={B.shape[1]}  E0={e0:.6e}  E1={e1:.6e}  dE={e0 - e1:.6e}")
    print(f"  max|f| q=0: {f0:.6e}  q*: {f1:.6e}")
    print(f"  max_div(q*)={max_div:.4e}")
    print(f"  tau_err rest→bent={tau_err_bent:.6e}  rest→star={tau_err_star:.6e}")
    print(f"  |V| rest/bent/star={v_rest:.4f}/{v_bent:.4f}/{v_star:.4f}")
    print(f"  top ring hoop_inward mean={hoop_mean:.3f}  disp max={stats['max']:.4f}")

    ok_vol = tau_err_star < tau_err_bent * 1.01 + 1e-9
    ok_gauge = f1 < f0 * 0.95
    ok_div = max_div < 1e-6
    ok_hoop = hoop_mean > 0.15
    ok_tau_abs = tau_err_star < 0.05
    ok_f_abs = f1 < 0.5
    status = (
        "PASS"
        if (ok_vol and ok_gauge and ok_div and ok_tau_abs and ok_f_abs)
        else "CHECK"
    )
    print(f"  incompressible check: {status}"
          f" (vol_vs_bent={ok_vol} gauge={ok_gauge} div={ok_div}"
          f" tau*<0.05={ok_tau_abs} |f|*<0.5={ok_f_abs} hoop={ok_hoop:.0f})")

    return {
        "hoop_mean": hoop_mean,
        "tau_err_bent": tau_err_bent,
        "tau_err_star": tau_err_star,
        "max_div": max_div,
        "f0": f0,
        "f1": f1,
        "v_star": v_star,
        "status": status,
    }


def run_at_rest(
    stencil: Stencil,
    pos_rest: np.ndarray,
    Q: np.ndarray,
    B: np.ndarray,
    vert_ids: list[int],
    foot_rest: np.ndarray,
    tau_rest: np.ndarray,
    wedge_weights: np.ndarray,
    gd_iters: int,
    nu_wedge: float,
    mu_gauge: float = 1.0,
) -> bool:
    """At rest, q=0: zero energy, GD returns ~0, no ring displacement."""
    q0 = np.zeros(B.shape[1])
    e0 = combined_energy(
        Q, pos_rest, stencil, B, vert_ids, foot_rest, q0,
        tau_rest, wedge_weights, 1.0, mu_gauge, nu_wedge,
    )
    ga = gauge_grad_analytic(
        Q, pos_rest, stencil, B, vert_ids, foot_rest, q0, 1.0, mu_gauge,
    )
    q_warm = np.array([0.02 * (i + 1) for i in range(B.shape[1])])
    q1 = resolve_gd(
        Q, pos_rest, stencil, B, vert_ids, foot_rest, q_warm, 1.0, gd_iters,
        tau_rest, wedge_weights, 0.25, mu_gauge, nu_wedge, use_fd=True,
    )
    e1 = combined_energy(
        Q, pos_rest, stencil, B, vert_ids, foot_rest, q1,
        tau_rest, wedge_weights, 1.0, mu_gauge, nu_wedge,
    )
    pos_star = apply_phi(pos_rest, stencil, B, vert_ids, q1)
    max_disp = max(
        np.linalg.norm(pos_star[vi] - pos_rest[vi]) for vi in gauge_sample_vert_ids(stencil)
    )
    ok = e0 < 1e-20 and e1 < 1e-16 and np.linalg.norm(q1) < 1e-4 and max_disp < 1e-6
    print(f"\n=== at rest (mu={mu_gauge}, nu={nu_wedge}) ===")
    print(f"  E(q=0)={e0:.4e}  |grad|={np.linalg.norm(ga):.4e}  E(q*)={e1:.4e}")
    print(f"  |q*|={np.linalg.norm(q1):.4e}  max_disp={max_disp:.4e}  PASS={ok}")
    return ok


def run_two_phase(
    pos_rest: np.ndarray,
    pos_bent: np.ndarray,
    stencil: Stencil,
    Q: np.ndarray,
    B: np.ndarray,
    vert_ids: list[int],
    foot_rest: np.ndarray,
    tau_rest: np.ndarray,
    wedge_weights: np.ndarray,
    gd_iters: int,
    nu_wedge: float,
    mu_gauge: float,
) -> dict:
    """Volume memory first (mu=0), then joint polish."""
    q0 = np.zeros(B.shape[1])
    q_vol = resolve_gd(
        Q, pos_bent, stencil, B, vert_ids, foot_rest, q0, 1.0, gd_iters,
        tau_rest, wedge_weights, 0.25, 0.0, nu_wedge, use_fd=True,
    )
    q1 = resolve_gd(
        Q, pos_bent, stencil, B, vert_ids, foot_rest, q_vol, 1.0, gd_iters,
        tau_rest, wedge_weights, 0.25, mu_gauge, nu_wedge, use_fd=True,
    )
    pos_star = apply_phi(pos_bent, stencil, B, vert_ids, q1)
    taus_star, _ = fan_wedge_taus_and_weights(stencil, pos_star)
    tau_err_star = float(np.sum(wedge_weights * (taus_star - tau_rest) ** 2))
    f1 = max_field_magnitude(Q, pos_bent, stencil, B, vert_ids, foot_rest, q1)
    v_star = total_fan_wedge_volume(stencil, pos_star)
    top_ids = [v for v in stencil.top_ids[1:]]
    hoop_mean = mean_hoop_inward(pos_bent, pos_bent, pos_star, top_ids)
    print(f"\n=== two-phase (vol nu={nu_wedge}, then mu={mu_gauge}+nu) ===")
    print(f"  tau_err*={tau_err_star:.4e}  max|f|*={f1:.4e}  |V|*={v_star:.4f}"
          f"  hoop={hoop_mean:.3f}")
    return {"tau_err_star": tau_err_star, "f1": f1, "v_star": v_star, "hoop_mean": hoop_mean}


def maybe_plot(
    pos_rest: np.ndarray,
    pos_bent: np.ndarray,
    pos_star: np.ndarray,
    stencil: Stencil,
) -> None:
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    except ImportError:
        print("(matplotlib not installed — skip --plot)")
        return

    fig = plt.figure(figsize=(10, 4))
    for ax_idx, (pos, title) in enumerate(
        [(pos_rest, "rest"), (pos_bent, "bent"), (pos_star, "after phi GD")]
    ):
        ax = fig.add_subplot(1, 3, ax_idx + 1, projection="3d")
        ids = ring_vert_ids(stencil)
        P = pos[ids]
        ax.scatter(P[:, 0], P[:, 1], P[:, 2], s=40)
        ax.scatter(
            pos[stencil.foot, 0],
            pos[stencil.foot, 1],
            pos[stencil.foot, 2],
            c="r",
            s=80,
            label="foot",
        )
        ax.scatter(
            pos[stencil.pierce_foot, 0],
            pos[stencil.pierce_foot, 1],
            pos[stencil.pierce_foot, 2],
            c="b",
            s=80,
            label="pierce",
        )
        ax.set_title(title)
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_zlim(-1.5, 1.5)
    plt.tight_layout()
    plt.savefig("slab_phi_sphere_caps.png", dpi=120)
    print("Wrote slab_phi_sphere_caps.png")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Incompressible phi on sphere caps — iterate E_gauge + nu*E_wedge"
    )
    parser.add_argument("--eps", type=float, default=0.05, help="radial stretch factor-1")
    parser.add_argument("--separation", type=float, default=0.08, help="cap hub separation")
    parser.add_argument("--nu-wedge", type=float, default=1.0, help="wedge volume memory weight")
    parser.add_argument("--mu-gauge", type=float, default=1.0, help="quadric gauge weight")
    parser.add_argument("--sweep-nu", action="store_true", help="sweep nu on cap separation")
    parser.add_argument("--sweep-mu", action="store_true", help="sweep mu at fixed nu-wedge")
    parser.add_argument("--two-phase", action="store_true", help="volume-first then joint GD")
    parser.add_argument("--at-rest-only", action="store_true", help="only run at-rest identity check")
    parser.add_argument("--gd-iters", type=int, default=30)
    parser.add_argument("--n-ring", type=int, default=6)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    stencil, pos_rest = make_unit_sphere_stencil(n_ring=args.n_ring)
    B, foot_rest, vert_ids = build_div_free_basis(stencil, pos_rest)
    Q = fit_quadric_on_stencil(stencil, pos_rest)
    tau_rest, wedge_weights = fan_wedge_taus_and_weights(stencil, pos_rest)

    print("Incompressible quadric phi — sphere cap diagnostic")
    print(f"  n_ring={args.n_ring}  dim_q={B.shape[1]}  foot_rest={foot_rest}")
    print(f"  E = mu*E_gauge + nu*E_wedge   (Q_r from foot-local LS fit)")

    ok_rest = run_at_rest(
        stencil, pos_rest, Q, B, vert_ids, foot_rest, tau_rest, wedge_weights,
        args.gd_iters, args.nu_wedge, args.mu_gauge,
    )
    if not ok_rest:
        raise SystemExit("at-rest identity failed")

    if args.at_rest_only:
        return

    pos_sep = cap_separation(pos_rest, stencil, args.separation)

    if args.two_phase:
        run_two_phase(
            pos_rest, pos_sep, stencil, Q, B, vert_ids, foot_rest,
            tau_rest, wedge_weights, args.gd_iters, args.nu_wedge, args.mu_gauge,
        )
        return

    if args.sweep_mu:
        print(f"\n--- mu sweep (cap separation, nu={args.nu_wedge}) ---")
        print(f"{'mu':>8} {'hoop':>8} {'tau_err*':>12} {'max|f|*':>10} {'|V|*':>8} status")
        for mu in (0.0, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0):
            r = run_scenario(
                f"cap sep delta={args.separation}",
                pos_rest, pos_sep, stencil, Q, B, vert_ids, foot_rest,
                tau_rest, wedge_weights, args.gd_iters,
                nu_wedge=args.nu_wedge, mu_gauge=mu,
            )
            print(
                f"{mu:8.3g} {r['hoop_mean']:8.3f} {r['tau_err_star']:12.4e}"
                f" {r['f1']:10.4e} {r['v_star']:8.4f} {r['status']}"
            )
        return

    if args.sweep_nu:
        print("\n--- nu sweep (cap separation) ---")
        print(f"{'nu':>8} {'hoop':>8} {'tau_err*':>12} {'max|f|*':>10} {'|V|*':>8} status")
        for nu in (0.0, 0.1, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0):
            r = run_scenario(
                f"cap sep delta={args.separation}",
                pos_rest, pos_sep, stencil, Q, B, vert_ids, foot_rest,
                tau_rest, wedge_weights, args.gd_iters, nu_wedge=nu,
                mu_gauge=args.mu_gauge,
            )
            print(
                f"{nu:8.3g} {r['hoop_mean']:8.3f} {r['tau_err_star']:12.4e}"
                f" {r['f1']:10.4e} {r['v_star']:8.4f} {r['status']}"
            )
        return

    run_scenario(
        f"cap separation delta={args.separation}",
        pos_rest, pos_sep, stencil, Q, B, vert_ids, foot_rest,
        tau_rest, wedge_weights, args.gd_iters,
        nu_wedge=args.nu_wedge, mu_gauge=args.mu_gauge,
    )

    pos_radial = radial_stretch(pos_rest, stencil, args.eps)
    run_scenario(
        f"radial stretch eps={args.eps}",
        pos_rest, pos_radial, stencil, Q, B, vert_ids, foot_rest,
        tau_rest, wedge_weights, args.gd_iters,
        nu_wedge=args.nu_wedge, mu_gauge=args.mu_gauge,
    )

    if args.plot:
        q0 = np.zeros(B.shape[1])
        q1 = resolve_gd(
            Q, pos_sep, stencil, B, vert_ids, foot_rest, q0, 1.0, args.gd_iters,
            tau_rest, wedge_weights, 0.25, 1.0, args.nu_wedge, use_fd=True,
        )
        pos_star = apply_phi(pos_sep, stencil, B, vert_ids, q1)
        maybe_plot(pos_rest, pos_sep, pos_star, stencil)


if __name__ == "__main__":
    main()
