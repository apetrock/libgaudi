#!/usr/bin/env python3
"""Incompressible phi material on a slab foot — symbolic E_gauge and dE/dq.

Model (see docs/slab_incompressible_phi_material.md):
  f_r(x; Q_r)     rest quadric (fixed)
  u_q(x) = B q    div-free velocity field (columns of B)
  phi_q(x) = x + u_q(x)
  delta D = f_r(phi_q(x))   target 0 on rest iso-surface
  E_gauge = sum_i |grad delta D|^2 + lam_gap (delta D)^2

Usage:
  python3 scripts/sympy/slab_phi_incompressible.py
  python3 scripts/sympy/slab_phi_incompressible.py --fd-check
  python3 scripts/sympy/slab_phi_incompressible.py --nr-steps 8
"""

from __future__ import annotations

import argparse

import numpy as np
import sympy as sp
from sympy import Matrix, symbols

from implicit import gradient, hessian, x, y, z
from quadric import D_expr, Q

# Rest quadric coefficients (fixed).
Qr = symbols("Qr0:10", real=True)

# Div-free control (placeholder: 3 constant solenoidal modes in foot frame).
q = symbols("q0:3", real=True)


def u_linear_solenoidal() -> Matrix:
    """Minimal div-free field: u = q0 * z_hat x r + q1 * x_hat x r + q2 * y_hat x r."""
    r = Matrix([x, y, z])
    ex, ey, ez = Matrix([1, 0, 0]), Matrix([0, 1, 0]), Matrix([0, 0, 1])
    return q[0] * ez.cross(r) + q[1] * ex.cross(r) + q[2] * ey.cross(r)


def f_r_at(pos: Matrix) -> sp.Expr:
    D = D_expr().subs({Q[i]: Qr[i] for i in range(10)})
    pt = {x: pos[0], y: pos[1], z: pos[2]}
    return sp.simplify(D.subs(pt))


def delta_D_at(pos: Matrix, u: Matrix) -> sp.Expr:
    phi = pos + u
    return sp.simplify(f_r_at(phi))


def gauge_energy_density(delta: sp.Expr, lam_gap: sp.Symbol) -> sp.Expr:
    g = gradient(delta)
    return g.dot(g) + lam_gap * delta * delta


def energy_wrt_q(lam_gap: sp.Symbol | None = None) -> tuple[sp.Expr, Matrix]:
    pos = Matrix([x, y, z])
    u = u_linear_solenoidal()
    delta = delta_D_at(pos, u)
    lam = lam_gap if lam_gap is not None else sp.Symbol("lam_gap", real=True, positive=True)
    E_pt = gauge_energy_density(delta, lam)
    E = sp.integrate(sp.integrate(sp.integrate(E_pt, (x, -1, 1)), (y, -1, 1)), (z, -1, 1))
    grad_q = Matrix([sp.diff(E, qi) for qi in q])
    return sp.simplify(E), sp.simplify(grad_q)


def sample_energy_at_point(lam_gap: sp.Symbol) -> tuple[sp.Expr, Matrix]:
    """Energy summed at a single sample (used for FD check on q)."""
    pos = Matrix([x, y, z])
    u = u_linear_solenoidal()
    delta = delta_D_at(pos, u)
    E = gauge_energy_density(delta, lam_gap)
    grad_q = Matrix([sp.diff(E, qi) for qi in q])
    return sp.simplify(E), sp.simplify(grad_q)


def fd_check_bent_q0() -> None:
    """At q=0 on a offset sample, gauge energy and grad should be nonzero."""
    lam = sp.Symbol("lam_gap", real=True, positive=True)
    E, grad_q = sample_energy_at_point(lam)
    subs = {
        x: 0.12,
        y: -0.08,
        z: 0.05,
        lam: 1.0,
        **{Qr[i]: 0.1 * (i + 1) for i in range(10)},
        **{q[i]: 0.0 for i in range(3)},
    }
    e0 = float(E.subs(subs))
    g0 = float(grad_q[0].subs(subs)) ** 2 + float(grad_q[1].subs(subs)) ** 2
    g0 += float(grad_q[2].subs(subs)) ** 2
    print(f"bent q=0: E={e0:.6e} |grad|^2={g0:.6e}")
    if e0 < 1e-20:
        raise SystemExit("expected E > 0 at q=0 on bent sample")
    if g0 < 1e-20:
        raise SystemExit("expected grad != 0 at q=0 on bent sample")


def fd_check_q(eps: float = 1e-6) -> float:
    lam = sp.Symbol("lam_gap", real=True, positive=True)
    E, grad_q = sample_energy_at_point(lam)
    subs = {
        x: 0.12,
        y: -0.08,
        z: 0.05,
        lam: 1.0,
        **{Qr[i]: 0.1 * (i + 1) for i in range(10)},
        **{q[i]: 0.01 * (i + 1) for i in range(3)},
    }
    err = 0.0
    for i, qi in enumerate(q):
        s_plus = dict(subs)
        s_minus = dict(subs)
        s_plus[qi] = float(subs[qi]) + eps
        s_minus[qi] = float(subs[qi]) - eps
        fd = (float(E.subs(s_plus)) - float(E.subs(s_minus))) / (2 * eps)
        sym = float(grad_q[i].subs(subs))
        err = max(err, abs(fd - sym) / max(abs(sym), 1e-12))
    return err


def nr_demo(steps: int = 6, step: float = 0.15) -> None:
    """Simple NR on q for a toy quadric; print energy descent."""
    lam_val = 1.0
    Qr_num = [0.1 * (i + 1) for i in range(10)]
    q_num = np.zeros(3)
    pt = {x: 0.1, y: -0.15, z: 0.02}
    lam = sp.Symbol("lam_gap", real=True, positive=True)

    def eval_E(qv: np.ndarray) -> float:
        subs = {lam: lam_val, **pt}
        subs.update({Qr[i]: Qr_num[i] for i in range(10)})
        subs.update({q[i]: float(qv[i]) for i in range(3)})
        E, _ = sample_energy_at_point(lam)
        return float(E.subs(subs))

    def eval_grad(qv: np.ndarray) -> np.ndarray:
        subs = {lam: lam_val, **pt}
        subs.update({Qr[i]: Qr_num[i] for i in range(10)})
        subs.update({q[i]: float(qv[i]) for i in range(3)})
        _, gq = sample_energy_at_point(lam)
        return np.array([float(gq[i].subs(subs)) for i in range(3)])

    print("NR on q (single sample point):")
    for k in range(steps):
        g = eval_grad(q_num)
        e = eval_E(q_num)
        print(f"  iter {k}: E={e:.6e} |g|={np.linalg.norm(g):.6e}")
        q_num -= step * g
    print(f"  final q = {q_num}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fd-check", action="store_true")
    parser.add_argument("--fd-check-bent-q0", action="store_true")
    parser.add_argument("--nr-steps", type=int, default=0)
    args = parser.parse_args()

    if args.fd_check_bent_q0:
        fd_check_bent_q0()
        return

    if args.fd_check:
        err = fd_check_q()
        print(f"FD check max rel err: {err:.3e}")
        if err > 1e-3:
            raise SystemExit(1)
        return

    if args.nr_steps > 0:
        nr_demo(args.nr_steps)
        return

    lam = sp.Symbol("lam_gap", real=True, positive=True)
    E, grad_q = sample_energy_at_point(lam)
    print("E_gauge sample (sympy):")
    print(E)
    print("grad_q:")
    for i in range(3):
        print(f"  dE/dq{i} = {grad_q[i]}")


if __name__ == "__main__":
    main()
