"""Self-advection of curvature (MATH.md 19.21): a Burgers equation on a torus.

A perturbation `q` of the mean curvature rides along the ring direction with
speed `c q`: `q_t + c q dq/ds = 0` on each ring circle (`s` arclength). From
`q0 = A sin(u)` a shock forms at `t* = rho / (c A)` (`rho` the ring radius),
so the outer equator creases later than the inner one: non-uniform creasing
from a uniform start.

Per scheme the steepest ring-direction gradient of `q` over time is reported
against the exact pre-shock value `(A / rho) / (1 - c A t / rho)`, with the
range of `q` (overshoot) and, under `pullback`, the anisotropy `|z|` the
compressive flow generates from an isotropic start.
"""

import time

import numpy as np

import cr_ops as ops
import curvature_advection as ca


def setup(R=2.0, r=0.7, nu=64, nv=16):
    X, F, _ = ops.torus(R, r, nu, nv)
    m = ops.build_mesh(F, len(X))
    adv = ca.SurfaceAdvector(m, X)
    u = np.arctan2(X[:, 1], X[:, 0])
    e_ring = np.column_stack([-np.sin(u), np.cos(u), np.zeros_like(u)])
    rho = np.linalg.norm(X[:, :2], axis=1)
    d = X[m.E[:, 1]] - X[m.E[:, 0]]
    L = np.linalg.norm(d, axis=1)
    mid = 0.5 * (X[m.E[:, 0]] + X[m.E[:, 1]])
    um = np.arctan2(mid[:, 1], mid[:, 0])
    em = np.column_stack([-np.sin(um), np.cos(um), np.zeros_like(um)])
    ring_edges = np.abs(np.sum(d * em, 1)) / L > 0.9
    return X, m, adv, u, e_ring, rho, L, ring_edges


def steepest(q, m, L, ring_edges, rho, band):
    g = np.abs(q[m.E[:, 1]] - q[m.E[:, 0]]) / L
    rm = 0.5 * (rho[m.E[:, 0]] + rho[m.E[:, 1]])
    sel = ring_edges & band(rm)
    return float(g[sel].max())


def run(scheme, transport="parallel", A=1.0, c=1.0, T=2.4, steps=24, geom=None, blowup=20.0):
    """Returns the sampled rows and the time `q` first exceeds `blowup A`
    (None if it never does)."""
    X, m, adv, u, e_ring, rho, L, ring_edges = geom or setup()
    q = A * np.sin(u)
    z = np.zeros(m.n, complex)
    dt = T / steps
    r_in, r_out = rho.min(), rho.max()
    inner = lambda rm: rm < r_in + 1e-6 + 0.05
    outer = lambda rm: rm > r_out - 1e-6 - 0.05
    rows = []
    for s in range(1, steps + 1):
        V = c * q[:, None] * e_ring
        with np.errstate(all="ignore"):
            q, z = adv.step(q, z, V, dt, scheme, transport)
        if not np.all(np.isfinite(q)) or np.abs(q).max() > blowup * A:
            return rows, s * dt
        if s % (steps // 6) == 0:
            t = s * dt
            ex = lambda rr: (A / rr) / (1 - c * A * t / rr) if c * A * t < rr else np.inf
            rows.append((t, steepest(q, m, L, ring_edges, rho, inner), ex(r_in),
                         steepest(q, m, L, ring_edges, rho, outer), ex(r_out),
                         float(q.min()), float(q.max()), float(np.abs(z).max())))
    return rows, None


if __name__ == "__main__":
    t0 = time.time()
    geom = setup()
    r_in, r_out = geom[5].min(), geom[5].max()
    print(f"Burgers on a torus ring field: q0 = A sin(u), shock at t* = rho / A: inner {r_in:.2f}, outer {r_out:.2f} (A = 1)")
    print("1. Scalar Burgers (parallel transport), steepest ring gradient against exact")
    for scheme in ("first", "bfecc", "maccormack"):
        rows, _ = run(scheme, geom=geom)
        print(f"  {scheme:10s}  t   |  inner grad (exact)  |  outer grad (exact)  |  q range", flush=True)
        for t, gi, ei, go, eo, lo, hi, zm in rows:
            print(f"     {t:4.2f} | {gi:6.2f} ({ei:6.2f})      | {go:6.2f} ({eo:6.2f})      | [{lo:+.3f}, {hi:+.3f}]", flush=True)
    print("2. Tensor Burgers (pullback): H rides its own compressive flow, J^T A J amplifies it")
    for A in (0.25, 0.5, 1.0):
        rows, tb = run("first", "pullback", A=A, T=2.4 / A, steps=48, geom=geom)
        last = rows[-1] if rows else None
        zm = f"{last[7]:.3f}" if last else "-"
        print(f"  A {A:4.2f}: inner shock t* {r_in / A:5.2f}, blow-up (|q| > 20 A) at "
              f"{'never' if tb is None else f'{tb:5.2f}'}   last sampled max |z| {zm}", flush=True)
    print(f"({time.time() - t0:.0f}s)")
