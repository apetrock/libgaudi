"""Checks for `curvature_advection.py` (MATH.md 19.21).

1. Rigid rotation of a sphere carrying a curvature bump (an `H` bump and a
   `z` pattern): after a full turn, error against the start per scheme and
   transport, at two resolutions (order and dissipation).
2. Torus, a bump transported along a loxodrome: where the peak lands against
   the exact path, and how much of its height is kept, per scheme.
3. A simple shear on a flat sheet: `pullback` against the exact `J^T A J`,
   `parallel` against the same (it does not deform).
4. The MacCormack limiter: no value beyond the start's range.
"""

import sys
import time

import numpy as np

import closed_targets as ct
import cr_ops as ops
import curvature_advection as ca
import dupin_poisson as dp


def sphere(k):
    X, F = ct.cube(k)
    X = X / np.linalg.norm(X, axis=1)[:, None]
    return X, F


def rel_err(H, z, H0, z0):
    num = np.sum((H - H0) ** 2) + np.sum(np.abs(z - z0) ** 2)
    den = np.sum(H0 ** 2) + np.sum(np.abs(z0) ** 2)
    return float(np.sqrt(num / den))


def sphere_rotation(k, scheme, transport, turns_steps):
    X, F = sphere(k)
    m = ops.build_mesh(F, len(X))
    adv = ca.SurfaceAdvector(m, X)
    x0 = np.array([1.0, 0.3, 0.2])
    x0 /= np.linalg.norm(x0)
    w = np.exp(-np.sum((X - x0) ** 2, 1) / 0.15)
    # z pattern: anisotropy along the projected z axis, times the bump.
    e = np.array([0, 0, 1.0]) - X * X[:, 2:3]
    e /= np.maximum(np.linalg.norm(e, axis=1), 1e-12)[:, None]
    A0 = w[:, None, None] * (np.einsum("ni,nj->nij", e, e) - np.einsum("ni,nj->nij", np.cross(adv.N, e), np.cross(adv.N, e)))
    H0 = 0.5 * w
    z0 = np.array([adv.read(v, A0[v])[1] for v in range(m.n)])
    V = np.cross(np.array([0, 0, 1.0]), X)  # rigid rotation, angular speed 1
    steps = turns_steps
    dt = 2 * np.pi / steps
    H, z = H0.copy(), z0.copy()
    for _ in range(steps):
        H, z = adv.step(H, z, V, dt, scheme, transport)
    return rel_err(H, z, H0, z0), float(H.max() / H0.max()), m.n


def torus_loxodrome(scheme, nu=48, nv=24, steps=24):
    R, r = 2.0, 0.7
    X, F, _ = ops.torus(R, r, nu, nv)
    m = ops.build_mesh(F, len(X))
    adv = ca.SurfaceAdvector(m, X)
    u = np.arctan2(X[:, 1], X[:, 0])
    rho = np.linalg.norm(X[:, :2], axis=1)
    v = np.arctan2(X[:, 2], rho - R)
    e_ring = np.column_stack([-np.sin(u), np.cos(u), np.zeros_like(u)])
    e_tube = np.column_stack([-np.sin(v) * np.cos(u), -np.sin(v) * np.sin(u), np.cos(v)])
    alpha = 0.6
    speed = 1.0
    V = speed * (np.cos(alpha) * e_tube + np.sin(alpha) * e_ring)
    u0, v0 = 0.0, 0.5
    du = np.angle(np.exp(1j * (u - u0)))
    dv = np.angle(np.exp(1j * (v - v0)))
    H0 = np.exp(-((du * R) ** 2 + (dv * r) ** 2) / 0.12)
    z0 = np.zeros(m.n, complex)
    T = 1.5
    dt = T / steps
    H, z = H0.copy(), z0.copy()
    for _ in range(steps):
        H, z = adv.step(H, z, V, dt, scheme, "pullback")
    # Exact path of the peak: du/dt = speed sin(alpha) / rho(v), dv/dt = speed cos(alpha) / r.
    uu, vv = u0, v0
    n_sub = 2000
    for _ in range(n_sub):
        rr = R + r * np.cos(vv)
        uu += speed * np.sin(alpha) / rr * T / n_sub
        vv += speed * np.cos(alpha) / r * T / n_sub
    exact = np.array([(R + r * np.cos(vv)) * np.cos(uu), (R + r * np.cos(vv)) * np.sin(uu), r * np.sin(vv)])
    wgt = np.maximum(H - 0.3 * H.max(), 0)
    peak = (wgt[:, None] * X).sum(0) / wgt.sum()
    h = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1).mean()
    return float(np.linalg.norm(peak - exact) / h), float(H.max() / H0.max())


def shear_flow(transport, gamma=0.5, t=0.4, steps=4):
    P, F = ops.strip_param(17, 17, 2.0, 2.0)
    X = ops.flat_from_param(P)  # plane y = 0
    m = ops.build_mesh(F, len(X))
    adv = ca.SurfaceAdvector(m, X)
    V = np.column_stack([gamma * X[:, 2], np.zeros(m.n), np.zeros(m.n)])  # v_x = gamma z
    ex, ez = np.array([1.0, 0, 0]), np.array([0, 0, 1.0])
    A0 = np.outer(ex, ex) - np.outer(ez, ez)
    Amb = np.repeat(A0[None], m.n, 0)
    H0 = np.array([adv.read(v, Amb[v])[0] for v in range(m.n)])
    z0 = np.array([adv.read(v, Amb[v])[1] for v in range(m.n)])
    H, z = H0.copy(), z0.copy()
    for _ in range(steps):
        H, z = adv.step(H, z, V, t / steps, "first", transport)
    # Exact: inverse flow Psi(x, z) = (x - gamma t z, z), J = [[1, -gamma t], [0, 1]] in (x, z).
    J = np.eye(3)
    J[0, 2] = -gamma * t
    Ae = J.T @ A0 @ J
    inner = np.where((np.abs(X[:, 0]) < 0.5) & (np.abs(X[:, 2]) < 0.5))[0]
    He = np.array([adv.read(v, Ae)[0] for v in inner])
    ze = np.array([adv.read(v, Ae)[1] for v in inner])
    return rel_err(H[inner], z[inner], He, ze), rel_err(H0[inner], z0[inner], He, ze)


if __name__ == "__main__":
    fails = 0

    def expect(ok, what):
        global fails
        fails += 0 if ok else 1
        print(f"  {'ok  ' if ok else 'FAIL'} {what}")

    t0 = time.time()
    print("1. Rigid rotation of a sphere, one full turn (relative error, peak kept)")
    res = {}
    for k, steps in ((6, 24), (12, 48)):
        for scheme in ("first", "bfecc", "maccormack", "maccormack_free"):
            e, peak, n = sphere_rotation(k, scheme, "pullback", steps)
            res[(k, scheme)] = e
            print(f"  cube {k:2d} ({n} verts, {steps} steps) {scheme:15s} pullback: error {e:.3f}  peak {peak:.3f}", flush=True)
        e, peak, n = sphere_rotation(k, "bfecc", "parallel", steps)
        print(f"  cube {k:2d} ({n} verts, {steps} steps) bfecc           parallel: error {e:.3f}  peak {peak:.3f}", flush=True)
    for scheme in ("first", "bfecc", "maccormack", "maccormack_free"):
        order = np.log2(res[(6, scheme)] / res[(12, scheme)])
        print(f"  observed order {scheme:10s}: {order:.2f}")
    expect(res[(12, "maccormack")] < res[(12, "first")], "MacCormack beats first order at the finer resolution")
    expect(res[(12, "bfecc")] < res[(12, "first")], "BFECC beats first order at the finer resolution")

    print("2. Torus, bump along a loxodrome (peak offset in edge lengths, height kept)")
    for scheme in ("first", "bfecc", "maccormack"):
        off, keep = torus_loxodrome(scheme)
        print(f"  {scheme:10s}: offset {off:.2f} edges, height kept {keep:.3f}", flush=True)

    print("3. Simple shear on a flat sheet (error against the exact pullback J^T A J)")
    ep, base = shear_flow("pullback")
    epar, _ = shear_flow("parallel")
    print(f"  pullback {ep:.3f}   parallel {epar:.3f}   (undeformed against exact {base:.3f})")
    expect(ep < 0.2 * base, "pullback deforms the tensor as the flow does")

    print("4. MacCormack limiter: no new extrema")
    X, F = sphere(6)
    m = ops.build_mesh(F, len(X))
    adv = ca.SurfaceAdvector(m, X)
    H0 = np.where(X[:, 0] > 0.3, 1.0, 0.0)  # a sharp step
    z0 = np.zeros(m.n, complex)
    V = np.cross(np.array([0, 0, 1.0]), X)
    H, z = H0.copy(), z0.copy()
    for _ in range(6):
        H, z = adv.step(H, z, V, 0.1, "maccormack", "parallel")
    Hb, _ = H0.copy(), None
    for _ in range(6):
        Hb, _ = adv.step(Hb, z0, V, 0.1, "bfecc", "parallel")
    print(f"  MacCormack range [{H.min():.3f}, {H.max():.3f}]   BFECC range [{Hb.min():.3f}, {Hb.max():.3f}]")
    expect(H.min() >= -1e-9 and H.max() <= 1 + 1e-9, "MacCormack stays within [0, 1]")
    print(f"{'FAILED' if fails else 'all passed'}  ({time.time() - t0:.0f}s)")
