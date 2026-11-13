"""Variational Dupin circle rows and the Moebius shift mode (MATH.md 19.20).

Variational form: at the frame start each group's target shape `T` is its
points projected onto their circle (geometric fit), then mapped by the
mode; each iterate the live points are matched to rigid copies of `T`
(Kabsch, size fixed at `T`'s), so the row is a projection onto a set fixed
for the frame and the local-global energy never rises.

Moebius shift: in the circle's unit-disk coordinates `z = (x - c) / r`,
`z -> (z - a) / (1 - conj(a) z)` with `a = amount / 2` times the unit axis
projected into the circle's plane: the conformal form of "rays from the
centre reprojected onto a circle shifted by `delta = amount r`" (equal to it
to first order; exact Moebius, so cross ratios are kept and shifts compose).

Groups: four consecutive vertices along a ring or a tube line (as
`circle4`), or a curvature-line cell `(i, j), (i+1, j), (i+1, j+1), (i, j+1)`
(cyclic, for kites). Axis fields: `tube` (the torus's mean-curvature
gradient direction), `diagonal` (`(ring + tube) / sqrt 2`, kites on cells),
`ring`.
"""

import sys
import time

import numpy as np

import cr_ops as ops
import dupin_circle_study as st

R0, r0 = st.R0, st.r0


def cells(nu, nv):
    vid = lambda i, j: (i % nu) * nv + (j % nv)
    return np.array([[vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)] for i in range(nu) for j in range(nv)])


def frames_field(X, kind):
    u = np.arctan2(X[:, 1], X[:, 0])
    rho = np.linalg.norm(X[:, :2], axis=1)
    v = np.arctan2(X[:, 2], rho - R0)
    e_ring = np.column_stack([-np.sin(u), np.cos(u), np.zeros_like(u)])
    e_tube = np.column_stack([-np.sin(v) * np.cos(u), -np.sin(v) * np.sin(u), np.cos(v)])
    if kind == "tube":
        return e_tube
    if kind == "ring":
        return e_ring
    return (e_ring + e_tube) / np.sqrt(2)


def circle_frame(P):
    """Geometric circle fit: centre, plane basis `(e1, e2)`, radius."""
    circ = st.fit_circle(P, "geometric")
    if circ is None:
        return None
    c, nu, r = circ
    _, _, e1, _ = st.plane_frame(P)
    e1 = e1 - (e1 @ nu) * nu
    e1 /= np.linalg.norm(e1)
    return c, e1, np.cross(nu, e1), r


def mobius_shift(P, axis, amount):
    """The group's points on its circle, shifted conformally along `axis`."""
    f = circle_frame(P)
    if f is None:
        return P.copy()
    c, e1, e2, r = f
    z = ((P - c) @ e1 + 1j * ((P - c) @ e2)) / r
    z /= np.abs(z)
    g = axis @ e1 + 1j * (axis @ e2)
    if abs(g) < 1e-12 or amount == 0:
        zp = z
    else:
        a = 0.5 * amount * g / abs(g)
        zp = (z - a) / (1 - np.conj(a) * z)
    return c + r * (np.outer(zp.real, e1) + np.outer(zp.imag, e2))


def kabsch(P, T):
    """Rigid copy of centred `T` closest to centred `P`."""
    H = T.T @ P
    U, S, Vt = np.linalg.svd(H)
    D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))])
    R = (U @ D @ Vt).T
    return T @ R.T


def cross_ratio(P):
    """Complex cross ratio of four (near-)coplanar points, in their plane."""
    m, nu, e1, e2 = st.plane_frame(P)
    z = (P - m) @ e1 + 1j * ((P - m) @ e2)
    return (z[0] - z[2]) * (z[1] - z[3]) / ((z[0] - z[3]) * (z[1] - z[2]))


def quad_area(P):
    return 0.5 * (np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0])) + np.linalg.norm(np.cross(P[2] - P[0], P[3] - P[0])))


def run(X0, F, G, axis_kind="tube", amount=0.1, frames=8, iters=20, mu=0.05, log_energy=False, sign=None,
        size="spread"):
    """`sign`: per-group multiplier of the axis (e.g. a checkerboard on cells:
    any quad tiles the plane by half turns about its edge midpoints, so
    alternating kites are compatible where identical ones are not)."""
    n = len(X0)
    w = np.ones(len(G))
    solve = st.system(n, G, w, mu)
    X = X0.copy()
    hist = [(*st.volume_area(X, F), *st.torus_radii(X))]
    rises = 0
    for f in range(frames):
        Xf = X.copy()
        A = frames_field(Xf, axis_kind)
        T = []
        for g_, ids in enumerate(G):
            P = Xf[ids]
            t = mobius_shift(P, A[ids].mean(0) * (1 if sign is None else sign[g_]), amount)
            t = t - t.mean(0)
            # Size fixed at the frame-start group's spread: the Moebius map
            # gives shape only (on a circle it also contracts or expands the
            # arc by |f'|, a scale the row should not inject).
            if size == "area":
                s0, st_ = np.sqrt(quad_area(P)), np.sqrt(quad_area(t))
            else:
                s0 = np.sqrt(np.sum((P - P.mean(0)) ** 2))
                st_ = np.sqrt(np.sum(t ** 2))
            T.append(t * (s0 / st_) if st_ > 0 else t)
        E_prev = np.inf
        for k in range(iters):
            Tm = np.array([kabsch(X[ids] - X[ids].mean(0), T[g]) for g, ids in enumerate(G)])
            E = sum(np.sum((X[ids] - X[ids].mean(0) - Tm[g]) ** 2) for g, ids in enumerate(G)) + mu * np.sum((X - Xf) ** 2)
            if E > E_prev * (1 + 1e-12) + 1e-14:
                rises += 1
            E_prev = E
            b = st.rhs(n, G, w, Tm) + mu * Xf
            X = np.column_stack([solve(b[:, c]) for c in range(3)])
        hist.append((*st.volume_area(X, F), *st.torus_radii(X)))
    return X, np.array(hist), rises


def kite_stats(X, C):
    """Per cell (a rectangle at the start): the mean deviation of its corner
    angles from 90 degrees, and how much the angles at corners 0 and 2
    differ from those at 1 and 3 (a kite on diagonal 0-2 keeps 1 and 3
    equal; a trapezoid pairs 0 with 1)."""
    P = X[C]
    a = np.zeros((len(C), 4))
    for k in range(4):
        u = P[:, (k - 1) % 4] - P[:, k]
        v = P[:, (k + 1) % 4] - P[:, k]
        a[:, k] = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", u, v) /
                                               (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)), -1, 1)))
    dev = np.abs(a - 90).mean()
    kite13 = np.abs(a[:, 1] - a[:, 3]).mean()
    trap01 = np.abs(a[:, 0] - a[:, 1]).mean()
    return dev, kite13, trap01


def checks(X, G, C):
    print("1. Checks")
    for name, P in (("line group", X[G[5]]), ("cell", X[C[5]])):
        T0 = mobius_shift(P, np.array([0.3, 0.2, 0.9]), 0.0)
        T1 = mobius_shift(P, np.array([0.3, 0.2, 0.9]), 0.3)
        print(f"  {name}: amount 0 moves {np.abs(T0 - P).max():.1e}; amount 0.3 moves {np.abs(T1 - P).max():.3f}, "
              f"cross ratio {cross_ratio(P):.6f} -> {cross_ratio(T1):.6f}")
    # Square cell, axis along a diagonal: a right kite (sides pair up about that diagonal).
    sq = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float) / np.sqrt(2) @ np.array(
        [[np.cos(np.pi / 4), -np.sin(np.pi / 4), 0], [np.sin(np.pi / 4), np.cos(np.pi / 4), 0], [0, 0, 1]]).T
    d02 = sq[2] - sq[0]
    kt = mobius_shift(sq, d02, 0.4)
    s = np.linalg.norm(np.roll(kt, -1, axis=0) - kt, axis=1)
    ang = lambda a, b, c: np.degrees(np.arccos((a - b) @ (c - b) / np.linalg.norm(a - b) / np.linalg.norm(c - b)))
    print(f"  square, axis on diagonal 0-2, amount 0.4: sides {np.round(s, 4)}, angles at 1 / 3 "
          f"{ang(kt[0], kt[1], kt[2]):.2f} / {ang(kt[2], kt[3], kt[0]):.2f} deg (right kite)")
    et = mobius_shift(sq, sq[1] - sq[0] + sq[2] - sq[3], 0.4)
    s = np.linalg.norm(np.roll(et, -1, axis=0) - et, axis=1)
    print(f"  square, axis along an edge direction, amount 0.4: sides {np.round(s, 4)} (isosceles trapezoid)")


if __name__ == "__main__":
    which = sys.argv[1:] or ["checks", "lines", "cells"]
    t0 = time.time()
    X, F, nu, nv = st.torus_grid()
    G, fam = st.groups(nu, nv)
    C = cells(nu, nv)
    ref = (*st.volume_area(X, F), *st.torus_radii(X))
    if "checks" in which:
        checks(X, G, C)
    if "lines" in which:
        print("2. Line groups (ring + tube lines), 8 frames x 20 iterates, mu 0.05")
        for axis_kind in ("tube", "diagonal"):
            for amount in (0.0, 0.05, 0.1):
                Xe, H, rises = run(X, F, G, axis_kind, amount)
                print(st.line(f"{axis_kind:8s} amount {amount}", H, ref) + f"  energy rises {rises}", flush=True)
    if "cells" in which:
        print("3. Cells (curvature-line quads), 8 frames")
        k0 = kite_stats(X, C)
        print(f"  start: corner angles off 90 by {k0[0]:.2f} deg")
        checker = np.array([(-1) ** ((g // nv) + (g % nv)) for g in range(len(C))])
        for label, sgn in (("uniform", None), ("checker", checker)):
            for axis_kind in ("diagonal", "tube"):
                for amount in (0.1, 0.2):
                    Xe, H, rises = run(X, F, C, axis_kind, amount, sign=sgn)
                    k = kite_stats(Xe, C)
                    print(st.line(f"{label} {axis_kind:8s} amount {amount}", H, ref) +
                          f"  rises {rises}  angles off 90 {k[0]:.2f} deg, |a1 - a3| {k[1]:.2f}, |a0 - a1| {k[2]:.2f}",
                          flush=True)
    print(f"({time.time() - t0:.0f}s)")
