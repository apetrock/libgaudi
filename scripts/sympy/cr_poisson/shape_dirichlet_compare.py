"""R2: curvature-tensor Dirichlet energies (edge CR element; face dual graph,
circumcentric and barycentric) against the continuous energy
`int |grad z|^2 + 1/2 |grad tau|^2`.

Planar field (exact reference per mesh by triangle quadrature):
    z = a(y) e^{2 i pi x}, a = 1 + 0.5 sin(2 pi y), tau = 1 + 0.3 x
    integrand = pi^2 cos^2(2 pi y) + 4 pi^2 a^2 + 0.045.
Families: right-triangle grid (cocircular diagonals), equilateral lattice,
jittered grid, 4:1 high aspect, Delaunay of random points; five levels.
Curved checks: cylinder shape operator (parallel intrinsically, exact
energy 0, tests the unfolding transport), umbilic sphere (trace only, exact
0), torus shape operator (no closed form: convergence and spread over
triangulations).
"""

from __future__ import annotations

import time

import numpy as np
from scipy.spatial import Delaunay

import cr_ops as ops
import shape_dirichlet as sd

METHODS = (
    ("edge", lambda m, X: sd.edge_dirichlet(m, X, 2), sd.sample_edges),
    ("face circum", lambda m, X: sd.face_dirichlet(m, X, "circum"), sd.sample_faces),
    ("face bary", lambda m, X: sd.face_dirichlet(m, X, "bary"), sd.sample_faces),
)


def orient_up(P, F):
    F = F.copy()
    for k, (a, b, c) in enumerate(F):
        if np.cross(P[b] - P[a], P[c] - P[a])[2] < 0:
            F[k] = (a, c, b)
    return F


def grid(nx, ny, w=1.0, h=1.0, jitter=0.0, seed=0):
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.linspace(0, w, nx), np.linspace(0, h, ny), indexing="ij")
    P = np.stack([xs.ravel(), ys.ravel(), np.zeros(xs.size)], 1)
    if jitter > 0:
        inner = (P[:, 0] > 0) & (P[:, 0] < w) & (P[:, 1] > 0) & (P[:, 1] < h)
        hx, hy = w / (nx - 1), h / (ny - 1)
        P[inner, 0] += jitter * hx * rng.uniform(-1, 1, inner.sum())
        P[inner, 1] += jitter * hy * rng.uniform(-1, 1, inner.sum())
    F = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            a, b, c, d = i * ny + j, (i + 1) * ny + j, (i + 1) * ny + j + 1, i * ny + j + 1
            F += [(a, b, c), (a, c, d)]
    return P, orient_up(P, np.array(F))


def lattice(n):
    """Equilateral triangular lattice clipped to about the unit square."""
    h = 1.0 / (n - 1)
    rows = int(round(1.0 / (h * np.sqrt(3) / 2))) + 1
    pts = []
    for r in range(rows):
        off = 0.5 * h if r % 2 else 0.0
        for k in range(n + 1):
            x = off + k * h
            if x <= 1.0 + 1e-9:
                pts.append((x, r * h * np.sqrt(3) / 2))
    P2 = np.array(pts)
    tri = Delaunay(P2)
    P = np.c_[P2, np.zeros(len(P2))]
    return P, orient_up(P, tri.simplices)


def random_delaunay(n, seed=0):
    rng = np.random.default_rng(seed)
    k = n * n
    border = np.linspace(0, 1, n)
    B = np.r_[np.c_[border, np.zeros(n)], np.c_[border, np.ones(n)],
              np.c_[np.zeros(n), border], np.c_[np.ones(n), border]]
    P2 = np.unique(np.r_[B, rng.uniform(0.02, 0.98, (k, 2))], axis=0)
    tri = Delaunay(P2)
    P = np.c_[P2, np.zeros(len(P2))]
    return P, orient_up(P, tri.simplices)


# Dunavant degree-6 rule (12 points) on the reference triangle.
_DUN = np.array([
    [0.873821971016996, 0.063089014491502, 0.050844906370207],
    [0.063089014491502, 0.873821971016996, 0.050844906370207],
    [0.063089014491502, 0.063089014491502, 0.050844906370207],
    [0.501426509658179, 0.249286745170910, 0.116786275726379],
    [0.249286745170910, 0.501426509658179, 0.116786275726379],
    [0.249286745170910, 0.249286745170910, 0.116786275726379],
    [0.636502499121399, 0.310352451033785, 0.082851075618374],
    [0.310352451033785, 0.636502499121399, 0.082851075618374],
    [0.636502499121399, 0.053145049844816, 0.082851075618374],
    [0.053145049844816, 0.636502499121399, 0.082851075618374],
    [0.310352451033785, 0.053145049844816, 0.082851075618374],
    [0.053145049844816, 0.310352451033785, 0.082851075618374],
])


def planar_field(p, n):
    a = 1 + 0.5 * np.sin(2 * np.pi * p[1])
    z = a * np.array([np.cos(2 * np.pi * p[0]), np.sin(2 * np.pi * p[0])])
    tau = 1 + 0.3 * p[0]
    return np.array([[tau / 2 + z[0], z[1], 0], [z[1], tau / 2 - z[0], 0], [0, 0, 0]])


def planar_reference(P, F):
    def g(x, y):
        a = 1 + 0.5 * np.sin(2 * np.pi * y)
        return np.pi ** 2 * np.cos(2 * np.pi * y) ** 2 + 4 * np.pi ** 2 * a * a + 0.045
    E = 0.0
    for t in F:
        A, B, C = P[t[0], :2], P[t[1], :2], P[t[2], :2]
        area = 0.5 * abs(np.cross(B - A, C - A))
        for l1, l2, wq in _DUN:
            q = l1 * A + l2 * B + (1 - l1 - l2) * C
            E += area * wq * g(*q)
    return E


def energy(method, m, X, field):
    _, Lfun, samp = method
    c = samp(m, X, field).ravel()
    return float(c @ (Lfun(m, X) @ c))


def planar_study(levels=(6, 11, 21, 41, 81)):
    fams = {
        "right grid": lambda n: grid(n, n),
        "lattice": lattice,
        "jittered": lambda n: grid(n, n, jitter=0.3),
        "aspect 4:1": lambda n: grid(max(3, (n - 1) // 4 + 1), n),
        "random Delaunay": random_delaunay,
    }
    results = {}
    for fam, build in fams.items():
        for n in levels:
            P, F = build(n)
            m = ops.build_mesh(F, len(P))
            Eref = planar_reference(P, F)
            row = {}
            for meth in METHODS:
                row[meth[0]] = energy(meth, m, P, planar_field)
            results[(fam, n)] = (Eref, row, len(F))
            print(f"{fam:16s} n {n:3d} faces {len(F):6d} ref {Eref:8.3f}  "
                  + "  ".join(f"{k} {(v - Eref) / Eref:+.2e}" for k, v in row.items()),
                  flush=True)
    return results


def summarize(results, levels):
    fams = sorted({k[0] for k in results})
    print("\nconvergence order (log2 of error ratio between the last two levels) and "
          "spread of relative error across families at the finest level")
    for meth in (m[0] for m in METHODS):
        orders, finest = [], []
        for fam in fams:
            errs = [abs(results[(fam, n)][1][meth] - results[(fam, n)][0]) / results[(fam, n)][0]
                    for n in levels]
            orders.append(np.log2(errs[-2] / max(errs[-1], 1e-300)))
            finest.append((results[(fam, levels[-1])][1][meth] - results[(fam, levels[-1])][0])
                          / results[(fam, levels[-1])][0])
        print(f"  {meth:12s} order " + " ".join(f"{o:5.2f}" for o in orders)
              + f"   finest rel err min {min(finest):+.2e} max {max(finest):+.2e}")


def curved_study():
    print("\ncurved surfaces")
    # Cylinder: S = (1/R) c c^T is parallel on the (flat) cylinder: exact 0.
    for ns in (17, 33, 65):
        Pp, F = ops.strip_param(ns, ns // 2, 2.0, 1.0)
        R = 1.0
        X = ops.cylinder_from_param(Pp, R)
        m = ops.build_mesh(F, len(X))
        def cyl(p, n, R=R):  # noqa: E306
            r = np.array([p[0], p[1] - R, 0.0])
            c = np.array([-r[1], r[0], 0.0]) / np.linalg.norm(r)
            return np.outer(c, c) / R
        print(f"  cylinder ns {ns:3d}: " + "  ".join(
            f"{meth[0]} {energy(meth, m, X, cyl):.2e}" for meth in METHODS))
    # Umbilic sphere (trace only).
    X, F, _ = ops.torus(2.0, 0.7, 4, 4)  # placeholder to import; sphere below
    from closed_targets import cube, radial_ellipsoid
    Xc, Fc = cube(8)
    Xs = radial_ellipsoid(Xc, (1.0, 1.0, 1.0))
    m = ops.build_mesh(Fc, len(Xs))
    def umb(p, n):  # noqa: E306
        return np.eye(3) - np.outer(n, n)
    print("  umbilic sphere: " + "  ".join(
        f"{meth[0]} {energy(meth, m, Xs, umb):.2e}" for meth in METHODS))
    # Torus shape operator: convergence and triangulation spread.
    def torus_S(p, n, R=2.0, r=0.7):  # noqa: E306
        u = np.arctan2(p[1], p[0])
        rho = np.hypot(p[0], p[1])
        v = np.arctan2(p[2], rho - R)
        eu = np.array([-np.sin(u), np.cos(u), 0.0])
        ev = np.array([-np.sin(v) * np.cos(u), -np.sin(v) * np.sin(u), np.cos(v)])
        ku = np.cos(v) / (R + r * np.cos(v))
        kv = 1.0 / r
        return ku * np.outer(eu, eu) + kv * np.outer(ev, ev)
    rng = np.random.default_rng(0)
    for nu, nv in ((24, 12), (48, 24), (96, 48)):
        X, F, _ = ops.torus(2.0, 0.7, nu, nv)
        Ff = F.copy()
        for q in range(0, len(F), 2):  # flip the diagonal of a random half of the quads
            if rng.uniform() < 0.5:
                a, b, c = F[q]
                d = F[q + 1][2]
                Ff[q], Ff[q + 1] = (a, b, d), (b, c, d)
        out = []
        for FF in (F, Ff):
            m = ops.build_mesh(FF, len(X))
            out.append([energy(meth, m, X, torus_S) for meth in METHODS])
        print(f"  torus {nu}x{nv}: " + "  ".join(
            f"{METHODS[k][0]} {out[0][k]:.4f} / flipped {out[1][k]:.4f}"
            for k in range(len(METHODS))), flush=True)


if __name__ == "__main__":
    t = time.time()
    levels = (6, 11, 21, 41, 81)
    res = planar_study(levels)
    summarize(res, levels)
    curved_study()
    print(f"\ntotal {time.time() - t:.0f}s")
