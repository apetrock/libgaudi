"""Checks for `lie_edges.py`: signs, exactness on circles, the twist
residual of the common sphere, saddles, and the hinge projection."""

import numpy as np

import closed_targets as ct
import cr_ops as ops
import dupin_poisson as dp
import lie_edges as le
import lie_ops as lo

FAILS = []


def expect(cond, what):
    print(("ok   " if cond else "FAIL ") + what)
    if not cond:
        FAILS.append(what)


def euler_is_linear():
    p, n = np.array([0.3, -0.2, 1.0]), np.array([0.0, 0.6, 0.8])
    k1, k2, phi = -1.7, 0.4, 0.6
    s = lo.curvature_sphere(p, n, k1 * np.cos(phi) ** 2 + k2 * np.sin(phi) ** 2)
    mix = np.cos(phi) ** 2 * lo.curvature_sphere(p, n, k1) + np.sin(phi) ** 2 * lo.curvature_sphere(p, n, k2)
    expect(np.abs(s - mix).max() < 1e-14, "Euler is linear on the contact line")


def sphere_readouts():
    print("sphere (outward normals, k = -1/R):")
    errs = []
    for res in (4, 8, 16):
        Xc, F = ct.cube(res)
        R = 1.3
        X = ct.radial_ellipsoid(Xc, (R, R, R))
        m = ops.build_mesh(F, len(X))
        h = le.hinges(m)
        kp, _ = le.hinge_curvature(h, X, m=m)
        w = le.hinge_weight(h, X)
        good = w > 0.2 * w.mean()
        kc = le.contact_curvature(m.E, X, X / R)  # analytic normals
        kcm = le.contact_curvature(m.E, X, dp.vertex_normals(m, X))
        e = (np.abs(kp[good] * R + 1).max(), np.average(np.abs(kp * R + 1), weights=w),
             np.abs(kc * R + 1).max(), np.abs(kcm * R + 1).max())
        errs.append(e)
        print(f"  cube {res:2d}: hinge (weight > 0.2 mean, {good.mean():.0%} of hinges) max rel err {e[0]:.1e}, "
              f"weighted mean {e[1]:.1e}  contact (analytic n) {e[2]:.1e}  contact (mesh n) {e[3]:.1e}")
    expect(errs[-1][2] < 1e-12, "contact is exact on a sphere with exact normals")
    expect(errs[-1][0] < errs[0][0] and errs[-1][0] < 0.05, "hinge (reliable hinges) converges on a sphere")


def torus_readouts():
    R, r = 2.0, 0.7
    X, F, _ = ops.torus(R, r, 48, 24)
    m = ops.build_mesh(F, len(X))
    Na = dp.torus_normals(X, R)
    nv = 24
    u_idx, v_idx = np.divmod(np.arange(m.n), nv)
    vang = 2 * np.pi * v_idx / nv
    kc = le.contact_curvature(m.E, X, Na)
    errs_u, errs_v, res_grid, res_diag = [], [], [], []
    for e, (i, j) in enumerate(m.E):
        du = (u_idx[j] - u_idx[i]) % 48
        dv = (v_idx[j] - v_idx[i]) % nv
        k_lie, res = le.common_sphere(X[i], Na[i], X[j], Na[j])
        if dv == 0:  # along the spine: a parallel circle
            errs_u.append(abs(kc[e] + np.cos(vang[i]) / (R + r * np.cos(vang[i]))))
            res_grid.append(res)
            errs_u.append(abs(k_lie - kc[e]))
        elif du == 0:  # along the tube
            errs_v.append(abs(kc[e] + 1 / r))
            res_grid.append(res)
            errs_v.append(abs(k_lie - kc[e]))
        else:
            res_diag.append(res)
    expect(max(errs_u) < 1e-12 and max(errs_v) < 1e-12,
           f"contact exact on the torus's curvature-line edges ({max(errs_u):.1e}, {max(errs_v):.1e}); "
           "the Lie common sphere agrees")
    expect(max(res_grid) < 1e-10 and min(res_diag) > 1e-6,
           f"common-sphere residual: 0 on curvature-line edges ({max(res_grid):.1e}), "
           f"nonzero on diagonals (min {min(res_diag):.1e}, median {np.median(res_diag):.1e})")
    h = le.hinges(m)
    kp, _ = le.hinge_curvature(h, X, m=m)
    # Across a tube edge (along v... no: an edge along the spine has its
    # across direction along the tube).
    across = []
    for k in range(len(h.e)):
        i, j = h.a[k], h.b[k]
        dv = (v_idx[j] - v_idx[i]) % nv
        if dv == 0:
            across.append(abs(kp[k] * r + 1))
    print(f"  hinge across spine edges (tube curvature): max rel err {max(across):.1e} "
          "(diagonal triangulation: across is not exactly the tube direction)")


def saddle_readouts():
    n = 21
    s = np.linspace(-0.5, 0.5, n)
    P, F = ops.strip_param(n, n, 1.0, 1.0)
    X = np.stack([P[:, 0], P[:, 1], 0.5 * (P[:, 0] ** 2 - P[:, 1] ** 2)], 1)
    m = ops.build_mesh(F, len(X))
    N = dp.vertex_normals(m, X)
    if N[:, 2].mean() < 0:
        N = -N
        F = F[:, ::-1]
        m = ops.build_mesh(F, len(X))
    kc = le.contact_curvature(m.E, X, N)
    dx = X[m.E[:, 1]] - X[m.E[:, 0]]
    along_x = (np.abs(dx[:, 1]) < 1e-12)
    along_y = (np.abs(dx[:, 0]) < 1e-12)
    centre = np.linalg.norm(X[m.E[:, 0], :2], axis=1) < 0.15
    kx, ky = kc[along_x & centre].mean(), kc[along_y & centre].mean()
    expect(kx > 0.8 and ky < -0.8, f"saddle: contact signs along x {kx:+.3f}, along y {ky:+.3f} (+1 / -1)")
    h = le.hinges(m)
    kp, _ = le.hinge_curvature(h, X, m=m)
    t = X[h.b] - X[h.a]
    cen = np.linalg.norm(X[h.a, :2], axis=1) < 0.15
    ax = (np.abs(t[:, 1]) < 1e-12) & cen  # edges along x: across is y
    ay = (np.abs(t[:, 0]) < 1e-12) & cen
    expect(kp[ax].mean() < -0.5 and kp[ay].mean() > 0.5,
           f"saddle: hinge across x-edges {kp[ax].mean():+.3f} (-1), across y-edges {kp[ay].mean():+.3f} (+1)")


def hinge_projection():
    X, F, _ = ops.torus(2.0, 0.7, 24, 12)
    m = ops.build_mesh(F, len(X))
    rng = np.random.default_rng(0)
    X = X + 0.02 * rng.standard_normal(X.shape)
    h = le.hinges(m)
    tgt = le.read_targets(m, h, X)
    tgt.k_perp = tgt.k_perp * 1.3
    sol = le.EdgeSphereSolver(m, X, tgt)
    _, P = sol.projections(X)
    worst = 0
    for k in range(len(h.e)):
        Q = P[k]  # a, b, c, d (centred)
        a, b, c, d = Q
        t = (b - a) / np.linalg.norm(b - a)
        nf = np.cross(b - a, c - a)
        nf /= np.linalg.norm(nf)
        ng = np.cross(a - b, d - b)
        ng /= np.linalg.norm(ng)
        th = -np.arctan2(t @ np.cross(nf, ng), nf @ ng)
        def cot(p, q, o):
            u, v = p - o, q - o
            return u @ v / np.linalg.norm(np.cross(u, v))
        ls = 0.5 * np.linalg.norm(b - a) * (cot(a, b, c) + cot(a, b, d))
        if ls > 1e-9:
            worst = max(worst, abs(2 * np.tan(0.5 * th) / ls - tgt.k_perp[k]))
    expect(worst < 1e-9, f"hinge projection reaches its target curvature ({worst:.1e})")


if __name__ == "__main__":
    euler_is_linear()
    sphere_readouts()
    torus_readouts()
    saddle_readouts()
    hinge_projection()
    print("all passed" if not FAILS else f"{len(FAILS)} failed")
