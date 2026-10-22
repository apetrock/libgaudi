"""R3: end to end per curvature site. Read the clean target tensor at the
site, add 10% noise, fair it with the site's Dirichlet energy
(`(M + 2 lambda L) h = M h_noisy`), then recover the shape with the
Weingarten rows (live read-out), the intrinsic vector Dirichlet memory and
restarts (8 frames of 30 iterates).

Sites: "face" (face_dof; barycentric dual energy, the only finite face
weight), "edge_avg" (per-edge mean of the two faces' readings; CR spin-2
energy), and "edge_pair" (each face's reading against the edge target) for
reference. Writes `shape_site_e2e.png`.
"""

from __future__ import annotations

import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402

import base_loop as bl  # noqa: E402
import closed_targets as ct  # noqa: E402
import cr_goal as cg  # noqa: E402
import cr_loop as cl  # noqa: E402
import cr_ops as ops  # noqa: E402
import face_loops as fl  # noqa: E402
import shape_dirichlet as sd  # noqa: E402
from edge_compare import align, ellipsoid_dist, torus_dist  # noqa: E402
from plot_failures import draw  # noqa: E402

LAMBDAS = (0.03, 0.1, 0.3, 1.0)  # times the mean edge length squared


def cases():
    X0, F, _ = ops.torus(2.0, 0.7, 24, 12)
    yield "torus tube 0.7->0.9", X0, F, ops.torus(2.0, 0.9, 24, 12)[0], torus_dist(2.0, 0.9)
    X0, F = ct.cube(6)
    R = np.sqrt(24.0 / (4 * np.pi))
    yield "cube->sphere", X0, F, ct.radial_ellipsoid(X0, (R, R, R)), ellipsoid_dist((R, R, R))


def site_data(m, Xt, site):
    """Clean target (n, 3), mass M (3n), Dirichlet L at the target."""
    D0, Me = ops.cr_derivative(m, Xt)
    if site == "face":
        h = cl.sym_goal(bl.readout(m, D0, Xt))
        _, area = ops.face_normals(m, Xt)
        return h, np.repeat(area, 3), sd.face_dirichlet(m, Xt, "bary")
    st = fl.build_stencils(m, D0)
    h = sd.edge_average_readings(m, st, Xt)
    return h, np.repeat(Me, 3), sd.edge_dirichlet(m, Xt, 2)


def fair(h_noisy, M, L, lam):
    if lam == 0:
        return h_noisy.copy()
    A = (sp.diags(M) + 2 * lam * L).tocsc()
    return spla.spsolve(A, M * h_noisy.ravel()).reshape(h_noisy.shape)


def rel(a, b, M):
    d = (a - b).ravel()
    return float(np.sqrt(d @ (M * d) / (b.ravel() @ (M * b.ravel()))))


def recover(m, X0, target, site, frames=8):
    goal = target
    if site == "face":  # goal_gn_loop takes (m, 4): k11, k12, k21, k22
        goal = np.stack([target[:, 0], target[:, 1], target[:, 1], target[:, 2]], 1)
    X = X0.copy()
    for _ in range(frames):
        r = cg.goal_gn_loop(m, X, X, goal, cg.GoalGNConfig(
            iterations=30, inner=30, w_arap=1e-3, rebuild="none",
            memory="dir_intrinsic", readout="live", site=site))
        X = r.X
    return X


def result_reading(m, X, site):
    D0, _ = ops.cr_derivative(m, X)
    if site == "face":
        return cl.sym_goal(bl.readout(m, D0, X))
    return sd.edge_average_readings(m, fl.build_stencils(m, D0), X)


if __name__ == "__main__":
    t0 = time.time()
    rng = np.random.default_rng(7)
    rows = list(cases())
    fig = plt.figure(figsize=(13, 3.4 * len(rows)))
    for k, (name, X0, F, Xt, dist) in enumerate(rows):
        m = ops.build_mesh(F, len(X0))
        size = np.ptp(Xt, axis=0).max()
        ell = np.mean(np.linalg.norm(Xt[m.E[:, 1]] - Xt[m.E[:, 0]], axis=1))
        Le_t = sd.edge_dirichlet(m, Xt, 2)
        st_t = fl.build_stencils(m, ops.cr_derivative(m, Xt)[0])
        he_clean = sd.edge_average_readings(m, st_t, Xt)
        E_clean = float(he_clean.ravel() @ (Le_t @ he_clean.ravel()))
        print(f"{name}: clean edge smoothness {E_clean:.3f}", flush=True)
        best = {}
        for site in ("face", "edge_avg", "edge_pair"):
            h, M, L = site_data(m, Xt, site if site == "face" else "edge")
            noise = 0.1 * np.sqrt(np.mean(h ** 2)) * rng.normal(size=h.shape)
            lams = (None, 0.0) + (LAMBDAS if site != "edge_pair" else (0.3,))
            for lam in lams:
                target = h if lam is None else fair(h + noise, M, L, 0 if lam == 0 else lam * ell ** 2)
                terr = rel(target, h, M)
                t = time.time()
                X = recover(m, X0, target, site)
                shape = float(np.abs(dist(align(X, Xt))).mean() / size)
                he = sd.edge_average_readings(
                    m, fl.build_stencils(m, ops.cr_derivative(m, X)[0]), X)
                E_res = float(he.ravel() @ (sd.edge_dirichlet(m, X, 2) @ he.ravel()))
                rerr = rel(result_reading(m, X, "face" if site == "face" else "edge"), h, M)
                label = "clean" if lam is None else ("noisy" if lam == 0 else f"lambda {lam}")
                print(f"  {site:9s} {label:12s} target err {terr:.3f}  shape {shape:.4f}  "
                      f"curvature err {rerr:.3f}  smoothness / clean {E_res / E_clean:.3f}  "
                      f"t {time.time() - t:.0f}s", flush=True)
                if lam is not None and lam > 0 and (site not in best or shape < best[site][0]):
                    best[site] = (shape, lam, X)
        for j, site in enumerate(("face", "edge_avg", "edge_pair")):
            if site not in best:
                continue
            shape, lam, X = best[site]
            ax = fig.add_subplot(len(rows), 4, 4 * k + j + 2, projection="3d")
            c = np.abs(dist(align(X, Xt)))[m.F].mean(1) / size
            draw(ax, X, F, c, 0.01, f"{name}\n{site} lambda {lam}: shape {shape:.4f}")
        ax = fig.add_subplot(len(rows), 4, 4 * k + 1, projection="3d")
        draw(ax, Xt, F, np.zeros(len(F)), 0.01, f"{name}\ntarget")
    fig.tight_layout()
    fig.savefig("shape_site_e2e.png", dpi=90)
    print(f"wrote shape_site_e2e.png  total {time.time() - t0:.0f}s")
