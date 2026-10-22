"""Edge rotations for the ARAP memory (CR-ARAP diamond Procrustes and
vector-Dirichlet-weighted frame Procrustes) against vertex and face
rotations, goal Gauss-Newton, frozen system. Writes `edge_compare.png`."""

import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import base_loop as bl  # noqa: E402
import closed_targets as ct  # noqa: E402
import cr_goal as cg  # noqa: E402
import cr_ops as ops  # noqa: E402
from plot_failures import draw, residual  # noqa: E402
from run_goal import isometric_cylinder  # noqa: E402
import cr_loop as cl  # noqa: E402

MEMORIES = ("vertex_tangent", "face_tangent", "edge_arap", "edge_tangent",
            "dir_tangent", "dir_arap")


def align(X, Y):
    cx, cy = X.mean(0), Y.mean(0)
    U, _, Vt = np.linalg.svd((Y - cy).T @ (X - cx))
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    return (X - cx) @ R.T + cy


def torus_dist(R, r):
    def f(X):
        rho = np.hypot(X[:, 0], X[:, 1])
        return np.hypot(rho - R, X[:, 2]) - r
    return f


def ellipsoid_dist(abc):
    abc = np.asarray(abc)
    return lambda X: (np.sqrt(np.sum((X / abc) ** 2, 1)) - 1) * abc.min()


def cases():
    P, F = ops.strip_param(17, 6, 2.0, 1.0)
    Xf = ops.flat_from_param(P)
    yield "isometric cylinder", Xf, F, isometric_cylinder(P, 1.5), None
    P, F = ops.strip_param(13, 9, 1.6, 1.2)
    Xf = ops.flat_from_param(P)
    rr = np.linalg.norm(P, axis=1)
    d = P / np.maximum(rr, 1e-12)[:, None]
    a = rr / 2.0
    Xs = np.stack([2 * np.sin(a) * d[:, 0], 2 * (1 - np.cos(a)),
                   2 * np.sin(a) * d[:, 1]], 1)
    yield "sphere cap", Xf, F, Xs, None
    X0, F, _ = ops.torus(2.0, 0.7, 24, 12)
    yield "torus tube 0.7->0.9", X0, F, ops.torus(2.0, 0.9, 24, 12)[0], torus_dist(2.0, 0.9)
    yield "torus R 2->2.5", X0, F, ops.torus(2.5, 0.7, 24, 12)[0], torus_dist(2.5, 0.7)
    X0, F = ct.cube(6)
    R = np.sqrt(24.0 / (4 * np.pi))
    yield "cube->sphere", X0, F, ct.radial_ellipsoid(X0, (R, R, R)), ellipsoid_dist((R, R, R))
    yield ("cube->ellipsoid", X0, F, ct.radial_ellipsoid(X0, (1.8, 1.2, 0.9)),
           ellipsoid_dist((1.8, 1.2, 0.9)))


if __name__ == "__main__":
    only = sys.argv[1:]
    rows = [c for c in cases() if not only or any(o in c[0] for o in only)]
    fig = plt.figure(figsize=(3.2 * (len(MEMORIES) + 1), 3.4 * len(rows)))
    for k, (name, X0, F, Xt, dist) in enumerate(rows):
        m = ops.build_mesh(F, len(X0))
        s4 = bl.readout(m, ops.cr_derivative(m, Xt)[0], Xt)
        s3g = cl.sym_goal(s4)
        size = np.ptp(Xt, axis=0).max()
        vmax = float(np.percentile(residual(m, X0, s3g), 95))
        ax = fig.add_subplot(len(rows), len(MEMORIES) + 1, k * (len(MEMORIES) + 1) + 1,
                             projection="3d")
        draw(ax, Xt, F, residual(m, Xt, s3g), vmax, f"{name}\ntarget")
        print(name, flush=True)
        for j, mem in enumerate(MEMORIES):
            t = time.time()
            try:
                r = cg.goal_gn_loop(m, X0, X0, s4, cg.GoalGNConfig(
                    iterations=200, w_arap=1e-2, rebuild="none", memory=mem))
                X = r.X
                Xa = align(X, Xt)
                shape = (np.abs(dist(Xa)).mean() / size if dist is not None
                         else bl.rigid_error(X, Xt) / size)
                c = residual(m, X, s3g)
                line = (f"its {len(r.energy):3d} E {r.energy[0]:.1e} -> {r.energy[-1]:.1e} "
                        f"{'surf' if dist is not None else 'err'}/size {shape:.4f} "
                        f"strain {ct.edge_strain(m, X, Xt).mean():.3f}")
                title = f"{mem}\nres {c.mean():.2f}  shape {shape:.4f}"
            except Exception as err:  # noqa: BLE001
                X, line, title = None, f"failed: {err}", f"{mem}\nfailed"
            print(f"  {mem:15s} {line} t {time.time() - t:.0f}s", flush=True)
            ax = fig.add_subplot(len(rows), len(MEMORIES) + 1,
                                 k * (len(MEMORIES) + 1) + j + 2, projection="3d")
            if X is None:
                ax.set_axis_off()
                ax.set_title(title, fontsize=8)
            else:
                draw(ax, X, F, c, vmax, title)
    fig.tight_layout()
    fig.savefig("edge_compare.png", dpi=90)
    print("wrote edge_compare.png")
