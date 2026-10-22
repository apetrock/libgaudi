"""Shapes reached on the non-isometric targets (goal GN, vertex-tangent
memory), frozen versus rebuilt every outer step. Saves `failures.npz` and
`failures.png` (rows: cases; columns: start, target, frozen, rebuilt; faces
colored by the curvature residual |S_f - s_g,f|)."""

import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import base_loop as bl  # noqa: E402
import closed_targets as ct  # noqa: E402
import cr_goal as cg  # noqa: E402
import cr_loop as cl  # noqa: E402
import cr_ops as ops  # noqa: E402


def cases():
    X0, F, _ = ops.torus(2.0, 0.7, 24, 12)
    yield "torus tube 0.7 -> 0.9", X0, F, ops.torus(2.0, 0.9, 24, 12)[0]
    yield "torus R 2 -> 2.5", X0, F, ops.torus(2.5, 0.7, 24, 12)[0]
    X0, F = ct.cube(6)
    R = np.sqrt(24.0 / (4 * np.pi))
    yield "cube -> sphere", X0, F, ct.radial_ellipsoid(X0, (R, R, R))
    yield "cube -> ellipsoid", X0, F, ct.radial_ellipsoid(X0, (1.8, 1.2, 0.9))


def residual(m, X, s3g):
    D0, _ = ops.cr_derivative(m, X)
    S = ops.face_shape(m, X, D0)
    T, B, _ = ops.face_frames(m, X)
    return np.linalg.norm(cl.sym_goal(ops.face_coeffs(S, T, B)) - s3g, axis=1)


def draw(ax, X, F, c, vmax, title):
    Xc = X - X.mean(0)
    surf = ax.plot_trisurf(Xc[:, 0], Xc[:, 1], Xc[:, 2], triangles=F,
                           cmap="viridis", linewidth=0.1, edgecolor="k",
                           antialiased=True)
    surf.set_array(c)
    surf.set_clim(0, vmax)
    r = np.abs(Xc).max()
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_zlim(-r, r)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.view_init(elev=30, azim=35)
    ax.set_title(title, fontsize=8)


if __name__ == "__main__":
    iterations = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    rows = list(cases())
    fig = plt.figure(figsize=(16, 4 * len(rows)))
    store = {}
    for k, (name, X0, F, Xt) in enumerate(rows):
        m = ops.build_mesh(F, len(X0))
        s4 = bl.readout(m, ops.cr_derivative(m, Xt)[0], Xt)
        s3g = cl.sym_goal(s4)
        size = np.ptp(Xt, axis=0).max()
        out = {}
        for rb in ("none", "D+G0"):
            t = time.time()
            try:
                r = cg.goal_gn_loop(m, X0, X0, s4, cg.GoalGNConfig(
                    iterations=iterations, w_arap=1e-2, rebuild=rb,
                    memory="vertex_tangent"))
                out[rb] = (r.X, r.energy[-1])
            except RuntimeError as err:
                out[rb] = (None, str(err))
            X = out[rb][0]
            msg = ("failed: " + out[rb][1]) if X is None else (
                f"E {out[rb][1]:.2e} err/size {bl.rigid_error(X, Xt) / size:.3f} "
                f"strain {ct.edge_strain(m, X, Xt).mean():.3f}")
            print(f"{name:22s} {rb:5s} {msg} t {time.time() - t:.0f}s", flush=True)
        c0 = residual(m, X0, s3g)
        vmax = float(np.percentile(c0, 95))
        panels = [(X0, c0, f"{name}\nstart  E-res mean {c0.mean():.2f}"),
                  (Xt, residual(m, Xt, s3g), "target")]
        for rb in ("none", "D+G0"):
            X, info = out[rb]
            if X is None:
                panels.append((None, None, f"{rb}: {info}"))
                continue
            c = residual(m, X, s3g)
            panels.append((X, c, f"{'frozen' if rb == 'none' else 'rebuilt (D+G0)'}"
                                 f"  res mean {c.mean():.2f}\n"
                                 f"err/size {bl.rigid_error(X, Xt) / size:.3f}"))
            store[f"{k}_{rb}"] = X
        store[f"{k}_start"], store[f"{k}_target"], store[f"{k}_F"] = X0, Xt, F
        for j, (X, c, title) in enumerate(panels):
            ax = fig.add_subplot(len(rows), 4, 4 * k + j + 1, projection="3d")
            if X is None:
                ax.set_axis_off()
                ax.set_title(title, fontsize=8)
                continue
            draw(ax, X, F, c, vmax, title)
    np.savez("failures.npz", **store)
    fig.tight_layout()
    fig.savefig("failures.png", dpi=110)
    print("wrote failures.png")
