"""Memory comparison under the working solver: live read-out, J^T W J on the
left, 30 iterates per frame, rest reseeded every frame (the restart
protocol). Memories: CR-ARAP (`edge_arap`), Stein's vector Dirichlet on the
intrinsic edge frames (`dir_intrinsic`), vertex tangents (baseline), and the
earlier extrinsic `dir_arap`. Writes `memory_restart_compare.npz` and
`memory_restart_compare.png`.

Run: python3 memory_restart_compare.py [frames]
"""

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
from edge_compare import cases  # noqa: E402
from live_compare import shape_error  # noqa: E402
from plot_failures import draw, residual  # noqa: E402

MEMORIES = ("edge_arap", "dir_intrinsic", "vertex_tangent", "dir_arap")
WEIGHTS = (1e-3, 1e-2, 1e-1)


def dihedrals(m, X):
    N, _ = ops.face_normals(m, X)
    th = np.zeros(len(m.E))
    for e, fs in enumerate(m.edge_faces):
        if len(fs) == 2:
            th[e] = np.arccos(np.clip(N[fs[0]] @ N[fs[1]], -1, 1))
    return th


def run(m, X0, Xt, dist, s4, memory, w, frames):
    size = np.ptp(Xt, axis=0).max()
    th_t = dihedrals(m, Xt)
    X = X0.copy()
    rows = []
    t0 = time.time()
    for _ in range(frames):
        r = cg.goal_gn_loop(m, X, X, s4, cg.GoalGNConfig(
            iterations=30, inner=30, w_arap=w, rebuild="none", memory=memory,
            readout="live"))
        X = r.X
        rows.append((shape_error(X, Xt, dist, size), r.energy[-1],
                     ct.edge_strain(m, X, Xt).mean(),
                     float(np.sqrt(np.mean((dihedrals(m, X) - th_t) ** 2)))))
    return X, np.array(rows), time.time() - t0


if __name__ == "__main__":
    frames = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    rows = list(cases())
    store, best = {}, {}
    for name, X0, F, Xt, dist in rows:
        m = ops.build_mesh(F, len(X0))
        s4 = bl.readout(m, ops.cr_derivative(m, Xt)[0], Xt)
        print(f"{name}: start shape {shape_error(X0, Xt, dist, np.ptp(Xt, 0).max()):.4f}",
              flush=True)
        for mem in MEMORIES:
            for w in WEIGHTS:
                try:
                    X, hist, t = run(m, X0, Xt, dist, s4, mem, w, frames)
                except Exception as err:  # noqa: BLE001
                    print(f"  {mem:14s} w {w:.0e} failed: {err}", flush=True)
                    continue
                store[f"{name}|{mem}|{w}"] = hist
                store[f"{name}|{mem}|{w}|X"] = X
                se = hist[:, 0]
                print(f"  {mem:14s} w {w:.0e} shape "
                      + " ".join(f"{v:.1e}" for v in se[[0, 1, 2, -1]])
                      + f"  E {hist[-1, 1]:.1e} strain {hist[-1, 2]:.4f} "
                        f"dihedral rms {hist[-1, 3]:.2e} t {t:.0f}s", flush=True)
                key = (name, mem)
                if key not in best or se[-1] < best[key][0]:
                    best[key] = (se[-1], w, X, hist)
    np.savez("memory_restart_compare.npz", **{k: v for k, v in store.items()})

    fig = plt.figure(figsize=(3.2 * (len(MEMORIES) + 1), 3.4 * len(rows)))
    for k, (name, X0, F, Xt, dist) in enumerate(rows):
        m = ops.build_mesh(F, len(X0))
        s3g = cl.sym_goal(bl.readout(m, ops.cr_derivative(m, Xt)[0], Xt))
        vmax = float(np.percentile(residual(m, X0, s3g), 95))
        ax = fig.add_subplot(len(rows), len(MEMORIES) + 1, k * (len(MEMORIES) + 1) + 1,
                             projection="3d")
        draw(ax, Xt, F, residual(m, Xt, s3g), vmax, f"{name}\ntarget")
        for j, mem in enumerate(MEMORIES):
            ax = fig.add_subplot(len(rows), len(MEMORIES) + 1,
                                 k * (len(MEMORIES) + 1) + j + 2, projection="3d")
            if (name, mem) not in best:
                ax.set_axis_off()
                ax.set_title(f"{mem}\nfailed", fontsize=8)
                continue
            se, w, X, hist = best[(name, mem)]
            draw(ax, X, F, residual(m, X, s3g), vmax,
                 f"{mem} w {w:.0e}\nshape {se:.1e} dih {hist[-1, 3]:.1e}")
    fig.tight_layout()
    fig.savefig("memory_restart_compare.png", dpi=90)

    fig, axes = plt.subplots(1, len(rows), figsize=(3.4 * len(rows), 3.0))
    for k, (name, *_rest) in enumerate(rows):
        for mem in MEMORIES:
            if (name, mem) in best:
                se, w, X, hist = best[(name, mem)]
                axes[k].semilogy(np.arange(1, len(hist) + 1),
                                 np.maximum(hist[:, 0], 1e-12), label=f"{mem} {w:.0e}")
        axes[k].set_title(name, fontsize=8)
        axes[k].set_xlabel("frame")
    axes[0].set_ylabel("shape error / size")
    axes[-1].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig("memory_restart_compare_curves.png", dpi=90)
    print("wrote memory_restart_compare.png, memory_restart_compare_curves.png")
