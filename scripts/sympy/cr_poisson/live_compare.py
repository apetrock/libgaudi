"""Live read-out (D follows the metric) with the restart protocol: solve the
frozen system to convergence, make the result the new rest, repeat. Writes
`live_compare.png` (start, target, after the first solve, after the last)."""

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
from edge_compare import align, cases  # noqa: E402
from plot_failures import draw, residual  # noqa: E402


def shape_error(X, Xt, dist, size):
    if dist is None:
        return bl.rigid_error(X, Xt) / size
    return float(np.abs(dist(align(X, Xt))).mean() / size)


if __name__ == "__main__":
    restarts = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    memory = sys.argv[2] if len(sys.argv) > 2 else "edge_arap"
    rows = list(cases())
    fig = plt.figure(figsize=(13, 3.3 * len(rows)))
    for k, (name, X0, F, Xt, dist) in enumerate(rows):
        m = ops.build_mesh(F, len(X0))
        s4 = bl.readout(m, ops.cr_derivative(m, Xt)[0], Xt)
        s3g = cl.sym_goal(s4)
        size = np.ptp(Xt, axis=0).max()
        print(f"{name}: start shape {shape_error(X0, Xt, dist, size):.4f}", flush=True)
        X = X0.copy()
        snaps = []
        for j in range(restarts):
            t = time.time()
            r = cg.goal_gn_loop(m, X, X, s4, cg.GoalGNConfig(
                iterations=150, w_arap=1e-2, rebuild="none", memory=memory,
                readout="live"))
            X = r.X
            snaps.append(X)
            print(f"  restart {j} its {len(r.energy):3d} E {r.energy[-1]:.2e} "
                  f"shape {shape_error(X, Xt, dist, size):.5f} "
                  f"strain vs target {ct.edge_strain(m, X, Xt).mean():.4f} "
                  f"t {time.time() - t:.0f}s", flush=True)
        c0 = residual(m, X0, s3g)
        vmax = float(np.percentile(c0, 95))
        panels = [(X0, c0, f"{name}\nstart"), (Xt, residual(m, Xt, s3g), "target")]
        for label, Xs in (("after solve 1", snaps[0]), (f"after solve {restarts}", snaps[-1])):
            panels.append((Xs, residual(m, Xs, s3g),
                           f"{label}\nshape {shape_error(Xs, Xt, dist, size):.4f}"))
        for j, (Xp, c, title) in enumerate(panels):
            ax = fig.add_subplot(len(rows), 4, 4 * k + j + 1, projection="3d")
            draw(ax, Xp, F, c, vmax, title)
    fig.tight_layout()
    fig.savefig("live_compare.png", dpi=90)
    print("wrote live_compare.png")
