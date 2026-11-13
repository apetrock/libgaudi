"""One Lie sync as the warm start of the contact net's Gauss-Newton
(MATH.md 19.18): iterates to each merit level, factorizations, time and the
final shape, from the start shape against from the sync's positions and
normals. The sync uses the start's metric (frame-start lengths), as a frame
would."""

import sys
import time

import numpy as np

import contact_net as cn
import cr_ops as ops
import dupin_poisson as dp
import lie_edge_net as lsn
import lie_edges as le
import lie_sync as ls
from lie_compare import cases
from live_compare import shape_error


def firsts(hist, levels=(1e-3, 1e-6, 1e-9)):
    h = np.array([q[1] for q in hist])
    out = []
    for t in levels:
        k = np.nonzero(h < t * h[0])[0]
        out.append(str(k[0]) if len(k) else "-")
    return " / ".join(out)


def run(m, X0, Xt, dist, size, tg_net, H, iterations=60, label=""):
    ring = dp.rings(m)
    pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
    g = cn.ContactNetGN(m, X0, tg_net, cn.NetConfig(), pairs)
    cold = g.run(X0.copy(), iterations=iterations)
    t0 = time.time()
    sync = ls.LieSync(m).solve(X0, H)
    ts = time.time() - t0
    g2 = cn.ContactNetGN(m, X0, tg_net, cn.NetConfig(), pairs)
    warm = g2.run(sync.X.copy(), sync.N.copy(), iterations=iterations)
    # Both merits are relative to the cold start's, so the levels compare.
    m0 = cold.history[0][1]
    hw = [(0.0, m0)] + warm.history
    f = lambda X: shape_error(X, Xt, dist, size) if Xt is not None else float("nan")
    print(f"{label:22s} sync {f(sync.X):.1e} ({ts:.2f}s)  merit at sync start {warm.history[0][1] / m0:.1e} of cold")
    print(f"    cold: to 1e-3/1e-6/1e-9 {firsts(cold.history)}  its {len(cold.history) - 1:3d}  facts {cold.factorizations:3d}  "
          f"{cold.seconds:5.2f}s  shape {f(cold.X):.1e}")
    print(f"    warm: to 1e-3/1e-6/1e-9 {firsts(hw)}  its {len(warm.history) - 1:3d}  facts {warm.factorizations:3d}  "
          f"{warm.seconds + ts:5.2f}s  shape {f(warm.X):.1e}", flush=True)


if __name__ == "__main__":
    only = [a for a in sys.argv[1:] if not a.startswith("-")]
    for name, X0, F, Xt, dist in cases():
        if only and not any(o in name for o in only):
            continue
        m = ops.build_mesh(F, len(X0))
        size = np.ptp(Xt, axis=0).max()
        ring = dp.rings(m)
        pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        run(m, X0, Xt, dist, size, cn.read_targets(m, pairs, Xt), lsn.read_spheres(m, le.hinges(m), Xt).H,
            label=name)
    # A true edit: the sphere's curvatures times 1.2 (answer: radius / 1.2).
    name, X0, F, Xt, dist = [c for c in cases() if c[0].startswith("cube->sphere")][0]
    m = ops.build_mesh(F, len(Xt))
    ring = dp.rings(m)
    pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
    lam = 1.2
    tg = cn.read_targets(m, pairs, Xt)
    tg = cn.NetTargets(tg.k * lam, tg.tau * lam, tg.c * lam)
    H = lsn.read_spheres(m, le.hinges(m), Xt).H * lam
    c = Xt.mean(0)
    Xs = c + (Xt - c) / lam
    run(m, Xt, Xs, None, np.ptp(Xs, axis=0).max(), tg, H, label="edit: sphere x 1.2")
