"""True edits (no target mesh) through `lie_edge_net` and `contact_net`'s
Gauss-Newton (MATH.md 19.14): the start's own read-out, edited, as targets.

- sphere x 1.2: every curvature of the cube -> sphere target times 1.2; the
  answer is the same sphere at radius / 1.2 (shape error with scale).
- torus umbilic: every edge of the tube torus told `a = b = H`, `tau = 0`
  (its mean curvature): incompatible; reported, not scored.
"""

import numpy as np

import contact_net as cn
import cr_ops as ops
import dupin_poisson as dp
import lie_edge_net as lsn
import lie_edges as le
from lie_compare import cases
from live_compare import shape_error


def roundness(X):
    """Radius spread over mean radius: 0 on a sphere, whatever the vertices'
    tangential placement (the shape error also charges sliding)."""
    r = np.linalg.norm(X - X.mean(0), axis=1)
    return r.std() / r.mean()


def sphere_scale(lam=1.2):
    name, X0, F, Xt, dist = [c for c in cases() if c[0].startswith("cube->sphere")][0]
    m = ops.build_mesh(F, len(Xt))
    c = Xt.mean(0)
    Xs = c + (Xt - c) / lam
    size = np.ptp(Xs, axis=0).max()
    hin = le.hinges(m)
    sph = lsn.read_spheres(m, hin, Xt)
    sph.H = sph.H * np.array([lam, 1.0, lam])
    sph.H[:, 1] *= lam
    net = lsn.LieEdgeNet(m, Xt, sph)
    r = net.run(Xt.copy(), iterations=60)
    print(f"sphere x {lam}  lie edge net : shape {shape_error(Xt, Xs, None, size):.1e} -> "
          f"{shape_error(r.X, Xs, None, size):.1e}  its {len(r.history) - 1}  factorizations {r.factorizations}  "
          f"radius {np.linalg.norm(r.X - r.X.mean(0), axis=1).mean():.4f} (want "
          f"{np.linalg.norm(Xs - c, axis=1).mean():.4f})  roundness {roundness(Xt):.1e} -> {roundness(r.X):.1e}")
    ring = dp.rings(m)
    pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
    tg = cn.read_targets(m, pairs, Xt)
    tg = cn.NetTargets(tg.k * lam, tg.tau * lam, tg.c * lam)
    g = cn.ContactNetGN(m, Xt, tg, cn.NetConfig(), pairs)
    rg = g.run(Xt.copy(), iterations=60)
    print(f"sphere x {lam}  contact net GN: shape {shape_error(Xt, Xs, None, size):.1e} -> "
          f"{shape_error(rg.X, Xs, None, size):.1e}  its {len(rg.history) - 1}  factorizations {rg.factorizations}  "
          f"radius {np.linalg.norm(rg.X - rg.X.mean(0), axis=1).mean():.4f}  roundness {roundness(rg.X):.1e}")


def torus_umbilic():
    name, X0, F, Xt, dist = [c for c in cases() if c[0].startswith("torus tube")][0]
    m = ops.build_mesh(F, len(X0))
    hin = le.hinges(m)
    sph = lsn.read_spheres(m, hin, X0)
    Hm = 0.5 * (sph.H[:, 0] + sph.H[:, 2])
    sph.H = np.column_stack([Hm, np.zeros_like(Hm), Hm])
    net = lsn.LieEdgeNet(m, X0, sph)
    r = net.run(X0.copy(), iterations=60)
    after = lsn.read_spheres(m, hin, r.X).H
    tg = sph.H
    print(f"torus umbilic edit (incompatible): merit {r.history[0][1]:.2e} -> {r.best_merit:.2e}  "
          f"its {len(r.history) - 1}  |a - b| mean {np.abs(lsn.read_spheres(m, hin, X0).H[:, 0] - lsn.read_spheres(m, hin, X0).H[:, 2]).mean():.3f} -> "
          f"{np.abs(after[:, 0] - after[:, 2]).mean():.3f}  |tau| {np.abs(after[:, 1]).mean():.3f}  "
          f"finite {np.all(np.isfinite(r.X))}")


if __name__ == "__main__":
    sphere_scale()
    torus_umbilic()
