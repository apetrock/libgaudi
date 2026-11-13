"""Checks for `contact_net.py`: the edge read-out reconstructs `dn` exactly,
the twist vanishes on curvature-line edges, and targets are fixed points."""

import numpy as np

import contact_net as cn
import cr_ops as ops
import dupin_poisson as dp
from lie_compare import cases

FAILS = []


def expect(cond, what):
    print(("ok   " if cond else "FAIL ") + what)
    if not cond:
        FAILS.append(what)


def reconstruction():
    rng = np.random.default_rng(0)
    X, F, _ = ops.torus(2.0, 0.7, 24, 12)
    X = X + 0.03 * rng.standard_normal(X.shape)
    m = ops.build_mesh(F, len(X))
    N = dp.vertex_normals(m, X)
    k, tau = cn.edge_readout(m.E, X, N)
    dx, nb, dxn, dxt, b = cn.edge_frames(m.E, X, N)
    dn = N[m.E[:, 1]] - N[m.E[:, 0]]
    rec = -k[:, None] * dxt - tau[:, None] * b
    expect(np.abs(rec - dn).max() < 1e-13, f"dn = -k dx_t - tau (nbar x dx) exactly ({np.abs(rec - dn).max():.1e})")


def torus_twist():
    R, r, nu, nv = 2.0, 0.7, 48, 24
    X, F, _ = ops.torus(R, r, nu, nv)
    m = ops.build_mesh(F, len(X))
    N = dp.torus_normals(X, R)
    k, tau = cn.edge_readout(m.E, X, N)
    u, v = np.divmod(np.arange(m.n), nv)
    i, j = m.E[:, 0], m.E[:, 1]
    grid = (u[i] == u[j]) | (v[i] == v[j])
    expect(np.abs(tau[grid]).max() < 1e-12 and np.abs(tau[~grid]).min() > 1e-3,
           f"twist 0 on curvature-line edges ({np.abs(tau[grid]).max():.1e}), "
           f"nonzero on diagonals (median {np.median(np.abs(tau[~grid])):.2f})")


def fixed_points():
    worst = 0
    for name, X0, F, Xt, dist in cases():
        m = ops.build_mesh(F, len(X0))
        ring = dp.rings(m)
        pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        tgt = cn.read_targets(m, pairs, Xt)
        s = cn.ContactNetSolver(m, X0, tgt, cn.NetConfig(), pairs)
        r = s.run(Xt.copy(), iterations=10, keep_best=False)
        d = np.abs(r.X - Xt).max() / np.ptp(Xt, 0).max()
        worst = max(worst, d)
        print(f"  {name:22s} fixed point drift {d:.1e}")
    expect(worst < 1e-9, f"every target is a fixed point ({worst:.1e}; the factor's conditioning floor)")


if __name__ == "__main__":
    reconstruction()
    torus_twist()
    fixed_points()
    print("all passed" if not FAILS else f"{len(FAILS)} failed")
