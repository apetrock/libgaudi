"""Checks for `lie_sync.py` (MATH.md 19.17)."""

import numpy as np

import cr_ops as ops
import lie_edge_net as lsn
import lie_edges as le
import lie_sync as ls
from lie_compare import cases
from live_compare import shape_error

fails = 0


def expect(ok, what, value=None):
    global fails
    if not ok:
        fails += 1
    print(f"  {'ok  ' if ok else 'FAIL'} {what}" + ("" if value is None else f"  ({value:.2e})"))


def rot(axis, ang):
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


for name, X0, F, Xt, dist in cases():
    if not any(t in name for t in ("tube", "R 2", "cube->sphere", "cylinder")):
        continue
    m = ops.build_mesh(F, len(X0))
    size = np.ptp(Xt, axis=0).max()
    tg = lsn.read_spheres(m, le.hinges(m), Xt)
    S = ls.LieSync(m)
    print(name)
    # The model's relative transforms against the target mesh's own frames.
    th = ls.polar_angles(m, ls.edge_lengths(m, Xt), S.fans, S.closed)
    G = ls.relative_transforms(m, ls.edge_lengths(m, Xt), tg.H, th)
    F0 = ls.initial_frames(m, Xt, S.fans)
    Ga = np.einsum("nij,njk->nik", ls._inv(F0[m.E[:, 0]]), F0[m.E[:, 1]])
    dR = np.median(np.abs(G[:, :3, :3] - Ga[:, :3, :3]).max(axis=(1, 2)))
    print(f"       relative rotations vs the target's own frames: median {dR:.1e}")
    r = S.solve(X0, tg.H, Xmetric=Xt)
    e = shape_error(r.X, Xt, dist, size)
    expect(e < 1e-2, "one sync from the start with the target metric reproduces the target", e)
    Rm, tv = rot([0.3, -0.5, 0.8], 0.9), np.array([0.4, -1.1, 2.0])
    r2 = S.solve(X0 @ Rm.T + tv, tg.H, Xmetric=Xt)
    d = np.abs(r2.X - (r.X @ Rm.T + tv)).max() / size
    expect(d < 1e-8, "rigid motion of the start moves the result rigidly", d)
    expect(r.factorizations <= 2, f"factorizations {r.factorizations}")

print("FAILED" if fails else "all passed")
