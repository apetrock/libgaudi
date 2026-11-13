"""Advected curvature as geometry (MATH.md 19.21): per frame, advect the
vertex `(H, z)` field on the live mesh, map the change to edges by `B` and
hand it to `ContactNetGN` as increment targets on the mesh's own read-out.

    k   <- k_read   + d k11,     tau <- tau_read + d k12,
    c   <- c_read   + |d_ij| / 2 * d kappa_n(i, d_ij)

(`kappa_n = H + Re(z conj(w)^2)`, `w` the spoke's unit direction in the
vertex frame), so the read-out's projection loss never enters the targets.
The field is state: carried to the next frame's vertex frames through its
ambient form.

Scenarios on a torus (R 2, r 0.7) with a normal bump on the outer equator:
`still` (no flow: the round trip's drift), `rigid` (rotation about the axis,
a symmetry: the bump should orbit at the flow's rate and the torus stay),
`shear` (constant speed along the rings: inner rings turn faster, the bump
should smear into a tilted streak), and `burgers` (a ring wave in `H` riding
its own amplitude: the geometry should steepen into a crease, inner equator
first).
"""

import sys
import time

import numpy as np

import contact_net as cn
import cr_ops as ops
import curvature_advection as ca
import dupin_poisson as dp
import face_loops as fl
import hz_store as hz
import shape_dirichlet as sd

R0, r0 = 2.0, 0.7


def cr_field(m, X, vs):
    """The vertex `(H, z)` fitted to the CR edge read-out (the contact
    read-out's hinge column spikes where a dual edge degenerates)."""
    st = fl.build_stencils(m, ops.cr_derivative(m, X)[0])
    return vs.fit(sd.edge_average_readings(m, st, X))


def torus_coords(X):
    rho = np.linalg.norm(X[:, :2], axis=1)
    u = np.arctan2(X[:, 1], X[:, 0])
    v = np.arctan2(X[:, 2], rho - R0)
    d = np.hypot(rho - R0, X[:, 2]) - r0  # offset from the ideal tube
    return u, v, rho, d


def make(nu=48, nv=24, a=0.12, s2=0.25, kind="bump"):
    X, F, _ = ops.torus(R0, r0, nu, nv)
    m = ops.build_mesh(F, len(X))
    N = dp.vertex_normals(m, X)
    u, v, rho, _ = torus_coords(X)
    if kind == "bump":
        du = np.angle(np.exp(1j * u))
        dv = np.angle(np.exp(1j * v))
        X = X + (a * np.exp(-((du * R0) ** 2 + (dv * r0) ** 2) / s2))[:, None] * N
    elif kind == "wave":
        X = X + (a * np.sin(u))[:, None] * N
    return m, X


def bump_stats(X):
    u, v, rho, d = torus_coords(X)
    w = np.maximum(d - 0.3 * d.max(), 0)
    up = float(np.angle((w * np.exp(1j * u)).sum()))
    vp = float(np.angle((w * np.exp(1j * v)).sum()))
    far = np.abs(np.angle(np.exp(1j * (u - up)))) > 1.2
    spread = float(np.sqrt((w * np.angle(np.exp(1j * (v - vp))) ** 2).sum() / w.sum()))
    return up, vp, float(d.max()), float(np.abs(d[far]).mean()), spread


def ring_gradient(X, m, row, nv=24):
    """Steepest ring-direction slope of the tube offset along grid ring `row`
    (vertex `i * nv + j`; row 0 the outer equator, `nv / 2` the inner)."""
    _, _, _, d = torus_coords(X)
    i, j = m.E[:, 0], m.E[:, 1]
    sel = (i % nv == row) & (j % nv == row)
    L = np.linalg.norm(X[j] - X[i], axis=1)
    return float((np.abs(d[j] - d[i]) / L)[sel].max())


def run(scenario, frames=12, dt=1.0, scheme="bfecc", transport="pullback", its=30, omega=0.25, speed=0.5,
        c_burgers=1.0, log=True):
    m, X = make(kind="wave" if scenario == "burgers" else "bump", a=0.08 if scenario == "burgers" else 0.12)
    ring = dp.rings(m)
    pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
    adv = ca.SurfaceAdvector(m, X)
    vs = hz.VertexStore(m, X)
    s = cr_field(m, X, vs)
    H, z = s[:, 0].copy(), s[:, 1] + 1j * s[:, 2]
    Hbase = None
    if scenario == "burgers":
        # The undisplaced torus's discrete H (same read-out and fit), so
        # read-out offsets carry no speed.
        Xp, _, _ = ops.torus(R0, r0, 48, 24)
        Hbase = cr_field(m, Xp, hz.VertexStore(m, Xp))[:, 0]
        A = float(np.abs(H - Hbase).max())
        c_burgers = (R0 - r0) / (6 * A)  # inner shock near frame 6
        if log:
            print(f"    wave amplitude in H {A:.4f}, c {c_burgers:.2f}: shock t* inner ~6, outer ~{6 * (R0 + r0) / (R0 - r0):.0f} frames")
    rows = []
    t_solve = 0.0
    for f in range(frames + 1):
        if f > 0:
            u, v, rho, _ = torus_coords(X)
            e_ring = np.column_stack([-np.sin(u), np.cos(u), np.zeros_like(u)])
            adv.set_reference(e_ring)
            if scenario == "still":
                V = np.zeros_like(X)
            elif scenario == "rigid":
                V = omega * np.cross(np.array([0, 0, 1.0]), X)
            elif scenario == "shear":
                V = speed * e_ring
            else:
                V = c_burgers * (H - Hbase)[:, None] * e_ring
            Hn, zn = adv.step(H, z, V, dt, scheme, transport) if scenario != "still" else (H, z)
            dh = vs.to_edges(np.column_stack([Hn - H, (zn - z).real, (zn - z).imag]))
            tg = cn.read_targets(m, pairs, X)
            p0 = pairs[:, 0]
            dvec = X[pairs[:, 1]] - X[p0]
            w = np.einsum("ij,ij->i", dvec, vs.T[p0]) + 1j * np.einsum("ij,ij->i", dvec, vs.Bv[p0])
            L = np.abs(w)
            w = w / L
            dkn = (Hn - H)[p0] + ((zn - z)[p0] * np.conj(w) ** 2).real
            tgt = cn.NetTargets(tg.k + dh[:, 0], tg.tau + dh[:, 1], tg.c + 0.5 * L * dkn)
            t = time.time()
            res = cn.ContactNetGN(m, X, tgt, cn.NetConfig(), pairs).run(X.copy(), iterations=its)
            t_solve += time.time() - t
            # Carry the field to the new frames through its ambient form.
            Amb = adv.ambient(Hn, zn)
            X = res.X
            adv = ca.SurfaceAdvector(m, X)
            vs = hz.VertexStore(m, X)
            hz_new = [adv.read(i, Amb[i]) for i in range(m.n)]
            H = np.array([a for a, _ in hz_new])
            z = np.array([b for _, b in hz_new])
        if scenario == "burgers":
            row = (f, ring_gradient(X, m, 12), ring_gradient(X, m, 0), float((H - Hbase).max()))
            if log:
                print(f"    frame {f:2d}: inner slope {row[1]:.3f}  outer slope {row[2]:.3f}  max dH {row[3]:.3f}", flush=True)
        else:
            row = (f, *bump_stats(X))
            if log:
                print(f"    frame {f:2d}: bump at u {row[1]:+.3f} v {row[2]:+.3f}  height {row[3]:.4f}  "
                      f"torus off {row[4]:.4f}  v spread {row[5]:.3f}", flush=True)
        rows.append(row)
    return rows, t_solve / max(frames, 1)


if __name__ == "__main__":
    only = sys.argv[1:] or ["still", "rigid", "shear", "burgers"]
    t0 = time.time()
    if "still" in only:
        print("still (no flow, 6 frames): the increment round trip's drift")
        run("still", frames=6)
    if "rigid" in only:
        print("rigid rotation about the axis, omega 0.25 per frame (exact: bump at u = 0.25 f, height kept)")
        rows, ts = run("rigid")
        print(f"    solve {ts:.2f} s per frame")
    if "shear" in only:
        print("constant speed 0.5 along the rings, pullback (not an isometry: shears the torus's own curvature)")
        rows, ts = run("shear")
    if "shear_parallel" in only:
        print("constant speed 0.5 along the rings, parallel (ambient components carried and projected)")
        rows, ts = run("shear", transport="parallel")
    if "rigid_frame" in only:
        print("rigid rotation, frame transport against the ring direction (exact: u = 0.25 f)")
        rows, ts = run("rigid", transport="frame")
    if "shear_frame" in only:
        print("constant speed 0.5 along the rings, frame transport (the bump's top turns 0.185 / frame at v = 0)")
        rows, ts = run("shear", transport="frame")
    if "burgers" in only:
        print("Burgers: ring wave in the tube offset, speed (H - H_torus) along the rings, MacCormack, frame")
        run("burgers", frames=12, scheme="maccormack", transport="frame")
    print(f"({time.time() - t0:.0f}s)")
