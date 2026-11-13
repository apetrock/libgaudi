"""Per-edge rotations as the curvature store (E2; MATH.md 19.21).

Each edge carries a unit quaternion `q_e`, the rotation taking vertex `i`'s
frame to vertex `j`'s. Its rotation vector `w = 2 log q_e`, split in the
edge's Darboux frame (`t` along the edge in the mean tangent plane, `nbar`
the mean normal, `b = nbar x t`), is

    w = tau l t - k l b + kappa_g l nbar

so it carries the bending `k = k11`, the twist `tau = k12` and the edge
path's geodesic curvature `kappa_g` (intrinsic; it replaces `k22`, which
the edges in other directions fix). `kappa_g` depends on how the vertex
frames are turned about their normals: the frames are gauge-fixed by a
tangent field (the stand-in for the Dirichlet edge tangents).

Contact net with the store as unknowns: positions, vertex normals, `q_e`.
Rows: rotation `n_j - q_e n_i conj(q_e)`; goals on `(k, tau)` read from
`q_e`; a weak gauge row holding `kappa_g` at its frame-start value;
tangency; unit normals and quaternions; the length memory. Compared with
the edge store (`hes.JointContact`) and vertex frames (E3,
`quat_store.ContactFrames`) on the suite, then frustrated twist targets
with the memory held. Holonomy: around a face the rotations close up to a
turn about the normal equal to the face's angle defect (`K A`, discrete
Gauss-Bonnet); its residual is the frustration diagnostic.
"""

import sys
import time

import numpy as np

import contact_net as cn
import cr_ops as ops
import dupin_poisson as dp
import hz_edge_solvers as hes
import lie_edge_net as lsn
import lie_edges as le
import quat_store as qs
from live_compare import shape_error


# ------------------------------------------------------------- quaternions


def qmul(a, b):
    w1, x1, y1, z1 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    w2, x2, y2, z2 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2], -1)


def qconj(q):
    return q * np.array([1.0, -1.0, -1.0, -1.0])


def qrot(q, v):
    """Rotate rows `v` by the (normalized) quaternions `q`."""
    q = q / np.linalg.norm(q, axis=-1, keepdims=True)
    p = np.concatenate([np.zeros(v.shape[:-1] + (1,)), v], -1)
    return qmul(qmul(q, p), qconj(q))[..., 1:]


def qlog(q):
    """Rotation vector `2 log q` (shortest: `w >= 0`)."""
    q = q / np.linalg.norm(q, axis=-1, keepdims=True)
    q = np.where(q[..., :1] < 0, -q, q)
    v = q[..., 1:]
    s = np.linalg.norm(v, axis=-1, keepdims=True)
    ang = 2 * np.arctan2(s, q[..., :1])
    return np.where(s > 1e-12, ang * v / np.maximum(s, 1e-300), 2 * v)


def quat_from_matrix(M):
    """Unit quaternions of rotation matrices (n, 3, 3)."""
    q = np.zeros((len(M), 4))
    for k, R in enumerate(M):
        tr = np.trace(R)
        if tr > 0:
            s = 2 * np.sqrt(tr + 1)
            q[k] = [0.25 * s, (R[2, 1] - R[1, 2]) / s, (R[0, 2] - R[2, 0]) / s, (R[1, 0] - R[0, 1]) / s]
        else:
            i = int(np.argmax(np.diag(R)))
            j, l = (i + 1) % 3, (i + 2) % 3
            s = 2 * np.sqrt(1 + R[i, i] - R[j, j] - R[l, l])
            v = np.zeros(3)
            v[i] = 0.25 * s
            v[j] = (R[j, i] + R[i, j]) / s
            v[l] = (R[l, i] + R[i, l]) / s
            q[k] = [(R[l, j] - R[j, l]) / s, *v]
    return q / np.linalg.norm(q, axis=1)[:, None]


def _nrm(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def shortest_arc(a, b):
    """Unit quaternions taking unit rows `a` to `b` by the shortest arc."""
    q = np.concatenate([1 + np.einsum("ij,ij->i", a, b)[:, None], np.cross(a, b)], 1)
    return q / np.linalg.norm(q, axis=1)[:, None]


def darboux_read(xi, xj, ni, nj, q):
    """`(k, tau, kappa_g)` per edge from the edge rotation. Exact split: the
    rotation's action on `n_i` gives the bending and twist by the Rodrigues
    chord form (`dn = -k dx_t - tau nbar x dx`, independent of the gauge);
    the rest of `q`, a turn about the image normal relative to the shortest
    arc, is `kappa_g l`. (`n_j` itself enters only through the Darboux frame
    when the rotation row is not yet satisfied: here the frame uses the
    rotated normal.)"""
    ni = _nrm(ni)
    nq = qrot(q, ni)
    dx = xj - xi
    nb = _nrm(ni + nq)
    dxt = dx - np.einsum("ij,ij->i", dx, nb)[:, None] * nb
    L = np.linalg.norm(dxt, axis=1)
    t = dxt / L[:, None]
    b = np.cross(nb, t)
    dn = nq - ni
    rel = qmul(_nrm(q), qconj(shortest_arc(ni, nq)))
    g = np.einsum("ij,ij->i", qlog(rel), nq)
    return -np.einsum("ij,ij->i", dn, t) / L, -np.einsum("ij,ij->i", dn, b) / L, g / L


def gauge_frames(m, X, N):
    """Vertex frames gauge-fixed by a smooth tangent field (the projected
    reference directions of `lie_edge_net`), as rotation matrices."""
    T, B, Nn = lsn.vertex_frames(N, lsn.reference_dirs(N))
    return np.stack([T, B, Nn], -1)  # columns


def seed_edge_rotations(m, X, N):
    R = gauge_frames(m, X, N)
    i, j = m.E[:, 0], m.E[:, 1]
    rel = np.einsum("eab,ecb->eac", R[j], R[i])  # R_j R_i^T: world rotation taking frame i to frame j
    return quat_from_matrix(rel)


# ------------------------------------------------------------- contact net, edge rotations


class ContactDarboux(qs.BlockGN):
    """Unknowns `[x (3n), n (3n), q (4 nE)]`."""

    def __init__(self, m, X0, k, tau, c, w_target=10.0, w_gauge=0.1, w_memory=0.01, fade=True):
        self.m, self.k, self.tau, self.c = m, k, tau, c
        self.w_t, self.w_g, self.w_m, self.fade = w_target, w_gauge, w_memory, fade
        ring = dp.rings(m)
        self.pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        self.hl = float(np.linalg.norm(X0[m.E[:, 1]] - X0[m.E[:, 0]], axis=1).mean())
        self.memory_scale = 1.0
        self.L = None
        self.g0 = None
        self.nE = len(m.E)
        self.nvar = 6 * m.n + 4 * self.nE

    def unpack(self, z):
        n = self.m.n
        return z[:3 * n].reshape(n, 3), z[3 * n:6 * n].reshape(n, 3), z[6 * n:].reshape(-1, 4)

    def row_blocks(self, z):
        X, N, Q = self.unpack(z)
        m, P = self.m, self.pairs
        i, j = m.E[:, 0], m.E[:, 1]
        e = np.arange(self.nE)
        xo = lambda v: 3 * v
        no = lambda v: 3 * m.n + 3 * v
        qo = lambda f: 6 * m.n + 4 * f
        k, tau, c, g0, hl = self.k, self.tau, self.c, self.g0, self.hl

        def rot(V):
            ni, nj, q = V[:, 0:3], V[:, 3:6], V[:, 6:10]
            return nj - qrot(q, ni)

        def goal(V):
            kk, tt, _ = darboux_read(V[:, 0:3], V[:, 3:6], V[:, 6:9], V[:, 9:12], V[:, 12:16])
            return hl * np.column_stack([kk - k, tt - tau])

        def gauge(V):
            _, _, gg = darboux_read(V[:, 0:3], V[:, 3:6], V[:, 6:9], V[:, 9:12], V[:, 12:16])
            return (hl * (gg - g0))[:, None]

        def tang(V):
            d = V[:, 3:6] - V[:, 0:3]
            return (np.einsum("ij,ij->i", _nrm(V[:, 6:9]), d) / np.linalg.norm(d, axis=1) - c)[:, None]

        def mem(V):
            return (np.linalg.norm(V[:, 3:6] - V[:, 0:3], axis=1) - self.L)[:, None]

        V_e = np.concatenate([X[i], X[j], N[i], N[j], Q], 1)
        c_e = [(xo(i), 3), (xo(j), 3), (no(i), 3), (no(j), 3), (qo(e), 4)]
        return [
            (1.0, rot, np.concatenate([N[i], N[j], Q], 1), [(no(i), 3), (no(j), 3), (qo(e), 4)]),
            (np.sqrt(self.w_t), goal, V_e, c_e),
            (np.sqrt(self.w_g), gauge, V_e, c_e),
            (1.0, tang, np.concatenate([X[P[:, 0]], X[P[:, 1]], N[P[:, 0]]], 1),
             [(xo(P[:, 0]), 3), (xo(P[:, 1]), 3), (no(P[:, 0]), 3)]),
            (1.0, lambda V: (np.linalg.norm(V, axis=1) - 1.0)[:, None], N.copy(), [(no(np.arange(m.n)), 3)]),
            (1.0, lambda V: (np.linalg.norm(V, axis=1) - 1.0)[:, None], Q.copy(), [(qo(e), 4)]),
            (np.sqrt(self.w_m * self.memory_scale) / hl, mem, np.concatenate([X[i], X[j]], 1),
             [(xo(i), 3), (xo(j), 3)]),
        ]

    def run(self, X0, iterations=40):
        m, E = self.m, self.m.E
        N0 = dp.vertex_normals(m, X0)
        Q0 = seed_edge_rotations(m, X0, N0)
        self.g0 = darboux_read(X0[E[:, 0]], X0[E[:, 1]], N0[E[:, 0]], N0[E[:, 1]], Q0)[2]
        z0 = np.concatenate([X0.ravel(), N0.ravel(), Q0.ravel()])
        self.L = np.linalg.norm(X0[E[:, 1]] - X0[E[:, 0]], axis=1)
        self.cost0 = None

        def on_iter(z, cost):
            X, _, _ = self.unpack(z)
            self.L = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1)
            if self.cost0 is None:
                self.cost0 = max(cost, 1e-300)
            if self.fade:
                self.memory_scale = min(1.0, cost / self.cost0)

        z, st = hes.lm(z0, self.residual_blocks, self.jacobian_blocks, iterations, on_iter=on_iter, nx=3 * m.n)
        return self.unpack(z), st


# ------------------------------------------------------------- holonomy


def face_holonomy(m, X, N, Q):
    """Per face: the composed rotation around the face, split into its turn
    about the first vertex normal (inconsistency of the `kappa_g` field) and
    its off-normal part (inconsistency of bending and twist). Rotations that
    come from vertex frames compose to the identity in any gauge, so both
    are zero unless the edge rotations stop being realizable by frames."""
    E = {tuple(e): k for k, e in enumerate(m.E)}
    ang, off = np.zeros(len(m.F)), np.zeros(len(m.F))
    for f, (a, b, c) in enumerate(m.F):
        q = np.array([1.0, 0, 0, 0])
        for u, v in ((a, b), (b, c), (c, a)):
            k = E[(min(u, v), max(u, v))]
            qe = Q[k] if u < v else qconj(Q[k])
            q = qmul(qe, q)  # frame a -> ... -> back to a
        w = qlog(q[None])[0]
        n = N[a] / np.linalg.norm(N[a])
        ang[f] = w @ n
        off[f] = np.linalg.norm(w - (w @ n) * n)
    return ang, off


def angle_defect_faces(m, X):
    """Per face the integrated Gaussian curvature (a third of each vertex's
    angle defect)."""
    defect = np.full(m.n, 2 * np.pi)
    for (a, b, c) in m.F:
        for p, q, r in ((a, b, c), (b, c, a), (c, a, b)):
            u, v = X[q] - X[p], X[r] - X[p]
            defect[p] -= np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))
    return defect[m.F].sum(1) / 3


# ------------------------------------------------------------- driver


def report(label, X, Xt, dist, size, st, t, nvar):
    return (f"    {label:7s}: unknowns {nvar:5d}  shape {shape_error(X, Xt, dist, size):.2e}  its {st['its']:2d}  "
            f"facts {st['factorizations']:2d}  fill {st['fill']:7d}  factor {st['t_factor']:.2f}s  total {t:.1f}s")


def compare(name, X0, F, Xt, dist, rng):
    m = ops.build_mesh(F, len(X0))
    size = np.ptp(Xt, axis=0).max()
    ring = dp.rings(m)
    pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
    tg = cn.read_targets(m, pairs, Xt)
    h = lsn.read_spheres(m, le.hinges(m), Xt).H
    h[:, 0], h[:, 1] = tg.k, tg.tau
    h0 = lsn.read_spheres(m, le.hinges(m), X0).H
    noise = 0.1 * np.sqrt(np.mean(h[:, :2] ** 2)) * rng.normal(size=h.shape)
    print(f"  {name}", flush=True)
    for tl, hh in (("clean", h), ("noisy", h + noise)):
        out = []
        t = time.time()
        js = hes.JointContact(m, X0, hh, tg.c, hes.C_K, hes.W_COMP)
        (X, _, _), st = js.run(X0.copy(), dp.vertex_normals(m, X0), h0)
        out.append(report("store", X, Xt, dist, size, st, time.time() - t, js.nvar))
        t = time.time()
        cf = qs.ContactFrames(m, X0, hh[:, 0], hh[:, 1], tg.c)
        (X, _), st = cf.run(X0.copy())
        out.append(report("E3", X, Xt, dist, size, st, time.time() - t, cf.nvar))
        t = time.time()
        cd = ContactDarboux(m, X0, hh[:, 0], hh[:, 1], tg.c)
        (X, N, Q), st = cd.run(X0.copy())
        out.append(report("E2", X, Xt, dist, size, st, time.time() - t, cd.nvar))
        ang, off = face_holonomy(m, X, N, Q)
        out.append(f"    E2 holonomy: turn mean {np.abs(ang).mean():.2e}, off-normal mean {off.mean():.2e}")
        print(f"   {tl}\n" + "\n".join(out), flush=True)


def frustrated(rng):
    from lie_compare import cases
    name, X0, F, Xt, dist = [c for c in cases() if c[0].startswith("torus tube")][0]
    m = ops.build_mesh(F, len(X0))
    ring = dp.rings(m)
    pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
    tg = cn.read_targets(m, pairs, Xt)
    h = lsn.read_spheres(m, le.hinges(m), Xt).H
    h[:, 0], h[:, 1] = tg.k, tg.tau
    h0 = lsn.read_spheres(m, le.hinges(m), X0).H
    E = m.E
    adj = np.array([(m.fe[f, a], m.fe[f, b]) for f in range(len(m.F)) for a in range(3) for b in range(3) if a != b])
    rough = lambda v: float(np.mean(np.abs(v[adj[:, 0]] - v[adj[:, 1]])))
    L0 = np.linalg.norm(X0[E[:, 1]] - X0[E[:, 0]], axis=1).mean()
    print("Frustrated twist (torus tube, length memory held at 0.05)")
    for tau0 in (0.3, 0.6):
        hh = h.copy()
        hh[:, 1] += tau0
        res = []
        js = hes.JointContact(m, X0, hh, tg.c, hes.C_K, hes.W_COMP, w_memory=0.05)
        js.fade = False
        (Xs, Ns, _), sts = js.run(X0.copy(), dp.vertex_normals(m, X0), h0)
        res.append(("store", Xs, Ns, sts, None))
        cf = qs.ContactFrames(m, X0, hh[:, 0], hh[:, 1], tg.c, w_memory=0.05, fade=False)
        (Xf, Y), stf = cf.run(X0.copy())
        res.append(("E3", Xf, qs.qnormal(Y), stf, None))
        cd = ContactDarboux(m, X0, hh[:, 0], hh[:, 1], tg.c, w_memory=0.05, fade=False)
        (Xd, Nd, Qd), std = cd.run(X0.copy())
        res.append(("E2", Xd, Nd, std, Qd))
        for lab, X, N, st, Q in res:
            k, tau = cn.edge_readout(E, X, N)
            Lm = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1).mean()
            line = (f"  twist {tau0:.1f} {lab:5s}: realized tau {tau.mean():+.3f} (asked {hh[:, 1].mean():+.3f}), "
                    f"tau roughness {rough(tau):.3f}, mean edge {Lm / L0:.3f} x, its {st['its']}")
            if Q is not None:
                kq, tq, _ = darboux_read(X[E[:, 0]], X[E[:, 1]], N[E[:, 0]], N[E[:, 1]], Q)
                ang, off = face_holonomy(m, X, N, Q)
                line += (f"; store tau {tq.mean():+.3f}, holonomy turn {np.abs(ang).mean():.2e}, "
                         f"off-normal {off.mean():.2e}")
            print(line, flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    rng = np.random.default_rng(7)
    t0 = time.time()
    if "frustrated" in args:
        frustrated(rng)
    else:
        from lie_compare import cases
        for name, X0, F, Xt, dist in cases():
            if args and not any(a in name for a in args):
                continue
            compare(name, X0, F, Xt, dist, rng)
    print(f"({time.time() - t0:.0f}s)")
