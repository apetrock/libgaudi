"""Quaternion frames as the curvature store (MATH.md 19.22).

Each vertex carries a quaternion `Y` (rotation: the vertex frame and its
normal `n = Y e3 conj(Y)`; scale `|Y|^2`, unused here). Across an edge the
frames' normals give the edge's curvature by the Rodrigues relation with
twist, `dn = -k dx_t - tau nbar x dx`:

    k11 = -dn . t / |dx_t|,   k12 = -dn . b / |dx_t|      (t lower -> higher, b = nbar x t)

so the frame field replaces the edge store `(k11, k12, k22)`: the edges in
three directions at a vertex fix the full tensor, and the rotation about
the normal is a gauge. Unlike the store, a frame field only represents
curvature some frame field realizes (no slack between store and mesh).

1. Contact net: store (`hz_edge_solvers.JointContact`: positions, normals,
   edge store, goal rows on the store) against frames (positions and `Y`;
   the goals act on the Rodrigues and tangency rows directly).
2. Weingarten: store (`hz_edge_solvers.JointWeingarten`: positions and edge
   store, the CR read-out coupled to the store) against frames (positions
   and `Y`; the CR read-out's `(k11, k12)` coupled to the frames' Rodrigues
   curvature, goal rows on the frames' curvature, tangency rows keeping the
   frames on the surface).
3. Frustrated targets (curvature-space buckling): a twist added to the
   torus target that no surface realizes; how each store spends it.
"""

import sys
import time

import numpy as np
import scipy.sparse as sp

import contact_net as cn
import cr_ops as ops
import dupin_poisson as dp
import face_loops as fl
import hz_edge_solvers as hes
import lie_edge_net as lsn
import lie_edges as le
import shape_dirichlet as sd
from live_compare import shape_error


# ------------------------------------------------------------- quaternions


def qnormal(Y):
    """`Y e3 conj(Y) / |Y|^2` for rows `(w, x, y, z)`."""
    q = Y / np.linalg.norm(Y, axis=1)[:, None]
    w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    return np.column_stack([2 * (x * z + w * y), 2 * (y * z - w * x), 1 - 2 * (x * x + y * y)])


def quat_from_normal(N):
    """The shortest-arc rotation taking `e3` to each `N`."""
    Y = np.column_stack([1 + N[:, 2], -N[:, 1], N[:, 0], np.zeros(len(N))])
    flip = Y[:, 0] < 1e-8
    Y[flip] = np.array([0.0, 1.0, 0.0, 0.0])
    return Y / np.linalg.norm(Y, axis=1)[:, None]


def _nrm(v):
    return v / np.linalg.norm(v, axis=1)[:, None]


def frame_curvature(xi, xj, ni, nj):
    """`(k11, k12)` per edge from endpoint positions and normals."""
    dx, dn = xj - xi, nj - ni
    nb = _nrm(ni + nj)
    dxt = dx - np.einsum("ij,ij->i", dx, nb)[:, None] * nb
    L = np.linalg.norm(dxt, axis=1)
    t = dxt / L[:, None]
    b = np.cross(nb, t)
    return np.column_stack([-np.einsum("ij,ij->i", dn, t) / L, -np.einsum("ij,ij->i", dn, b) / L])


# ------------------------------------------------------------- block FD


class BlockGN:
    """Row blocks (`row_blocks(z)`) are `(w, f, V, cols)`: `V` the rows'
    local unknowns side by side, `cols` a list of `(start_array, width)`
    (each group's global start per row)."""

    fd = 1e-7

    def residual_blocks(self, z):
        return np.concatenate([(w * f(V)).ravel() for w, f, V, _ in self.row_blocks(z)])

    def jacobian_blocks(self, z):
        rows, cols, vals = [], [], []
        r0 = 0
        for w, f, V, groups in self.row_blocks(z):
            base = f(V)
            nr, nc = base.shape
            off = 0
            for start, width in groups:
                for c in range(width):
                    Vp, Vm = V.copy(), V.copy()
                    Vp[:, off + c] += self.fd
                    Vm[:, off + c] -= self.fd
                    D = w * (f(Vp) - f(Vm)) / (2 * self.fd)
                    for q in range(nc):
                        rows.append(r0 + np.arange(nr) * nc + q)
                        cols.append(start + c)
                        vals.append(D[:, q])
                off += width
            r0 += nr * nc
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(r0, self.nvar))


# ------------------------------------------------------------- contact net, frames


class ContactFrames(BlockGN):
    """Unknowns `[x (3n), Y (4n)]`. Rows: Rodrigues with the goal `(k, tau)`
    per edge, tangency with the goal offsets `c` per spoke, the fading
    length memory (as `JointStore`)."""

    def __init__(self, m, X0, k, tau, c, w_memory=0.01, fade=True):
        self.m, self.k, self.tau, self.c = m, k, tau, c
        self.fade = fade
        ring = dp.rings(m)
        self.pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        self.hl = float(np.linalg.norm(X0[m.E[:, 1]] - X0[m.E[:, 0]], axis=1).mean())
        self.w_m = w_memory
        self.memory_scale = 1.0
        self.L = None
        self.nvar = 7 * m.n

    def unpack(self, z):
        n = self.m.n
        return z[:3 * n].reshape(n, 3), z[3 * n:].reshape(n, 4)

    def row_blocks(self, z):
        X, Y = self.unpack(z)
        m, P = self.m, self.pairs
        i, j = m.E[:, 0], m.E[:, 1]
        xo = lambda v: 3 * v
        yo = lambda v: 3 * m.n + 4 * v
        k, tau, c = self.k, self.tau, self.c

        def rod(V):
            xi, xj, ni, nj = V[:, 0:3], V[:, 3:6], qnormal(V[:, 6:10]), qnormal(V[:, 10:14])
            dx, dn = xj - xi, nj - ni
            nb = _nrm(ni + nj)
            dxt = dx - np.einsum("ij,ij->i", dx, nb)[:, None] * nb
            return k[:, None] * dxt + dn + tau[:, None] * np.cross(nb, dx)

        def tang(V):
            d = V[:, 3:6] - V[:, 0:3]
            return (np.einsum("ij,ij->i", qnormal(V[:, 6:10]), d) / np.linalg.norm(d, axis=1) - c)[:, None]

        def mem(V):
            return (np.linalg.norm(V[:, 3:6] - V[:, 0:3], axis=1) - self.L)[:, None]

        return [
            (1.0, rod, np.concatenate([X[i], X[j], Y[i], Y[j]], 1),
             [(xo(i), 3), (xo(j), 3), (yo(i), 4), (yo(j), 4)]),
            (1.0, tang, np.concatenate([X[P[:, 0]], X[P[:, 1]], Y[P[:, 0]]], 1),
             [(xo(P[:, 0]), 3), (xo(P[:, 1]), 3), (yo(P[:, 0]), 4)]),
            (np.sqrt(self.w_m * self.memory_scale) / self.hl, mem, np.concatenate([X[i], X[j]], 1),
             [(xo(i), 3), (xo(j), 3)]),
        ]

    def run(self, X0, iterations=40):
        m, E = self.m, self.m.E
        Y0 = quat_from_normal(dp.vertex_normals(m, X0))
        z0 = np.concatenate([X0.ravel(), Y0.ravel()])
        self.L = np.linalg.norm(X0[E[:, 1]] - X0[E[:, 0]], axis=1)
        self.cost0 = None

        def on_iter(z, cost):
            X, _ = self.unpack(z)
            self.L = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1)
            if self.cost0 is None:
                self.cost0 = max(cost, 1e-300)
            if self.fade:
                self.memory_scale = min(1.0, cost / self.cost0)

        z, st = hes.lm(z0, self.residual_blocks, self.jacobian_blocks, iterations, on_iter=on_iter, nx=3 * m.n)
        return self.unpack(z), st


# ------------------------------------------------------------- Weingarten, frames


class WeingartenFrames(BlockGN):
    """Unknowns `[x (3n), Y (4n)]`. Rows: the CR read-out's `(k11, k12)`
    coupled to the frames' Rodrigues curvature, goal rows on the frames'
    curvature, tangency (goal offsets `c`), fading length memory."""

    def __init__(self, m, X0, g, c, w_c=1.0, w_t=10.0, w_m=0.01, increments=True):
        self.m, self.g, self.c = m, g, c
        self.increments = increments
        self.delta = np.zeros((len(m.E), 2))
        self.w_c, self.w_t, self.w_m = w_c, w_t, w_m
        self.st = fl.build_stencils(m, ops.cr_derivative(m, X0)[0])
        A = __import__("hz_store").edge_dual_area(m, X0)
        self.sa = np.sqrt(A / A.mean())
        self.cnt = np.zeros(len(m.E))
        for q in range(3):
            np.add.at(self.cnt, m.fe[:, q], 1)
        ring = dp.rings(m)
        self.pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        self.hl = float(np.linalg.norm(X0[m.E[:, 1]] - X0[m.E[:, 0]], axis=1).mean())
        self.memory_scale = 1.0
        self.L = None
        self.nvar = 7 * m.n

    def unpack(self, z):
        n = self.m.n
        return z[:3 * n].reshape(n, 3), z[3 * n:].reshape(n, 4)

    def readout(self, X):
        Rf = sd.face_edge_readings(self.st, X[self.st.ids])
        acc = np.zeros((len(self.m.E), 3))
        for q in range(3):
            np.add.at(acc, self.m.fe[:, q], Rf[:, q])
        return acc / self.cnt[:, None]

    def blocks(self, z):
        X, Y = self.unpack(z)
        m, P = self.m, self.pairs
        i, j = m.E[:, 0], m.E[:, 1]
        xo = lambda v: 3 * v
        yo = lambda v: 3 * m.n + 4 * v
        sa, g, c = self.sa, self.g, self.c

        def fk(V):
            return frame_curvature(V[:, 0:3], V[:, 3:6], qnormal(V[:, 6:10]), qnormal(V[:, 10:14]))

        delta = self.delta

        # Goals in the frames' own discretization: the CR goal minus the
        # frame-start offset between the CR and Rodrigues read-outs.
        def goal(V):
            return sa[:, None] * (fk(V) - (g[:, :2] - delta))

        def tang(V):
            d = V[:, 3:6] - V[:, 0:3]
            return (np.einsum("ij,ij->i", qnormal(V[:, 6:10]), d) / np.linalg.norm(d, axis=1) - c)[:, None]

        def mem(V):
            return (np.linalg.norm(V[:, 3:6] - V[:, 0:3], axis=1) - self.L)[:, None]

        V_edge = np.concatenate([X[i], X[j], Y[i], Y[j]], 1)
        cols_edge = [(xo(i), 3), (xo(j), 3), (yo(i), 4), (yo(j), 4)]
        return [
            (np.sqrt(self.w_t), goal, V_edge, cols_edge),
            (1.0, tang, np.concatenate([X[P[:, 0]], X[P[:, 1]], Y[P[:, 0]]], 1),
             [(xo(P[:, 0]), 3), (xo(P[:, 1]), 3), (yo(P[:, 0]), 4)]),
            (np.sqrt(self.w_m * self.memory_scale) / self.hl, mem, np.concatenate([X[i], X[j]], 1),
             [(xo(i), 3), (xo(j), 3)]),
        ], V_edge, cols_edge, fk

    # The coupling row `sqrt(w_c) sa (R(x)[:, :2] - fk)`: the read-out part by
    # slot-wise differences over the face stencils (as `JointWeingarten`), the
    # frame part per edge.
    def residual(self, z):
        X, Y = self.unpack(z)
        blk, V_edge, _, fk = self.blocks(z)
        rc = np.sqrt(self.w_c) * self.sa[:, None] * (self.readout(X)[:, :2] - fk(V_edge) - self.delta)
        return np.concatenate([rc.ravel()] + [(w * f(V)).ravel() for w, f, V, _ in blk])

    def jacobian(self, z):
        X, Y = self.unpack(z)
        m, st, nE = self.m, self.st, len(self.m.E)
        blk, V_edge, cols_edge, fk = self.blocks(z)
        rows, cols, vals = [], [], []
        sc = np.sqrt(self.w_c)
        P = X[st.ids]
        for s in range(6):
            live = st.mask[:, s] > 0
            for c in range(3):
                Pp, Pm = P.copy(), P.copy()
                Pp[:, s, c] += self.fd
                Pm[:, s, c] -= self.fd
                dR = (sd.face_edge_readings(st, Pp) - sd.face_edge_readings(st, Pm)) / (2 * self.fd)
                for q in range(3):
                    e = m.fe[:, q]
                    coef = sc * self.sa[e] / self.cnt[e]
                    for r in range(2):
                        rows.append((2 * e + r)[live])
                        cols.append((3 * st.ids[:, s] + c)[live])
                        vals.append((coef * dR[:, q, r])[live])
        off = 0
        for start, width in cols_edge:
            for c in range(width):
                Vp, Vm = V_edge.copy(), V_edge.copy()
                Vp[:, off + c] += self.fd
                Vm[:, off + c] -= self.fd
                D = -sc * self.sa[:, None] * (fk(Vp) - fk(Vm)) / (2 * self.fd)
                for r in range(2):
                    rows.append(2 * np.arange(nE) + r)
                    cols.append(start + c)
                    vals.append(D[:, r])
            off += width
        Jc = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                           shape=(2 * nE, self.nvar))
        Jr = self.jacobian_blocks(z)
        return sp.vstack([Jc, Jr]).tocsr()

    def row_blocks(self, z):
        return self.blocks(z)[0]

    def run(self, X0, iterations=40):
        m, E = self.m, self.m.E
        Y0 = quat_from_normal(dp.vertex_normals(m, X0))
        z0 = np.concatenate([X0.ravel(), Y0.ravel()])
        self.L = np.linalg.norm(X0[E[:, 1]] - X0[E[:, 0]], axis=1)
        if self.increments:
            N0 = qnormal(Y0)
            self.delta = self.readout(X0)[:, :2] - frame_curvature(X0[E[:, 0]], X0[E[:, 1]], N0[E[:, 0]], N0[E[:, 1]])
        self.cost0 = None

        def on_iter(z, cost):
            X, _ = self.unpack(z)
            self.L = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1)
            if self.cost0 is None:
                self.cost0 = max(cost, 1e-300)
            self.memory_scale = min(1.0, cost / self.cost0)

        z, st = hes.lm(z0, self.residual, self.jacobian, iterations, on_iter=on_iter, nx=3 * m.n)
        return self.unpack(z), st


# ------------------------------------------------------------- driver


def report(label, X, Xt, dist, size, st, t, nvar):
    return (f"    {label:7s}: unknowns {nvar:5d}  shape {shape_error(X, Xt, dist, size):.2e}  its {st['its']:2d}  "
            f"facts {st['factorizations']:2d}  fill {st['fill']:7d}  factor {st['t_factor']:.2f}s  total {t:.1f}s")


def run_case(name, X0, F, Xt, dist, rng, which):
    m = ops.build_mesh(F, len(X0))
    size = np.ptp(Xt, axis=0).max()
    ring = dp.rings(m)
    pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
    tg = cn.read_targets(m, pairs, Xt)
    print(f"  {name}", flush=True)
    if "contact" in which:
        h = lsn.read_spheres(m, le.hinges(m), Xt).H
        h[:, 0], h[:, 1] = tg.k, tg.tau
        h0 = lsn.read_spheres(m, le.hinges(m), X0).H
        noise = 0.1 * np.sqrt(np.mean(h[:, :2] ** 2)) * rng.normal(size=h.shape)
        for tl, hh in (("clean", h), ("noisy", h + noise)):
            t = time.time()
            js = hes.JointContact(m, X0, hh, tg.c, hes.C_K, hes.W_COMP)
            (X, N, S), st = js.run(X0.copy(), dp.vertex_normals(m, X0), h0)
            a = report("store", X, Xt, dist, size, st, time.time() - t, js.nvar)
            t = time.time()
            cf = ContactFrames(m, X0, hh[:, 0], hh[:, 1], tg.c)
            (X, Y), st = cf.run(X0.copy())
            b = report("frames", X, Xt, dist, size, st, time.time() - t, cf.nvar)
            print(f"   contact net, {tl}\n{a}\n{b}", flush=True)
    if "weingarten" in which:
        g = sd.edge_average_readings(m, fl.build_stencils(m, ops.cr_derivative(m, Xt)[0]), Xt)
        g0 = sd.edge_average_readings(m, fl.build_stencils(m, ops.cr_derivative(m, X0)[0]), X0)
        noise = 0.1 * np.sqrt(np.mean(g ** 2)) * rng.normal(size=g.shape)
        for tl, gg in (("clean", g), ("noisy", g + noise)):
            t = time.time()
            jw = hes.JointWeingarten(m, X0, gg, hes.C_K, hes.W_COMP)
            (X, S), st = jw.run(X0.copy(), g0)
            a = report("store", X, Xt, dist, size, st, time.time() - t, 3 * m.n + 3 * len(m.E))
            t = time.time()
            wf = WeingartenFrames(m, X0, gg, tg.c)
            (X, Y), st = wf.run(X0.copy())
            b = report("frames", X, Xt, dist, size, st, time.time() - t, wf.nvar)
            print(f"   weingarten, {tl}\n{a}\n{b}", flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    which = [a for a in args if a in ("contact", "weingarten")] or ["contact", "weingarten"]
    only = [a for a in args if a not in ("contact", "weingarten")]
    rng = np.random.default_rng(7)
    t0 = time.time()
    from lie_compare import cases
    for name, X0, F, Xt, dist in cases():
        if only and not any(o in name for o in only):
            continue
        run_case(name, X0, F, Xt, dist, rng, which)
    print(f"({time.time() - t0:.0f}s)")
