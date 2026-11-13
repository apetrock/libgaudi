"""The vertex (H, z) store under the contact net (MATH.md 19.21).

1. `ContactNetGN` with its `(k, tau)` targets raw against projected onto the
   vertex store (the edge tensor `(a, tau, b)` of `lie_edge_net.read_spheres`).
2. A joint solve with the store as unknowns, as the demo has it (positions,
   normals and store in one factor): Rodrigues per edge on the store's edge
   values (`dn = -(H dx + z conj(dx))` along the edge: `k = k11`, `tau = k12`),
   tangency, unit normals, store target rows, the fading length memory.
   Edge store (3 per edge) against vertex store (3 per vertex, through `B`):
   unknowns, factor fill, factorization time, accuracy.
"""

import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import contact_net as cn
import cr_ops as ops
import dupin_poisson as dp
import hz_store as hz
import lie_edges as le
import lie_edge_net as lsn
from lie_compare import cases
from live_compare import shape_error


def _nrm(v):
    return v / np.linalg.norm(v, axis=1)[:, None]


class JointStore:
    """Unknowns `[x (3n), n (3n), s]`, `s` the edge store (3 nE) or the vertex
    store (3 nV, mapped to edges by `vs.B`)."""

    def __init__(self, m, X0, h, c, vs=None, w_target=10.0, w_memory=0.01, fd=1e-7):
        self.m, self.h, self.c, self.vs = m, h, c, vs
        self.w_t, self.w_m, self.fd = w_target, w_memory, fd
        ring = dp.rings(m)
        self.pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        self.hl = float(np.linalg.norm(X0[m.E[:, 1]] - X0[m.E[:, 0]], axis=1).mean())
        self.ns = 3 * (m.n if vs is not None else len(m.E))
        self.nvar = 6 * m.n + self.ns
        self.memory_scale = 1.0
        if vs is not None:
            # Per edge the phase factors of its two endpoints.
            self.cph = np.exp(-2j * vs.phi)

    def edge_kt(self, S, e_idx, V):
        """`(k, tau)` per edge from the store columns in `V` (edge: 3, vertex: 6)."""
        if self.vs is None:
            return V[:, 0], V[:, 1]
        c = self.cph[e_idx]
        zi = (V[:, 1] + 1j * V[:, 2]) * c[:, 0]
        zj = (V[:, 4] + 1j * V[:, 5]) * c[:, 1]
        H = 0.5 * (V[:, 0] + V[:, 3])
        z = 0.5 * (zi + zj)
        return H + z.real, z.imag

    def blocks(self, X, N, S, L):
        m, P = self.m, self.pairs
        i, j = m.E[:, 0], m.E[:, 1]
        e = np.arange(len(m.E))
        xo = lambda v: 3 * v
        no = lambda v: 3 * m.n + 3 * v
        so = lambda v: 6 * m.n + 3 * v
        if self.vs is None:
            Ve, starts_s = S[e], [so(e)]
        else:
            Ve, starts_s = np.concatenate([S[i], S[j]], 1), [so(i), so(j)]
        ns_cols = Ve.shape[1]

        def rod(V):
            xi, xj, ni, nj = V[:, 0:3], V[:, 3:6], V[:, 6:9], V[:, 9:12]
            k, tau = self.edge_kt(S, e, V[:, 12:12 + ns_cols])
            dx, dn = xj - xi, nj - ni
            nb = _nrm(ni + nj)
            dxt = dx - np.einsum("ij,ij->i", dx, nb)[:, None] * nb
            return k[:, None] * dxt + dn + tau[:, None] * np.cross(nb, dx)

        def target(V):
            if self.vs is None:
                return self.hl * (V - self.h)
            k11, k12 = self.edge_kt(S, e, V)
            c = self.cph
            zi = (V[:, 1] + 1j * V[:, 2]) * c[:, 0]
            zj = (V[:, 4] + 1j * V[:, 5]) * c[:, 1]
            H = 0.5 * (V[:, 0] + V[:, 3])
            z = 0.5 * (zi + zj)
            return self.hl * (np.column_stack([H + z.real, z.imag, H - z.real]) - self.h)

        def tang(V):
            d = V[:, 3:6] - V[:, 0:3]
            return (np.einsum("ij,ij->i", _nrm(V[:, 6:9]), d) / np.linalg.norm(d, axis=1) - self.c)[:, None]

        return [
            (1.0, rod, np.concatenate([X[i], X[j], N[i], N[j], Ve], 1), [xo(i), xo(j), no(i), no(j)] + starts_s),
            (1.0, tang, np.concatenate([X[P[:, 0]], X[P[:, 1]], N[P[:, 0]]], 1), [xo(P[:, 0]), xo(P[:, 1]), no(P[:, 0])]),
            (1.0, lambda V: (np.linalg.norm(V, axis=1) - 1.0)[:, None], N.copy(), [no(np.arange(m.n))]),
            (np.sqrt(self.w_t), target, Ve.copy(), starts_s),
            (np.sqrt(self.w_m * self.memory_scale) / self.hl,
             lambda V: (np.linalg.norm(V[:, 3:6] - V[:, 0:3], axis=1) - L)[:, None],
             np.concatenate([X[i], X[j]], 1), [xo(i), xo(j)]),
        ]

    def unpack(self, z):
        n = self.m.n
        return z[:3 * n].reshape(n, 3), z[3 * n:6 * n].reshape(n, 3), z[6 * n:].reshape(-1, 3)

    def residual(self, z, L):
        return np.concatenate([(w * f(V)).ravel() for w, f, V, _ in self.blocks(*self.unpack(z), L)])

    def jacobian(self, z, L):
        rows, cols, vals = [], [], []
        r0 = 0
        for w, f, V, starts in self.blocks(*self.unpack(z), L):
            base = f(V)
            nr, nc = base.shape
            for s, start in enumerate(starts):
                for c in range(3):
                    col = 3 * s + c
                    Vp, Vm = V.copy(), V.copy()
                    Vp[:, col] += self.fd
                    Vm[:, col] -= self.fd
                    D = w * (f(Vp) - f(Vm)) / (2 * self.fd)
                    for q in range(nc):
                        rows.append(r0 + np.arange(nr) * nc + q)
                        cols.append(start + c)
                        vals.append(D[:, q])
            r0 += nr * nc
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(r0, self.nvar))

    def run(self, X, N, S0, iterations=40, lam0=1e-4):
        m, E = self.m, self.m.E
        z = np.concatenate([X.ravel(), N.ravel(), S0.ravel()])
        L = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1)
        r = self.residual(z, L)
        cost0 = 0.5 * r @ r
        lam = lam0
        stats = dict(factorizations=0, t_factor=0.0, fill=0, its=0)
        for it in range(iterations):
            Xc = self.unpack(z)[0]
            L = np.linalg.norm(Xc[E[:, 1]] - Xc[E[:, 0]], axis=1)
            r = self.residual(z, L)
            cost = 0.5 * r @ r
            if getattr(self, "fade", True):
                self.memory_scale = min(1.0, cost / max(cost0, 1e-300))
            J = self.jacobian(z, L)
            H = (J.T @ J).tocsc()
            g = J.T @ r
            ok = False
            for _ in range(8):
                dg = H.diagonal()
                t = time.time()
                lu = spla.splu((H + sp.diags(lam * dg + 1e-12 * dg.mean())).tocsc())
                stats["t_factor"] += time.time() - t
                stats["factorizations"] += 1
                stats["fill"] = lu.L.nnz + lu.U.nnz
                dz = -lu.solve(g)
                rn = self.residual(z + dz, L)
                if np.all(np.isfinite(rn)) and 0.5 * rn @ rn < cost:
                    ok = True
                    break
                lam *= 4
            if not ok:
                break
            z = z + dz
            lam = max(lam / 3, 1e-12)
            stats["its"] = it + 1
            if np.abs(dz[:3 * m.n]).max() < 1e-12:
                break
        return self.unpack(z), stats


if __name__ == "__main__":
    only = [a for a in sys.argv[1:] if not a.startswith("-")]
    print("1. ContactNetGN targets: raw edge read-out against projected onto the vertex store")
    for name, X0, F, Xt, dist in cases():
        if only and not any(o in name for o in only):
            continue
        m = ops.build_mesh(F, len(X0))
        size = np.ptp(Xt, axis=0).max()
        ring = dp.rings(m)
        pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        tg = cn.read_targets(m, pairs, Xt)
        h = lsn.read_spheres(m, le.hinges(m), Xt).H
        vs = hz.VertexStore(m, Xt)
        hp = vs.project(h)
        out = []
        for label, k, tau in (("raw", tg.k, tg.tau), ("vertex", hp[:, 0], hp[:, 1])):
            g = cn.ContactNetGN(m, X0, cn.NetTargets(k, tau, tg.c), cn.NetConfig(), pairs)
            r = g.run(X0.copy(), iterations=60)
            out.append(f"{label} {shape_error(r.X, Xt, dist, size):.1e} ({len(r.history) - 1} its, "
                       f"{r.factorizations} facts)")
        print(f"  {name:22s} " + "   ".join(out), flush=True)

    print("2. Joint solve with the store as unknowns: edge store against vertex store")
    for name, X0, F, Xt, dist in cases():
        if only and not any(o in name for o in only):
            continue
        m = ops.build_mesh(F, len(X0))
        size = np.ptp(Xt, axis=0).max()
        ring = dp.rings(m)
        pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        tg = cn.read_targets(m, pairs, Xt)
        h = lsn.read_spheres(m, le.hinges(m), Xt).H
        h[:, 0], h[:, 1] = tg.k, tg.tau
        N0 = dp.vertex_normals(m, X0)
        h0 = lsn.read_spheres(m, le.hinges(m), X0).H
        res = []
        for label, vsx in (("edge", None), ("vertex", hz.VertexStore(m, X0))):
            js = JointStore(m, X0, h, tg.c, vsx)
            S0 = h0 if vsx is None else vsx.fit(h0)
            t = time.time()
            (X, N, S), st = js.run(X0.copy(), N0.copy(), S0)
            res.append(f"{label}: unknowns {js.nvar} shape {shape_error(X, Xt, dist, size):.1e} its {st['its']} "
                       f"facts {st['factorizations']} fill {st['fill']} factor {st['t_factor']:.2f}s "
                       f"total {time.time() - t:.1f}s")
        print(f"  {name:22s}\n    " + "\n    ".join(res), flush=True)
