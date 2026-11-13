"""Lie edge net: two Lie spheres per edge, composed at the vertices
(MATH.md 19.14).

State: one contact element per vertex (positions `x` and normals `n`, both
unknown, as in `contact_net.py`), and per edge its two spheres on the edge's
contact element plus the twist, `H_e = (a, tau, b)`:

- `a`: the along sphere `a P + T`, the edge's normal curvature along itself
  (the contact sphere of its endpoints when `tau = 0`);
- `b`: the across sphere `b P + T`, the curvature across the edge (the hinge
  sphere, tangent to both face planes);
- `tau`: the twist (geodesic torsion), 0 exactly on curvature-line edges.

That is the CR edge store `(k11, k12, k22)` in Lie form. On a triangle mesh
most edges are not curvature lines, so `a` and `b` are Euler samples, not
curvature spheres: the edge needs `tau` too.

Vertices compose. Each vertex holds an auxiliary shape tensor `S_i` (in a
tangent frame built from its live normal); every incident edge's spheres
must be its Euler samples, `kappa(phi) = cos^2 phi k_1 + sin^2 phi k_2`
(on the contact line, `s(phi) = cos^2 phi s_1 + sin^2 phi s_2`: linear in
R^{4,2}). Coupling acts on the composed `S` (the curvature spheres), never
sphere to sphere across edges of different directions.

Rows (Levenberg-Marquardt, finite-difference Jacobians per row block, the
`ContactNetGN` pattern):

- Rodrigues per edge: `a dx_t + dn + tau (nbar x dx) = 0` (pins the edge's
  direction; the sphere rows alone let a neighbour slide on its sphere);
- Lie incidence per directed spoke: `<a P(x_i) + T(x_i, n_i), P(x_j)> = 0`,
  the neighbour on the along sphere at `x_i` (both ends: their sum is the
  along part of Rodrigues, their difference the symmetric tangency);
- hinge per interior edge: `2 tan(theta / 2) = b l*` (the across sphere
  tangent to both face planes, `l*` the circumcentric dual length);
- vertex Euler consistency per (vertex, incident edge);
- Poisson coupling of the composed tensors across each edge (ambient, so no
  frame gauge);
- unit normals; targets / edits on `H`; the fading length memory.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import contact_net as cn
import cr_ops as ops
import dupin_poisson as dp
import lie_edges as le

# --------------------------------------------------------------- Lie rows


def lie_dot_rows(a, b):
    """Row-wise `<a, b>` in R^{4,2}, eta = diag(-1, 1, 1, 1, 1, -1)."""
    return -a[:, 0] * b[:, 0] + np.einsum("ij,ij->i", a[:, 1:5], b[:, 1:5]) - a[:, 5] * b[:, 5]


def point_rows(x):
    """`P(x) = ((1 + |x|^2) / 2, (1 - |x|^2) / 2, x, 0)` per row."""
    x2 = np.einsum("ij,ij->i", x, x)
    return np.column_stack([(1 + x2) / 2, (1 - x2) / 2, x, np.zeros(len(x))])


def plane_rows(n, h):
    """`T(n, h) = (h, -h, n, 1)` per row (oriented plane `n . x = h`)."""
    return np.column_stack([h, -h, n, np.ones(len(n))])


def curvature_sphere_rows(p, n, k):
    """`k P(p) + T(n, n . p)`: the sphere of curvature `k` tangent at `(p, n)`
    (centre `p + n / k`; a plane at `k = 0`)."""
    return k[:, None] * point_rows(p) + plane_rows(n, np.einsum("ij,ij->i", n, p))


def _nrm(v):
    return v / np.linalg.norm(v, axis=1)[:, None]


def vertex_frames(N, R):
    """Tangent frame `(T, B)` at each live normal from a fixed reference
    direction `R` (chosen at the start, least aligned with the normal)."""
    n = _nrm(N)
    T = _nrm(R - np.einsum("ij,ij->i", R, n)[:, None] * n)
    return T, np.cross(n, T), n


def reference_dirs(N):
    I = np.eye(3)
    return I[np.argmin(np.abs(_nrm(N)), axis=1)]


def tensor_in_dir(S, u, w):
    """`u^T S w` for 2-vectors `u`, `w` and `S = (s11, s12, s22)` per row."""
    return S[:, 0] * u[:, 0] * w[:, 0] + S[:, 1] * (u[:, 0] * w[:, 1] + u[:, 1] * w[:, 0]) + S[:, 2] * u[:, 1] * w[:, 1]


# ------------------------------------------------------------- row blocks


def _rodrigues(V):
    xi, xj, ni, nj, H = V[:, 0:3], V[:, 3:6], V[:, 6:9], V[:, 9:12], V[:, 12:15]
    dx, dn = xj - xi, nj - ni
    nb = _nrm(ni + nj)
    dxt = dx - np.einsum("ij,ij->i", dx, nb)[:, None] * nb
    return H[:, 0:1] * dxt + dn + H[:, 1:2] * np.cross(nb, dx)


def _incidence(V):
    xi, xj, ni, H = V[:, 0:3], V[:, 3:6], V[:, 6:9], V[:, 9:12]
    s = curvature_sphere_rows(xi, _nrm(ni), H[:, 0])
    d = np.linalg.norm(xj - xi, axis=1)
    return (lie_dot_rows(s, point_rows(xj)) / d)[:, None]


def _hinge(V):
    """`2 tan(theta / 2) - b l*`: `a -> b` along face f, `c` f's third
    vertex, `d` the other face's; the sign convention of `le.hinge_angles`."""
    a, b, c, d, H = V[:, 0:3], V[:, 3:6], V[:, 6:9], V[:, 9:12], V[:, 12:15]
    nf = _nrm(np.cross(b - a, c - a))
    ng = _nrm(np.cross(a - b, d - b))
    t = _nrm(b - a)
    th = -np.arctan2(np.einsum("ij,ij->i", t, np.cross(nf, ng)), np.einsum("ij,ij->i", nf, ng))
    le_ = np.linalg.norm(b - a, axis=1)

    def cot(p, q, o):
        u, v = p - o, q - o
        return np.einsum("ij,ij->i", u, v) / np.linalg.norm(np.cross(u, v), axis=1)

    ls = 0.5 * le_ * (cot(a, b, c) + cot(a, b, d))
    return (2 * np.tan(0.5 * th) - H[:, 2] * ls)[:, None]


def _vertex(V, R, h):
    """Edge `(a, tau, b)` against vertex `i`'s tensor sampled along and across
    the spoke `i -> j` (projected into `i`'s tangent plane)."""
    xi, xj, ni, H, S = V[:, 0:3], V[:, 3:6], V[:, 6:9], V[:, 9:12], V[:, 12:15]
    T, B, n = vertex_frames(ni, R)
    d = xj - xi
    t = _nrm(d - np.einsum("ij,ij->i", d, n)[:, None] * n)
    u = np.column_stack([np.einsum("ij,ij->i", t, T), np.einsum("ij,ij->i", t, B)])
    w = np.column_stack([-u[:, 1], u[:, 0]])
    k = np.column_stack([tensor_in_dir(S, u, u), tensor_in_dir(S, u, w), tensor_in_dir(S, w, w)])
    return h * (k - H)


def _ambient(S, T, B):
    """The tangent tensor as a symmetric 3 x 3, upper entries (off-diagonals
    times sqrt 2, so the row norm is the Frobenius norm)."""
    F = np.stack([T, B], 2)  # (m, 3, 2)
    S2 = np.stack([np.stack([S[:, 0], S[:, 1]], 1), np.stack([S[:, 1], S[:, 2]], 1)], 1)
    A = F @ S2 @ np.transpose(F, (0, 2, 1))
    r2 = np.sqrt(2.0)
    return np.column_stack([A[:, 0, 0], A[:, 1, 1], A[:, 2, 2], r2 * A[:, 0, 1], r2 * A[:, 0, 2], r2 * A[:, 1, 2]])


def _couple(V, Ri, Rj, h):
    ni, nj, Si, Sj = V[:, 0:3], V[:, 3:6], V[:, 6:9], V[:, 9:12]
    Ti, Bi, _ = vertex_frames(ni, Ri)
    Tj, Bj, _ = vertex_frames(nj, Rj)
    return h * (_ambient(Si, Ti, Bi) - _ambient(Sj, Tj, Bj))


def _unit(V):
    return (np.linalg.norm(V, axis=1) - 1.0)[:, None]


def _length(V, L):
    return (np.linalg.norm(V[:, 3:6] - V[:, 0:3], axis=1) - L)[:, None]


def _target(V, Ht, Wt, h):
    return h * Wt * (V - Ht)


# ---------------------------------------------------------------- targets


@dataclass
class EdgeSpheres:
    H: np.ndarray  # (nE, 3): a (along), tau (twist), b (across)
    W: np.ndarray  # (nE, 3): per-component target weights (0 where unknown)


def read_spheres(m: ops.Mesh, hin: le.Hinges, X, N=None) -> EdgeSpheres:
    """`a`, `tau` from the contact read-out (`contact_net.edge_readout`), `b`
    from the hinge (interior edges, weighted by the hinge's reliability; 0 on
    boundary edges, which the vertices compose)."""
    N = dp.vertex_normals(m, X) if N is None else N
    k, tau = cn.edge_readout(m.E, X, N)
    b = np.zeros(len(m.E))
    wb = np.zeros(len(m.E))
    kp, _ = le.hinge_curvature(hin, X, m=m)
    w = le.hinge_weight(hin, X)
    b[hin.e] = kp
    wb[hin.e] = np.sqrt(w / max(w.mean(), 1e-300))
    return EdgeSpheres(np.column_stack([k, tau, b]), np.column_stack([np.ones_like(k), np.ones_like(k), wb]))


def fit_vertex_tensors(m: ops.Mesh, X, N, H, R):
    """Per vertex the least-squares tensor whose Euler samples best match its
    incident edges' `(a, tau, b)` (the initial `S`)."""
    T, B, n = vertex_frames(N, R)
    S = np.zeros((m.n, 3))
    rows = [[] for _ in range(m.n)]
    for e, (i, j) in enumerate(m.E):
        rows[i].append((e, j))
        rows[j].append((e, i))
    for i in range(m.n):
        A, y = [], []
        for e, j in rows[i]:
            d = X[j] - X[i]
            t = d - (d @ n[i]) * n[i]
            t /= np.linalg.norm(t)
            u = np.array([t @ T[i], t @ B[i]])
            w = np.array([-u[1], u[0]])
            for p, q, val in ((u, u, H[e, 0]), (u, w, H[e, 1]), (w, w, H[e, 2])):
                A.append([p[0] * q[0], p[0] * q[1] + p[1] * q[0], p[1] * q[1]])
                y.append(val)
        S[i] = np.linalg.lstsq(np.array(A), np.array(y), rcond=None)[0]
    return S


# ---------------------------------------------------------------- solver


@dataclass
class EdgeNetConfig:
    iterations: int = 200
    w_rodrigues: float = 1.0
    w_incidence: float = 1.0
    w_hinge: float = 1.0
    w_vertex: float = 1.0
    w_couple: float = 1e-3
    w_unit: float = 1.0
    w_target: float = 10.0
    # Not faded: with targets that are not exactly consistent (the read-outs'
    # biases), a faded memory lets open sheets stretch instead of bending.
    w_memory: float = 0.1
    memory_fade: bool = False
    # Incidence and vertex rows at boundary vertices (incomplete rings,
    # one-sided normals, no across sphere on boundary edges) scaled by this.
    w_boundary: float = 0.1
    gn_refresh: int = 5
    lam0: float = 1e-4
    fd: float = 1e-7
    tol: float = 1e-12


@dataclass
class EdgeNetResult:
    X: np.ndarray
    N: np.ndarray
    H: np.ndarray
    S: np.ndarray
    history: list = field(default_factory=list)  # (move, merit)
    seconds: float = 0.0
    best_merit: float = np.inf
    factorizations: int = 0
    backsolves: int = 0


class LieEdgeNet:
    """Unknowns `z = [x (3n), n (3n), H (3 nE), S (3n)]`. `targets` holds the
    edge spheres to reach (read off a target mesh, or edited)."""

    def __init__(self, m: ops.Mesh, X0, targets: EdgeSpheres, cfg: EdgeNetConfig = EdgeNetConfig(), N0=None):
        self.m, self.cfg, self.tgt = m, cfg, targets
        self.hin = le.hinges(m)
        N0 = dp.vertex_normals(m, X0) if N0 is None else N0
        self.R = reference_dirs(N0)
        self.h = float(np.linalg.norm(X0[m.E[:, 1]] - X0[m.E[:, 0]], axis=1).mean())
        i, j = m.E[:, 0], m.E[:, 1]
        e = np.arange(len(m.E))
        self.spoke = np.column_stack([np.concatenate([i, j]), np.concatenate([j, i]), np.concatenate([e, e])])
        bnd = np.zeros(m.n, bool)
        for k, fs in enumerate(m.edge_faces):
            if len(fs) == 1:
                bnd[m.E[k]] = True
        self.spoke_w = np.where(bnd[self.spoke[:, 0]], np.sqrt(cfg.w_boundary), 1.0)[:, None]
        self.memory_scale = 1.0
        nE = len(m.E)
        self.off = {"x": 0, "n": 3 * m.n, "H": 6 * m.n, "S": 6 * m.n + 3 * nE}
        self.nvar = 9 * m.n + 3 * nE

    # -- packing
    def pack(self, X, N, H, S):
        return np.concatenate([X.ravel(), N.ravel(), H.ravel(), S.ravel()])

    def unpack(self, z):
        n, nE, o = self.m.n, len(self.m.E), self.off
        return (z[:3 * n].reshape(n, 3), z[o["n"]:o["n"] + 3 * n].reshape(n, 3),
                z[o["H"]:o["H"] + 3 * nE].reshape(nE, 3), z[o["S"]:].reshape(n, 3))

    def _blocks(self, X, N, H, S, L):
        m, cfg, o, hin, sp_ = self.m, self.cfg, self.off, self.hin, self.spoke
        xo = lambda v: o["x"] + 3 * v
        no = lambda v: o["n"] + 3 * v
        ho = lambda e: o["H"] + 3 * e
        so = lambda v: o["S"] + 3 * v
        i, j = m.E[:, 0], m.E[:, 1]
        si, sj, se = sp_[:, 0], sp_[:, 1], sp_[:, 2]
        h, R = self.h, self.R
        blocks = [
            (np.sqrt(cfg.w_rodrigues), _rodrigues,
             np.concatenate([X[i], X[j], N[i], N[j], H], 1), [xo(i), xo(j), no(i), no(j), ho(np.arange(len(m.E)))]),
            (np.sqrt(cfg.w_incidence), lambda V: self.spoke_w * _incidence(V),
             np.concatenate([X[si], X[sj], N[si], H[se]], 1), [xo(si), xo(sj), no(si), ho(se)]),
            (np.sqrt(cfg.w_vertex), lambda V: self.spoke_w * _vertex(V, R[si], h),
             np.concatenate([X[si], X[sj], N[si], H[se], S[si]], 1), [xo(si), xo(sj), no(si), ho(se), so(si)]),
            (np.sqrt(cfg.w_couple), lambda V: _couple(V, R[i], R[j], h),
             np.concatenate([N[i], N[j], S[i], S[j]], 1), [no(i), no(j), so(i), so(j)]),
            (np.sqrt(cfg.w_unit), _unit, N.copy(), [no(np.arange(m.n))]),
            (np.sqrt(cfg.w_target), lambda V: _target(V, self.tgt.H, self.tgt.W, h), H.copy(),
             [ho(np.arange(len(m.E)))]),
            (np.sqrt(cfg.w_memory * self.memory_scale) / h, lambda V: _length(V, L),
             np.concatenate([X[i], X[j]], 1), [xo(i), xo(j)]),
        ]
        if cfg.w_hinge > 0 and len(hin.e):
            wh = np.sqrt(np.maximum(self.tgt.W[hin.e, 2], 0.0))
            blocks.append((np.sqrt(cfg.w_hinge), lambda V: wh[:, None] * _hinge(V),
                           np.concatenate([X[hin.a], X[hin.b], X[hin.c], X[hin.d], H[hin.e]], 1),
                           [xo(hin.a), xo(hin.b), xo(hin.c), xo(hin.d), ho(hin.e)]))
        return blocks

    def residual(self, z, L):
        return np.concatenate([(w * f(V)).ravel() for w, f, V, _ in self._blocks(*self.unpack(z), L)])

    def jacobian(self, z, L):
        rows, cols, vals = [], [], []
        r0 = 0
        fd = self.cfg.fd
        for w, f, V, starts in self._blocks(*self.unpack(z), L):
            base = f(V)
            nr, nc = base.shape
            for s, start in enumerate(starts):
                for c in range(3):
                    col = 3 * s + c
                    Vp, Vm = V.copy(), V.copy()
                    Vp[:, col] += fd
                    Vm[:, col] -= fd
                    D = w * (f(Vp) - f(Vm)) / (2 * fd)
                    for q in range(nc):
                        rows.append(r0 + np.arange(nr) * nc + q)
                        cols.append(start + c)
                        vals.append(D[:, q])
            r0 += nr * nc
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(r0, self.nvar))

    def merit(self, z, L):
        r = self.residual(z, L)
        return 0.5 * float(r @ r)

    def initial(self, X, N=None):
        """`H` read off the start mesh, `S` fitted to it."""
        N = dp.vertex_normals(self.m, X) if N is None else N
        H = read_spheres(self.m, self.hin, X, N).H
        S = fit_vertex_tensors(self.m, X, N, H, self.R)
        return self.pack(X, N, H, S)

    def run(self, X, N=None, iterations=None, keep_best=True, log=None, z0=None) -> EdgeNetResult:
        t0 = time.time()
        cfg, m, E = self.cfg, self.m, self.m.E
        iterations = cfg.iterations if iterations is None else iterations
        z = self.initial(X, N) if z0 is None else z0.copy()
        L = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1)
        self.memory_scale = 1.0
        mer0 = max(self.merit(z, L), 1e-300)
        out = EdgeNetResult(*self.unpack(z))
        out.history.append((0.0, mer0))
        best, best_m = z.copy(), mer0
        lam, solve, age = cfg.lam0, None, 0
        for it in range(iterations):
            if cfg.memory_fade and (solve is None or age >= cfg.gn_refresh):
                self.memory_scale = min(1.0, out.history[-1][1] / mer0)
            Xc = self.unpack(z)[0]
            L = np.linalg.norm(Xc[E[:, 1]] - Xc[E[:, 0]], axis=1)
            r = self.residual(z, L)
            J = self.jacobian(z, L)
            g = J.T @ r
            cost = 0.5 * r @ r
            if solve is None or age >= cfg.gn_refresh:
                Hm = (J.T @ J).tocsc()
                dg = Hm.diagonal()
                solve = spla.factorized((Hm + sp.diags(lam * dg + 1e-12 * dg.mean())).tocsc())
                out.factorizations += 1
                age = 0
            ok = False
            for _ in range(10):
                dz = -solve(g)
                out.backsolves += 1
                zn = z + dz
                rn = self.residual(zn, L)
                if np.all(np.isfinite(rn)) and 0.5 * rn @ rn < cost:
                    ok = True
                    break
                lam *= 4.0
                Hm = (J.T @ J).tocsc()
                dg = Hm.diagonal()
                solve = spla.factorized((Hm + sp.diags(lam * dg + 1e-12 * dg.mean())).tocsc())
                out.factorizations += 1
                age = 0
            if not ok:
                break
            lam = max(lam / 3.0, 1e-12)
            age += 1
            move = float(np.abs(dz[:3 * m.n]).max())
            z = zn
            Xn = self.unpack(z)[0]
            mer = self.merit(z, np.linalg.norm(Xn[E[:, 1]] - Xn[E[:, 0]], axis=1))
            out.history.append((move, mer))
            if mer < best_m:
                best, best_m = z.copy(), mer
            if log:
                log(it, Xn, move, mer)
            if move < cfg.tol:
                break
        out.best_merit = best_m
        zf = best if keep_best else z
        out.X, out.N, out.H, out.S = self.unpack(zf)
        out.seconds = time.time() - t0
        return out
