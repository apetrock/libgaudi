"""Lie group synchronization: global vertex positions from the contact-sphere
network (MATH.md 19.17).

Each vertex holds a transformation `F_i` placing a canonical contact element
(origin, normal e_z, tangent e_x toward its first fan neighbour) on its own;
each edge prescribes the relative transformation `G_ij = F_i^-1 F_j` from
the network; one sparse least-squares solve synchronizes them,

    min sum_ij w_ij |F_j - F_i G_ij|^2 + mu sum_i M_i |F_i - F_i^0|^2,

and positions and normals are read off the `F_i` (translations are group
elements, so no second Poisson for positions). Row-vector form: each row of
the 3 x 4 block `[R_i | x_i]` satisfies `f_j = f_i G_ij`, so the three rows
share one matrix: one factorization, three right-hand sides.

Stage `euclidean` (this file's first stage): `F_i` in SE(3), the Euclidean
subgroup of the Lie sphere group (4 x 4 here; the 6 x 6 form is only needed
once inversions enter). `G_ij` in the edge's midpoint frame `(t, b, nbar)`:
the chord `d = (l, 0, d_n)`, the endpoint frames turned by `-+ omega / 2`
with `omega = tau l t - a l b` (the Rodrigues relation `dn = -a dx_t -
tau nbar x dx` as a rotation of the contact element, magnitude made exact),
and each frame's tangent axis at the spoke's polar angle `theta_ij` (from
the metric: corner angles scaled to 2 pi). So `G_ij` uses only the metric
(lengths and polar angles) and the edge's along sphere and twist: one
factorization per frame.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_ops as ops
import dupin_poisson as dp
import lie_edge_net as lsn


def _rot(axis_angle):
    """Rodrigues rotation matrices for rows of rotation vectors."""
    th = np.linalg.norm(axis_angle, axis=1)
    k = axis_angle / np.maximum(th, 1e-300)[:, None]
    K = np.zeros((len(th), 3, 3))
    K[:, 0, 1], K[:, 0, 2], K[:, 1, 2] = -k[:, 2], k[:, 1], -k[:, 0]
    K[:, 1, 0], K[:, 2, 0], K[:, 2, 1] = k[:, 2], -k[:, 1], k[:, 0]
    s, c = np.sin(th)[:, None, None], np.cos(th)[:, None, None]
    return np.eye(3)[None] + s * K + (1 - c) * K @ K


def _rz(a):
    c, s = np.cos(a), np.sin(a)
    R = np.zeros((len(a), 3, 3))
    R[:, 0, 0], R[:, 0, 1], R[:, 1, 0], R[:, 1, 1], R[:, 2, 2] = c, -s, s, c, 1
    return R


def _hom(R, t):
    T = np.zeros((len(R), 4, 4))
    T[:, :3, :3] = R
    T[:, :3, 3] = t
    T[:, 3, 3] = 1
    return T


def _inv(T):
    R, t = T[:, :3, :3], T[:, :3, 3]
    Rt = np.transpose(R, (0, 2, 1))
    return _hom(Rt, -np.einsum("nij,nj->ni", Rt, t))


# ------------------------------------------------------------- topology


def vertex_fans(m: ops.Mesh):
    """Per vertex its neighbours in counter-clockwise order about the outward
    normal (face winding), and whether the fan is closed."""
    nxt = [dict() for _ in range(m.n)]
    for f in m.F:
        for k in range(3):
            i, a, b = f[k], f[(k + 1) % 3], f[(k + 2) % 3]
            nxt[i][a] = b
    fans, closed = [], []
    for i in range(m.n):
        d = nxt[i]
        if not d:
            fans.append([])
            closed.append(False)
            continue
        starts = set(d) - set(d.values())
        v = next(iter(starts)) if starts else next(iter(d))
        order = [v]
        while v in d and d[v] != order[0] and len(order) <= len(d):
            v = d[v]
            order.append(v)
        fans.append(order)
        closed.append(not starts)
    return fans, closed


def edge_lengths(m: ops.Mesh, X):
    return np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1)


def _face_edge_lengths(m: ops.Mesh, lengths):
    index = {(min(a, b), max(a, b)): e for e, (a, b) in enumerate(m.E)}
    Lf = np.zeros((len(m.F), 3))  # Lf[f, k]: the edge opposite corner k
    for f, (p, q, r) in enumerate(m.F):
        Lf[f] = [lengths[index[(min(q, r), max(q, r))]], lengths[index[(min(p, r), max(p, r))]],
                 lengths[index[(min(p, q), max(p, q))]]]
    return Lf, index


def corner_angles(m: ops.Mesh, lengths):
    """Per face the angle at each corner, from edge lengths alone."""
    Lf, _ = _face_edge_lengths(m, lengths)
    a, b, c = Lf[:, 0], Lf[:, 1], Lf[:, 2]
    A0 = np.arccos(np.clip((b * b + c * c - a * a) / (2 * b * c), -1, 1))
    A1 = np.arccos(np.clip((a * a + c * c - b * b) / (2 * a * c), -1, 1))
    return np.column_stack([A0, A1, np.pi - A0 - A1])


def intrinsic(m: ops.Mesh, lengths):
    """Angle defects (interior; boundary vertices `pi - sum`), vertex areas
    (a third of each face's Heron area), and the cotan Laplacian (positive
    semidefinite), all from edge lengths."""
    ang = corner_angles(m, lengths)
    Lf, index = _face_edge_lengths(m, lengths)
    s = 0.5 * Lf.sum(1)
    area = np.sqrt(np.maximum(s * (s - Lf[:, 0]) * (s - Lf[:, 1]) * (s - Lf[:, 2]), 0.0))
    total = np.zeros(m.n)
    A = np.zeros(m.n)
    for k in range(3):
        np.add.at(total, m.F[:, k], ang[:, k])
        np.add.at(A, m.F[:, k], area / 3)
    bnd = np.zeros(m.n, bool)
    for e, fs in enumerate(m.edge_faces):
        if len(fs) == 1:
            bnd[m.E[e]] = True
    omega = np.where(bnd, np.pi, 2 * np.pi) - total
    w = np.zeros(len(m.E))
    for f, (p, q, r) in enumerate(m.F):
        for k, (u, v) in enumerate(((q, r), (p, r), (p, q))):
            w[index[(min(u, v), max(u, v))]] += 0.5 / np.tan(max(ang[f, k], 1e-6))
    i, j = m.E[:, 0], m.E[:, 1]
    L = sp.csr_matrix((np.concatenate([-w, -w, w, w]),
                       (np.concatenate([i, j, i, j]), np.concatenate([j, i, i, j]))), shape=(m.n, m.n))
    return omega, A, L, bnd


def polar_angles(m: ops.Mesh, lengths, fans, closed):
    """`theta[(i, j)]`: spoke `j`'s angle in `i`'s tangent frame from the
    metric (cumulative corner angles, scaled to 2 pi on closed fans), from
    edge lengths alone."""
    ang_f = corner_angles(m, lengths)
    corner = {}
    for f, tri in enumerate(m.F):
        for k in range(3):
            corner[(tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3])] = ang_f[f, k]
    theta = {}
    for i, (fan, cl) in enumerate(zip(fans, closed)):
        if not fan:
            continue
        pairs = list(zip(fan, fan[1:] + ([fan[0]] if cl else [])))
        ang = [corner.get((i, a, b), 0.0) for a, b in pairs]
        s = 2 * np.pi / sum(ang) if cl and sum(ang) > 0 else 1.0
        acc = 0.0
        for q, v in enumerate(fan):
            theta[(i, v)] = acc
            if q < len(ang):
                acc += s * ang[q]
    return theta


def initial_frames(m: ops.Mesh, X, fans):
    """`F_i^0`: origin `x_i`, normal the vertex normal, tangent axis toward the
    first fan neighbour (the polar angles' zero)."""
    N = dp.vertex_normals(m, X)
    R = np.zeros((m.n, 3, 3))
    for i, fan in enumerate(fans):
        z = N[i]
        d = X[fan[0]] - X[i] if fan else np.array([1.0, 0, 0])
        x = d - (d @ z) * z
        x /= max(np.linalg.norm(x), 1e-300)
        R[i] = np.column_stack([x, np.cross(z, x), z])
    return _hom(R, X)


# ------------------------------------------------------- relative transforms


def relative_transforms(m: ops.Mesh, lengths, H, theta, d_n=None):
    """`G_e = F_i^-1 F_j` per edge `(i, j) = E[e]` from the metric (edge
    lengths and polar angles) and the edge's along sphere `a = H[:, 0]` and
    twist `tau = H[:, 1]`."""
    i, j = m.E[:, 0], m.E[:, 1]
    l = lengths
    dn = np.zeros_like(l) if d_n is None else d_n
    lt = np.sqrt(np.maximum(l * l - dn * dn, 0.0))
    a, tau = H[:, 0], H[:, 1]
    om = np.column_stack([tau * lt, -a * lt, np.zeros_like(lt)])
    small = np.linalg.norm(om, axis=1)
    # Exact magnitude: |dn| = 2 sin(|omega| / 2).
    phi = 2 * np.arcsin(np.clip(0.5 * small, 0, 1))
    om *= (phi / np.maximum(small, 1e-300))[:, None]
    th_ij = np.array([theta.get((p, q), 0.0) for p, q in zip(i, j)])
    th_ji = np.array([theta.get((q, p), 0.0) for p, q in zip(i, j)])
    half = np.column_stack([0.5 * lt, np.zeros_like(lt), 0.5 * dn])
    Li = _hom(_rot(-0.5 * om) @ _rz(-th_ij), -half)
    Lj = _hom(_rot(0.5 * om) @ _rz(-(th_ji + np.pi)), half)
    return np.einsum("nij,njk->nik", _inv(Li), Lj)


# ---------------------------------------------------------------- solver


@dataclass
class SyncConfig:
    stage: str = "euclidean"
    mu: float = 1e-6  # pull toward the initial frames, relative
    pin: int = 0  # vertex whose frame is held exactly (the gauge)
    translation_resolve: bool = True  # re-solve positions with the projected rotations
    conformal_iterations: int = 5  # stage similarity: scale-factor updates
    conformal_step: float = 1.0
    weights: str = "uniform"  # uniform | length


@dataclass
class SyncResult:
    X: np.ndarray
    N: np.ndarray
    frames: np.ndarray
    residual: float = 0.0
    factorizations: int = 0
    backsolves: int = 0
    seconds: float = 0.0


class LieSync:
    def __init__(self, m: ops.Mesh, cfg: SyncConfig = SyncConfig()):
        self.m, self.cfg = m, cfg
        self.fans, self.closed = vertex_fans(m)

    def solve(self, X0, H, Xmetric=None) -> SyncResult:
        """`X0`: the start (initial frames, gauge); `H`: edge spheres
        `(a, tau, b)`; `Xmetric`: the shape whose lengths and angles are kept
        (default `X0`, the frame-start memory)."""
        t0 = time.time()
        m = self.m
        lengths = edge_lengths(m, X0 if Xmetric is None else Xmetric)
        facts = backs = 0
        X = X0
        if self.cfg.stage == "euclidean":
            r = self._sync(X0, H, lengths)
        else:
            # Conformal stretch: the target Gauss curvature per vertex from the
            # target spheres' vertex tensors; a scale factor u per vertex from
            # the linearized discrete-conformal equation (angle defects grow
            # with u: d Omega / d u = L, Springborn-Schroeder-Pinkall 2008),
            # Omega + L u = K_t A, iterated (Picard on the area factor):
            # (L + delta M) u = K_t A - Omega; lengths times
            # exp((u_i + u_j) / 2); re-synchronize.
            for it in range(self.cfg.conformal_iterations):
                r = self._sync(X0, H, lengths)
                facts += r.factorizations
                backs += r.backsolves
                N = r.N
                R = lsn.reference_dirs(N)
                S = lsn.fit_vertex_tensors(m, r.X, N, H, R)
                Kt = S[:, 0] * S[:, 2] - S[:, 1] ** 2
                omega, A, L, bnd = intrinsic(m, lengths)
                Kt = np.where(bnd, omega / np.maximum(A, 1e-300), Kt)  # boundary: keep
                delta = 1e-3 * L.diagonal().mean() / max(A.mean(), 1e-300)
                M_ = (L + sp.diags(delta * A)).tocsc()
                u = spla.spsolve(M_, Kt * A - omega)
                facts += 1
                u = np.clip(self.cfg.conformal_step * u, -0.5, 0.5)
                lengths = lengths * np.exp(0.5 * (u[m.E[:, 0]] + u[m.E[:, 1]]))
            r = self._sync(X0, H, lengths)
        r.factorizations += facts
        r.backsolves += backs
        r.seconds = time.time() - t0
        return r

    def _sync(self, X0, H, lengths) -> SyncResult:
        t0 = time.time()
        m = self.m
        h = float(lengths.mean())
        theta = polar_angles(m, lengths, self.fans, self.closed)
        G = relative_transforms(m, lengths / h, H * h, theta)  # unit edge length
        F0 = initial_frames(m, X0, self.fans)
        F0s = F0.copy()
        F0s[:, :3, 3] /= h
        n, nE = m.n, len(m.E)
        w = np.ones(nE) if self.cfg.weights == "uniform" else lengths / h
        # Residual (column form) per edge: f_j^T - G^T f_i^T, 4 rows.
        rows, cols, vals = [], [], []
        for e in range(nE):
            p, q = m.E[e]
            sw = np.sqrt(w[e])
            Gt = G[e].T
            for r in range(4):
                rows.append(4 * e + r)
                cols.append(4 * q + r)
                vals.append(sw)
                for c in range(4):
                    if Gt[r, c] != 0:
                        rows.append(4 * e + r)
                        cols.append(4 * p + c)
                        vals.append(-sw * Gt[r, c])
        Aop = sp.csr_matrix((vals, (rows, cols)), shape=(4 * nE, 4 * n))
        # The residual is homogeneous (all frames to zero minimizes it), so a
        # weak pull lets the frames shrink: pin one vertex's frame exactly
        # (its columns move to the right-hand side) and keep only a small pull
        # elsewhere.
        pin = self.cfg.pin
        keep = np.ones(4 * n, bool)
        keep[4 * pin:4 * pin + 4] = False
        Ak, Ap = Aop[:, keep].tocsc(), Aop[:, ~keep].tocsc()
        K = (Ak.T @ Ak).tocsc()
        mu = self.cfg.mu * K.diagonal().mean()
        solve = spla.factorized((K + mu * sp.eye(K.shape[0])).tocsc())
        # Three rows of [R | x], each a 4-vector per vertex.
        Fs = np.zeros((n, 3, 4))
        for r in range(3):
            fp = F0s[pin, r, :]
            rhs = -Ak.T @ (Ap @ fp) + mu * F0s[:, r, :].ravel()[keep]
            full = np.zeros(4 * n)
            full[keep] = solve(rhs)
            full[~keep] = fp
            Fs[:, r, :] = full.reshape(n, 4)
        res = 0.0
        for r in range(3):
            v = Aop @ Fs[:, r, :].ravel()
            res += float(v @ v)
        # Project the rotation part onto SO(3).
        R = Fs[:, :, :3]
        U, _, Vt = np.linalg.svd(R)
        D = np.ones((n, 3))
        D[:, 2] = np.sign(np.linalg.det(U @ Vt))
        Rp = np.einsum("nij,nj,njk->nik", U, D, Vt)
        X = Fs[:, :, 3] * h
        facts, backs = 1, 3
        if self.cfg.translation_resolve:
            # With the projected rotations, the translations are a Laplacian:
            # min sum w |x_j - x_i - R_i t_ij|^2, the pinned vertex held (the
            # linear relaxation shrinks the rotation blocks, and the
            # translations with them).
            i, j = m.E[:, 0], m.E[:, 1]
            d = np.einsum("nij,nj->ni", Rp[i], G[:, :3, 3]) * h
            r_ = np.arange(nE)
            D_ = sp.csr_matrix((np.concatenate([-np.ones(nE), np.ones(nE)]),
                                (np.concatenate([r_, r_]), np.concatenate([i, j]))), shape=(nE, n))
            Wd = sp.diags(w)
            keepv = np.ones(n, bool)
            keepv[pin] = False
            Dk, Dp = D_[:, keepv].tocsc(), D_[:, ~keepv].tocsc()
            Ksolve = spla.factorized((Dk.T @ Wd @ Dk).tocsc())
            Xn = X.copy()
            for c in range(3):
                rhs = Dk.T @ (w * (d[:, c] - Dp @ X[~keepv, c]))
                Xn[keepv, c] = Ksolve(rhs)
            X = Xn
            facts, backs = 2, 6
        return SyncResult(X, Rp[:, :, 2], _hom(Rp, X), res, facts, backs, time.time() - t0)
