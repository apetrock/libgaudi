"""Edge spheres: the Lie sphere surface solver (MATH.md section 19.7).

A face holds a plane; an edge holds one sphere. Two per edge, read two ways:

- `hinge` (primal, across the edge, from the two face planes): the sphere
  tangent to both planes at half the face heights,
  `k_perp = tan(theta / 2) / d`, `d = (A_f + A_g) / (2 |e|)` (exact for the
  incircle of a polygon). A member of the 2-parameter family
  `span{T_f, T_g}^perp`, picked by the metric.
- `contact` (dual, along the edge, from the endpoint contact elements): the
  sphere common to `L_i` and `L_j`, which exists when the normal lines meet
  (the principal contact element net condition) and then satisfies the
  discrete Rodrigues relation `x_j - x_i = (n_i - n_j) / k`. Read in closed
  form `k_par = -(dn . dx) / |dx|^2` (equal to the common sphere's curvature
  on principal edges, exact on circles); `common_sphere` gives the Lie
  least-squares version and its residual (the twist the hinge cannot see).

Both are Euclidean invariants, so targets read off a target mesh are pose
free and the merit `sum (k(X) - k*)^2` is a fixed energy. Projections: the
hinge turns its wings about the edge to the target dihedral; the contact
edge takes the Rodrigues length for the live normal turning (the stretch
outlet a closed tube needs). Both gather through one constant matrix.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_ops as ops
import dupin_poisson as dp
import lie_ops as lo

# ------------------------------------------------------------- topology


@dataclass
class Hinges:
    """Interior edges as diamonds: `a -> b` runs along face `f`'s winding,
    `c` is `f`'s third vertex, `d` the other face `g`'s."""
    e: np.ndarray
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    d: np.ndarray
    f: np.ndarray
    g: np.ndarray


def hinges(m: ops.Mesh) -> Hinges:
    rows = []
    for e, fs in enumerate(m.edge_faces):
        if len(fs) != 2:
            continue
        f, g = fs
        i = int(np.nonzero(m.fe[f] == e)[0][0])
        a, b, c = m.F[f, i], m.F[f, (i + 1) % 3], m.F[f, (i + 2) % 3]
        d = [v for v in m.F[g] if v != a and v != b][0]
        rows.append((e, a, b, c, d, f, g))
    r = np.array(rows, dtype=int).reshape(-1, 7)
    return Hinges(*[r[:, k] for k in range(7)])


# ------------------------------------------------------------- read-outs


def face_geometry(m: ops.Mesh, X: np.ndarray):
    N, A = ops.face_normals(m, X)
    return N, A


def hinge_angles(h: Hinges, X: np.ndarray, Nf: np.ndarray) -> np.ndarray:
    """Signed dihedral: positive where the surface bends toward its normal
    (the curvature sign convention of `lie_ops`)."""
    t = X[h.b] - X[h.a]
    t /= np.linalg.norm(t, axis=1)[:, None]
    nf, ng = Nf[h.f], Nf[h.g]
    return -np.arctan2(np.einsum("ij,ij->i", t, np.cross(nf, ng)), np.einsum("ij,ij->i", nf, ng))


def _cot(p, q, o):
    u, v = p - o, q - o
    return np.einsum("ij,ij->i", u, v) / np.linalg.norm(np.cross(u, v), axis=1)


def hinge_dual(h: Hinges, X: np.ndarray) -> np.ndarray:
    """Circumcentric dual length `l* = (|e| / 2)(cot alpha + cot beta)`: the
    distance between the two faces' circumcentres across the edge. A sphere
    through a face's vertices touches the face plane's parallel at its
    circumcentre, so the across-edge sphere is tangent to both planes at
    `l* / 2` from the edge (exact for a polygon's incircle). It vanishes on
    the diagonal of a cyclic quad (four concyclic points are coplanar): such
    a hinge measures nothing across itself."""
    le = np.linalg.norm(X[h.b] - X[h.a], axis=1)
    return 0.5 * le * (_cot(X[h.a], X[h.b], X[h.c]) + _cot(X[h.a], X[h.b], X[h.d]))


def hinge_curvature(h: Hinges, X: np.ndarray, Nf=None, Af=None, m=None, floor=1e-12):
    """`k_perp = 2 tan(theta / 2) / l*` and the hinge's reliability weight
    `max(l*, 0) |e|` (its dual area); hinges with `l*` below `floor` read 0."""
    if Nf is None:
        Nf, Af = ops.face_normals(m, X)
    th = hinge_angles(h, X, Nf)
    ls = hinge_dual(h, X)
    k = np.where(ls > floor, 2 * np.tan(0.5 * th) / np.maximum(ls, floor), 0.0)
    return k, th


def contact_curvature(E: np.ndarray, X: np.ndarray, N: np.ndarray) -> np.ndarray:
    """`k_par = -(dn . dx) / |dx|^2` per edge, normals `N` at the vertices."""
    dx = X[E[:, 1]] - X[E[:, 0]]
    dn = N[E[:, 1]] - N[E[:, 0]]
    return -np.einsum("ij,ij->i", dn, dx) / np.einsum("ij,ij->i", dx, dx)


def common_sphere(p_i, n_i, p_j, n_j):
    """The sphere in `L_i` closest to also lying in `L_j`: `s = a P_i + b T_i`
    with rows `<s, w> = 0` for `w` spanning `L_j^perp` (4 x 2, smallest
    singular vector). Returns (k, relative residual); the residual is 0
    exactly when the normal lines meet."""
    Pi, Ti = lo.contact(p_i, n_i)
    _, perp = lo.contact_frame(p_j, n_j)
    Pn, Tn = Pi / np.linalg.norm(Pi), Ti / np.linalg.norm(Ti)
    M = np.stack([perp.T @ (lo.ETA @ Pn), perp.T @ (lo.ETA @ Tn)], 1)
    _, s, Vt = np.linalg.svd(M)
    a, b = Vt[1]
    a, b = a / np.linalg.norm(Pi), b / np.linalg.norm(Ti)
    k = a / b if abs(b) > 1e-300 else np.inf
    return k, s[1] / s[0]


def hinge_weight(h: Hinges, X: np.ndarray) -> np.ndarray:
    le = np.linalg.norm(X[h.b] - X[h.a], axis=1)
    return np.maximum(hinge_dual(h, X), 0.0) * le


@dataclass
class EdgeTargets:
    k_perp: np.ndarray  # per hinge
    k_par: np.ndarray  # per edge
    w_perp: np.ndarray | None = None  # per hinge reliability (dual area), mean 1


def read_targets(m: ops.Mesh, h: Hinges, X: np.ndarray) -> EdgeTargets:
    Nf, Af = ops.face_normals(m, X)
    kp, _ = hinge_curvature(h, X, Nf, Af)
    w = hinge_weight(h, X)
    return EdgeTargets(kp, contact_curvature(m.E, X, dp.vertex_normals(m, X)), w / max(w.mean(), 1e-300))


# ---------------------------------------------------------------- solver


@dataclass
class EdgeConfig:
    iterations: int = 200
    w_hinge: float = 1.0
    w_contact: float = 1.0
    eps: float = 1e-4  # proximal anchor (translation null space), relative
    max_stretch: float = 0.25  # cap on the length change against the memory
    kmin_rel: float = 1e-3  # |k*| below this (relative to the mean |k*|): keep the length
    # Length memory (a soft material term): each edge's target length blends
    # the Rodrigues length (trusted as (|k*| / mean |k*|)^2, capped at 1: flat
    # directions say nothing about length) with its rest length. The rest is
    # the start (`restart` 0, elastic) or reseeded every `restart` iterates
    # (plastic, so the metric can grow).
    w_memory: float = 1.0
    restart: int = 20
    tension: float = 0.0
    tol: float = 1e-12
    patience: int = 40


@dataclass
class EdgeResult:
    X: np.ndarray
    history: list = field(default_factory=list)  # (move, merit)
    seconds: float = 0.0
    best_merit: float = np.inf


class EdgeSphereSolver:
    """Positions from edge-sphere targets by projective dynamics: each
    iterate projects every hinge (wings turned about the edge to the target
    dihedral) and every contact edge (Rodrigues length for the live normal
    turning), then one constant factor gathers them."""

    def __init__(self, m: ops.Mesh, X0: np.ndarray, tgt: EdgeTargets, cfg: EdgeConfig = EdgeConfig()):
        self.m, self.cfg, self.tgt = m, cfg, tgt
        self.h = hinges(m)
        nE, nH, n = len(m.E), len(self.h.e), m.n
        rows, cols, vals = [], [], []
        r = 0
        # Contact + memory rows: sqrt(w_c + w_m) (x_j - x_i), the target length
        # blended per iterate.
        wc = np.sqrt(cfg.w_contact + cfg.w_memory)
        self.sw_contact = wc
        for e, (i, j) in enumerate(m.E):
            rows += [r, r]
            cols += [i, j]
            vals += [-wc, wc]
            r += 1
        self.r_hinge = r
        # Hinge rows: sqrt(w_h w_k) (I - 11^T / 4) over the diamond (a, b, c, d),
        # `w_k` the target's reliability (dual area; 0 on cyclic-quad diagonals).
        wk = tgt.w_perp if tgt.w_perp is not None else np.ones(nH)
        self.sw_hinge = np.sqrt(cfg.w_hinge * wk)
        for k in range(nH):
            ids = [self.h.a[k], self.h.b[k], self.h.c[k], self.h.d[k]]
            for p in range(4):
                for q in range(4):
                    rows.append(r + p)
                    cols.append(ids[q])
                    vals.append(self.sw_hinge[k] * ((1.0 if p == q else 0.0) - 0.25))
            r += 4
        self.A = sp.csr_matrix((vals, (rows, cols)), shape=(r, n))
        self.Mv = ops.vertex_areas(m, X0)
        K = self.A.T @ self.A
        self.eps = cfg.eps * K.diagonal().mean() / self.Mv.mean()
        K = K + self.eps * sp.diags(self.Mv)
        if cfg.tension > 0:
            K = K + cfg.tension * dp.positive_laplacian(m, X0)
        self.solve = spla.factorized(K.tocsc())
        kscale = np.mean(np.abs(np.concatenate([tgt.k_perp, tgt.k_par]))) + 1e-300
        self.kmin = cfg.kmin_rel * kscale
        self.kbar = np.mean(np.abs(tgt.k_par)) + 1e-300
        self.L_rest = np.linalg.norm(X0[m.E[:, 1]] - X0[m.E[:, 0]], axis=1)

    def merit(self, X: np.ndarray) -> float:
        cur = read_targets(self.m, self.h, X)
        wk = self.sw_hinge ** 2 / self.cfg.w_hinge
        return float(np.sqrt(np.mean(wk * (cur.k_perp - self.tgt.k_perp) ** 2) +
                             np.mean((cur.k_par - self.tgt.k_par) ** 2)))

    def projections(self, X: np.ndarray, k_par: np.ndarray | None = None, k_perp: np.ndarray | None = None):
        """Contact spokes (per edge) and centred hinge diamonds (per hinge),
        toward `k_par` / `k_perp` (default: the targets)."""
        m, h, cfg = self.m, self.h, self.cfg
        k_par = self.tgt.k_par if k_par is None else k_par
        k_perp = self.tgt.k_perp if k_perp is None else k_perp
        N = dp.vertex_normals(m, X)
        # Contact: Rodrigues length along the live edge direction.
        dx = X[m.E[:, 1]] - X[m.E[:, 0]]
        le = np.linalg.norm(dx, axis=1)
        t = dx / le[:, None]
        turn = -np.einsum("ij,ij->i", N[m.E[:, 1]] - N[m.E[:, 0]], t)  # = k le
        k = k_par
        ok = (np.abs(k) > self.kmin) & (turn * k > 0)
        L_rod = self.L_rest.copy()
        L_rod[ok] = turn[ok] / k[ok]
        conf = np.where(ok, np.minimum((np.abs(k) / self.kbar) ** 2, 1.0), 0.0) * cfg.w_contact
        L = (conf * L_rod + cfg.w_memory * self.L_rest) / np.maximum(conf + cfg.w_memory, 1e-300)
        L = np.clip(L, self.L_rest * (1 - cfg.max_stretch), self.L_rest * (1 + cfg.max_stretch))
        d_contact = L[:, None] * t
        # Hinge: turn the wings about the edge to the target dihedral
        # (turning keeps both triangles, so `l*` is unchanged).
        Nf, Af = ops.face_normals(m, X)
        th = hinge_angles(h, X, Nf)
        ls = hinge_dual(h, X)
        th_t = np.where(ls > 0, 2 * np.arctan(0.5 * k_perp * np.maximum(ls, 0)), th)
        delta = th_t - th
        a, b, c, d = X[h.a], X[h.b], X[h.c], X[h.d]
        axis = (b - a) / np.linalg.norm(b - a, axis=1)[:, None]
        c2 = rotate_about(c, a, axis, 0.5 * delta)
        d2 = rotate_about(d, a, axis, -0.5 * delta)
        P = np.stack([a, b, c2, d2], 1)
        P -= P.mean(1, keepdims=True)
        return d_contact, P

    def step(self, X: np.ndarray) -> np.ndarray:
        d_contact, P = self.projections(X)
        rhs_rows = np.concatenate([self.sw_contact * d_contact,
                                   (self.sw_hinge[:, None, None] * P).reshape(-1, 3)])
        rhs = self.A.T @ rhs_rows + self.eps * self.Mv[:, None] * X
        return np.stack([self.solve(rhs[:, c]) for c in range(3)], 1)

    def run(self, X: np.ndarray, iterations: int | None = None, log=None, keep_best: bool = True) -> EdgeResult:
        t0 = time.time()
        out = EdgeResult(X.copy())
        best, best_m, since = X.copy(), self.merit(X), 0
        out.history.append((0.0, best_m))
        E = self.m.E
        self.L_rest = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1)
        for k in range(iterations or self.cfg.iterations):
            if self.cfg.restart > 0 and k > 0 and k % self.cfg.restart == 0:
                self.L_rest = np.linalg.norm(out.X[E[:, 1]] - out.X[E[:, 0]], axis=1)
            Xn = self.step(out.X)
            move = np.abs(Xn - out.X).max()
            out.X = Xn
            mer = self.merit(Xn)
            out.history.append((move, mer))
            if mer < best_m:
                best, best_m, since = Xn.copy(), mer, 0
            else:
                since += 1
            if log:
                log(k, Xn, move, mer)
            if move < self.cfg.tol or not np.isfinite(mer) or (self.cfg.patience and since >= self.cfg.patience):
                break
        out.best_merit = best_m
        if keep_best:
            out.X = best
        out.seconds = time.time() - t0
        return out


def rotate_about(p, origin, axis, angle):
    """Rodrigues rotation of points `p` about the line (origin, unit axis)."""
    v = p - origin
    c, s = np.cos(angle)[:, None], np.sin(angle)[:, None]
    kxv = np.cross(axis, v)
    kdv = np.einsum("ij,ij->i", axis, v)[:, None]
    return origin + v * c + kxv * s + axis * kdv * (1 - c)
