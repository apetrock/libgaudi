"""Lie shape projection: edge targets projected onto what Lie sphere elements
express, then positions (MATH.md 19.16).

Targets are the edge shape functions `H_e = (a, tau, b)` (along, twist,
across; `lie_edge_net.read_spheres`). Each frame:

1. Vertex curvature spheres. Each vertex's ring-fitted tensor (its Euler
   fit to the incident edges) gives its two curvature spheres `s_u` (larger
   |k|) and `s_v`. (Per-edge principal decompositions were tried: each edge
   mixes the contact and hinge read-outs and swaps families on diagonals.)
2. Implicit vertex cyclides. A Dupin cyclide's U-family curvature spheres
   all lie in one 3-dim subspace `U` (signature (2, 1)), the other family in
   `U^perp`. Per vertex, the one ring's sphere moments `M_u = sum w s_u s_u^T`
   and `M_v = sum w (eta s_v)(eta s_v)^T` (weights: anisotropy; the vertex's
   own spheres `w_self` times a neighbour's), smoothed across vertices by the
   Poisson `(w_c L + M)`, then `U_i` = the top three eigenvectors of
   `M_u - M_v` in the vertex's local frame (translated to the vertex and
   scaled by the mean edge length), and `sigma_i = from_subspace(U_i)`. No
   stored state, no fit. Exact on exact data (5e-14 on the analytic torus);
   on read-out data the third-order part is poorly conditioned (one ring's
   spheres are nearly collinear) but the read-back at the vertex is faithful.
3. Central spheres `Y_i = H_i P(x_i) + T(x_i, n_i)` (`H_i` the ring mean of
   `(a + b) / 2`; on any contact line `s_a + s_b = s_u + s_v = 2 Y`),
   smoothed by the same Poisson (the conformal Gauss map: its Dirichlet
   energy is the Willmore energy, Bryant 1984).
4. Read back: `H*` from `Y` projected onto each contact line; the cyclide's
   `k_u - k_v` and principal direction where it exists (else the target's
   own deviatoric part); per-edge `(a, tau, b)` by Euler from the endpoint
   tensors.

Positions from the projected targets (`positions_*`): the one-ring spoke
gather (spokes set on their along spheres, one prefactored Laplacian),
envelope landing on the central spheres (local), the mean-curvature Poisson
(`Delta x = 2 H n`), or `contact_net.ContactNetGN` as the reference.
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
import lie_edge_net as lsn
import lie_net as ln
import lie_ops as lo

ETA = lo.ETA


def _nrm(v):
    return v / np.maximum(np.linalg.norm(v, axis=1), 1e-300)[:, None]


# -------------------------------------------------------------- spheres


def edge_frames(m: ops.Mesh, X, N):
    i, j = m.E[:, 0], m.E[:, 1]
    mid = 0.5 * (X[i] + X[j])
    nb = _nrm(N[i] + N[j])
    d = X[j] - X[i]
    t = _nrm(d - np.einsum("ij,ij->i", d, nb)[:, None] * nb)
    return mid, nb, t, np.cross(nb, t)


def edge_principal(H):
    """Per edge `(k_u, k_v, anisotropy)` of `[[a, tau], [tau, b]]`."""
    a, tau, b = H[:, 0], H[:, 1], H[:, 2]
    mean, dev = 0.5 * (a + b), np.sqrt(0.25 * (a - b) ** 2 + tau ** 2)
    k1, k2 = mean + dev, mean - dev
    big = np.abs(k1) >= np.abs(k2)
    ku, kv = np.where(big, k1, k2), np.where(big, k2, k1)
    return ku, kv, np.abs(ku - kv) / np.maximum(np.abs(ku) + np.abs(kv), 1e-300)


def central_curvature(m: ops.Mesh, X, H):
    """Ring mean of `(a + b) / 2`, weighted by edge length."""
    le = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1)
    s = np.zeros(m.n)
    w = np.zeros(m.n)
    for k in range(2):
        np.add.at(s, m.E[:, k], le * 0.5 * (H[:, 0] + H[:, 2]))
        np.add.at(w, m.E[:, k], le)
    return s / np.maximum(w, 1e-300)


def local_frame(x, h):
    """`A = scale(1 / h) translation(-x)`: the vertex at the origin, unit
    edge length (conditioning of the eigen step; rigid + scale equivariant)."""
    return lo.scale(1.0 / h) @ lo.translation(-x)


# ---------------------------------------------------------- projection


@dataclass
class ProjectionConfig:
    w_couple: float = 1.0  # Poisson coupling strength (dimensionless; 0: raw targets)
    min_aniso: float = 0.05  # edges / vertices below: no cyclide information
    beta_v: float = 1.0  # weight of the U^perp family's rows in the eigen step
    w_self: float = 10.0  # the vertex's own curvature spheres against each ring neighbour's
    use_cyclide: bool = True
    land_iterations: int = 6


@dataclass
class Projection:
    H: np.ndarray  # (nE, 3) projected (a, tau, b)
    A: np.ndarray  # (n, 3, 3) projected ambient vertex tensors
    Hmean: np.ndarray  # (n,) filtered mean curvature
    ok: np.ndarray  # (n,) vertex has a cyclide
    sigma: np.ndarray  # (n, 6, 6)
    dupin_defect: np.ndarray  # (n,) leftover eigenvalue share
    seconds: float = 0.0


class LieShapeProjection:
    def __init__(self, m: ops.Mesh, cfg: ProjectionConfig = ProjectionConfig()):
        self.m, self.cfg = m, cfg
        self.ring = dp.rings(m)
        self.incident = [[] for _ in range(m.n)]
        for e, (i, j) in enumerate(m.E):
            self.incident[i].append(e)
            self.incident[j].append(e)

    def _poisson(self, X):
        L = dp.positive_laplacian(self.m, X)
        Mv = ops.vertex_areas(self.m, X)
        scale = Mv.mean() / max(L.diagonal().mean(), 1e-300)
        return spla.factorized((self.cfg.w_couple * scale * L + sp.diags(Mv)).tocsc()), Mv

    def __call__(self, X, N, H) -> Projection:
        t0 = time.time()
        m, cfg = self.m, self.cfg
        n, nE = m.n, len(m.E)
        N = _nrm(N)
        h = float(np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1).mean())
        Hm = central_curvature(m, X, H)
        if cfg.w_couple <= 0:
            A = self.ambient_from_edges(X, N, H)
            return Projection(H.copy(), A, Hm, np.zeros(n, bool), np.zeros((n, 6, 6)), np.zeros(n),
                              time.time() - t0)
        solve, Mv = self._poisson(X)
        # Central spheres, smoothed, read back on each contact line.
        Y = lsn.curvature_sphere_rows(X, N, Hm)
        Ys = np.column_stack([solve(Mv * Y[:, c]) for c in range(6)])
        Hs = np.zeros(n)
        for i in range(n):
            T = lo.translation(-X[i])
            y = T @ Ys[i]
            P0 = np.array([0.5, 0.5, 0, 0, 0, 0.0])
            T0 = np.concatenate([[0.0, 0.0], N[i], [1.0]])
            ab = np.linalg.lstsq(np.stack([P0, T0], 1), y, rcond=None)[0]
            Hs[i] = ab[0] / ab[1] if abs(ab[1]) > 1e-12 else Hm[i]
        # Implicit cyclides from the ring's curvature spheres.
        ok = np.zeros(n, bool)
        sig = np.zeros((n, 6, 6))
        defect = np.zeros(n)
        A = self.ambient_from_edges(X, N, H)
        if cfg.use_cyclide:
            # Curvature spheres at the vertices from each vertex's ring-fitted
            # tensor (per-edge principal decompositions mix the contact and
            # hinge read-outs and swap families on diagonals).
            kuv = np.zeros(n)
            kvv = np.zeros(n)
            for i in range(n):
                ev = np.linalg.eigvalsh(0.5 * (A[i] + A[i].T))
                ev = ev[np.argsort(-np.abs(ev))]
                kuv[i], kvv[i] = ev[0], ev[1]
            anv = np.abs(kuv - kvv) / np.maximum(np.abs(kuv) + np.abs(kvv), 1e-300)
            su = lsn.curvature_sphere_rows(X, N, kuv)
            esv = lsn.curvature_sphere_rows(X, N, kvv) @ ETA
            w = np.where(anv > cfg.min_aniso, anv, 0.0)
            Ou = w[:, None, None] * np.einsum("ni,nj->nij", su, su)
            Ov = w[:, None, None] * np.einsum("ni,nj->nij", esv, esv)
            Mu = cfg.w_self * Ou
            Mvv = cfg.w_self * Ov
            for a_, b_ in m.E:
                Mu[a_] += Ou[b_]
                Mu[b_] += Ou[a_]
                Mvv[a_] += Ov[b_]
                Mvv[b_] += Ov[a_]
            Mu = np.stack([solve(Mv * Mu.reshape(n, 36)[:, c]) for c in range(36)], 1).reshape(n, 6, 6)
            Mvv = np.stack([solve(Mv * Mvv.reshape(n, 36)[:, c]) for c in range(36)], 1).reshape(n, 6, 6)
            wsum = cfg.w_self * w.copy()
            for a_, b_ in m.E:
                wsum[a_] += w[b_]
                wsum[b_] += w[a_]
            for i in range(n):
                if wsum[i] <= 0:
                    continue
                Al = local_frame(X[i], h)
                Ali = lo.inverse(Al)
                Ml = Al @ Mu[i] @ Al.T - cfg.beta_v * Ali.T @ Mvv[i] @ Ali
                Ml = 0.5 * (Ml + Ml.T)
                ev, V = np.linalg.eigh(Ml)
                order = np.argsort(-ev)
                B = Ali @ V[:, order[:3]]
                G = B.T @ ETA @ B
                gev = np.linalg.eigvalsh(0.5 * (G + G.T))
                if (gev < 0).sum() != 1 or (gev > 0).sum() != 2:
                    continue
                s = lo.from_subspace(B)
                if not np.all(np.isfinite(s)) or lo.defect(s) > 1e-6:
                    continue
                sig[i] = dp.canonical_sign(s, X[i], N[i])
                tot = np.abs(ev).sum()
                defect[i] = np.abs(ev[order[3:]]).sum() / max(tot, 1e-300)
                ok[i] = True
            if ok.any():
                vr = ln.read_vertices(sig, X, N, cfg.land_iterations)
                ok &= vr.ok
                for i in np.nonzero(ok)[0]:
                    if vr.aniso[i] < cfg.min_aniso:
                        ok[i] = False
                        continue
                    nn = N[i]
                    e1 = vr.eU[i] - (vr.eU[i] @ nn) * nn
                    e1 /= max(np.linalg.norm(e1), 1e-300)
                    e2 = np.cross(nn, e1)
                    dev = 0.5 * (vr.ku[i] - vr.kv[i])
                    A[i] = Hs[i] * (np.eye(3) - np.outer(nn, nn)) + dev * (np.outer(e1, e1) - np.outer(e2, e2))
        # Vertices without a cyclide: the filtered mean, the target's deviatoric.
        for i in np.nonzero(~ok)[0]:
            nn = N[i]
            Pt = np.eye(3) - np.outer(nn, nn)
            tr = np.trace(Pt @ A[i] @ Pt)
            A[i] = Pt @ A[i] @ Pt + (Hs[i] - 0.5 * tr) * Pt
        Hp = self.edges_from_ambient(X, N, A)
        return Projection(Hp, A, Hs, ok, sig, defect, time.time() - t0)

    def ambient_from_edges(self, X, N, H):
        """Per vertex the ambient tensor of its least-squares Euler fit to the
        incident edges (`lie_edge_net.fit_vertex_tensors`)."""
        R = lsn.reference_dirs(N)
        S = lsn.fit_vertex_tensors(self.m, X, N, H, R)
        T, B, _ = lsn.vertex_frames(N, R)
        F = np.stack([T, B], 2)
        S2 = np.stack([np.stack([S[:, 0], S[:, 1]], 1), np.stack([S[:, 1], S[:, 2]], 1)], 1)
        return F @ S2 @ np.transpose(F, (0, 2, 1))

    def edges_from_ambient(self, X, N, A):
        mid, nb, t, b = edge_frames(self.m, X, N)
        Ae = 0.5 * (A[self.m.E[:, 0]] + A[self.m.E[:, 1]])
        return np.column_stack([np.einsum("ei,eij,ej->e", t, Ae, t), np.einsum("ei,eij,ej->e", t, Ae, b),
                                np.einsum("ei,eij,ej->e", b, Ae, b)])


# ------------------------------------------------------------ positions


@dataclass
class PositionResult:
    X: np.ndarray
    iterations: int = 0
    factorizations: int = 0
    backsolves: int = 0
    seconds: float = 0.0
    history: list = field(default_factory=list)


class SpokeGather:
    """One-ring spokes set on their along spheres, one prefactored Laplacian.

    Local: from `x_i` along the spoke's tangent direction `t` (in `n_i`'s
    plane), the neighbour on the circle of curvature `kappa = t^T A_i t` at
    length `l`: `d = (sin(kappa l) / kappa) t + ((1 - cos(kappa l)) / kappa) n_i`.
    Global: `(A^T W A + eps M) x = A^T W d + eps M x_prev`. Lengths: the rest
    (reseeded every `restart` iterates) blended with the Rodrigues length
    `turning / |kappa|` by `rodrigues`."""

    def __init__(self, m: ops.Mesh, X0, eps=1e-3, restart=10, rodrigues=0.5, anderson=5):
        self.m, self.restart, self.rodrigues, self.anderson = m, restart, rodrigues, anderson
        ring = dp.rings(m)
        self.pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])
        wE = np.maximum(ops.cotan_weights(m, X0), 1e-3)
        index = {(min(a, b), max(a, b)): k for k, (a, b) in enumerate(m.E)}
        self.w = np.array([wE[index[(min(i, j), max(i, j))]] for i, j in self.pairs])
        P = len(self.pairs)
        r = np.arange(P)
        self.Aop = sp.csr_matrix((np.concatenate([-np.ones(P), np.ones(P)]),
                                  (np.concatenate([r, r]), np.concatenate([self.pairs[:, 0], self.pairs[:, 1]]))),
                                 shape=(P, m.n))
        K = self.Aop.T @ sp.diags(self.w) @ self.Aop
        self.Mv = ops.vertex_areas(m, X0)
        self.eps = eps * K.diagonal().mean() / self.Mv.mean()
        self.solve = spla.factorized((K + self.eps * sp.diags(self.Mv)).tocsc())

    def targets(self, X, A, L):
        i, j = self.pairs[:, 0], self.pairs[:, 1]
        N = dp.vertex_normals(self.m, X)
        d = X[j] - X[i]
        n = N[i]
        t = _nrm(d - np.einsum("ij,ij->i", d, n)[:, None] * n)
        kap = np.einsum("pi,pij,pj->p", t, A[i], t)
        turn = np.arccos(np.clip(np.einsum("ij,ij->i", N[i], N[j]), -1.0, 1.0))
        lr = np.where(np.abs(kap) > 1e-6, turn / np.maximum(np.abs(kap), 1e-300), L)
        lr = np.clip(lr, 0.5 * L, 2.0 * L)
        l = (1 - self.rodrigues) * L + self.rodrigues * lr
        x = kap * l
        small = np.abs(x) < 1e-6
        s1 = np.where(small, l * (1 - x * x / 6), np.sin(x) / np.where(small, 1.0, kap))
        s2 = np.where(small, 0.5 * kap * l * l, (1 - np.cos(x)) / np.where(small, 1.0, kap))
        return s1[:, None] * t + s2[:, None] * n

    def step(self, X, A, L, Xprev):
        d = self.targets(X, A, L)
        rhs = self.Aop.T @ (self.w[:, None] * d) + self.eps * self.Mv[:, None] * Xprev
        return np.column_stack([self.solve(rhs[:, c]) for c in range(3)])

    def residual(self, X, A, L):
        d = self.targets(X, A, L)
        r = self.Aop @ X - d
        return float((self.w[:, None] * r * r).sum())

    def run(self, X, A, iterations=60, tol=1e-10) -> PositionResult:
        t0 = time.time()
        out = PositionResult(X.copy(), factorizations=1)
        Xc = X.copy()
        L = np.linalg.norm(Xc[self.pairs[:, 1]] - Xc[self.pairs[:, 0]], axis=1)
        G_hist, F_hist = [], []
        for k in range(iterations):
            if self.restart and k > 0 and k % self.restart == 0:
                L = np.linalg.norm(Xc[self.pairs[:, 1]] - Xc[self.pairs[:, 0]], axis=1)
                G_hist, F_hist = [], []
            G = self.step(Xc, A, L, Xc)
            out.backsolves += 3
            Fk = (G - Xc).ravel()
            Xn = G
            if self.anderson:
                G_hist.append(G.ravel())
                F_hist.append(Fk)
                if len(F_hist) > self.anderson + 1:
                    G_hist.pop(0)
                    F_hist.pop(0)
                if len(F_hist) >= 2:
                    dF = np.stack([F_hist[q + 1] - F_hist[q] for q in range(len(F_hist) - 1)], 1)
                    dG = np.stack([G_hist[q + 1] - G_hist[q] for q in range(len(G_hist) - 1)], 1)
                    gam = np.linalg.lstsq(dF, Fk, rcond=None)[0]
                    cand = (G.ravel() - dG @ gam).reshape(-1, 3)
                    if np.all(np.isfinite(cand)) and self.residual(cand, A, L) <= self.residual(G, A, L):
                        Xn = cand
            move = float(np.abs(Xn - Xc).max())
            Xc = Xn
            out.history.append(move)
            if move < tol:
                break
        out.X, out.iterations, out.seconds = Xc, k + 1, time.time() - t0
        return out


def envelope_landing(m: ops.Mesh, X, Hs, iterations=60, w_self=10.0, eps=1e-2, gn=3, tol=1e-10,
                     relax=0.5) -> PositionResult:
    """Each vertex moves to the point on its central sphere (target mean
    curvature `Hs`, tangent at the current contact element) that best lies on
    its ring's central spheres: `min w_self f_i^2 + sum f_j^2 + eps |x - x_i|^2 / h^2`,
    `f_j(x) = n_j . (x - x_j) - H_j |x - x_j|^2 / 2` (the sphere's power, a
    signed distance near it). Local 3 x 3 Gauss-Newton per vertex, Jacobi
    updates, no global solve."""
    t0 = time.time()
    ring = dp.rings(m)
    deg = max(len(r) for r in ring)
    nbr = np.full((m.n, deg + 1), -1)
    for i, r in enumerate(ring):
        nbr[i, 0] = i
        nbr[i, 1:len(r) + 1] = r
    valid = nbr >= 0
    wts = np.where(valid, 1.0, 0.0)
    wts[:, 0] = w_self
    idx = np.where(valid, nbr, 0)
    h = float(np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1).mean())
    out = PositionResult(X.copy())
    Xc = X.copy()
    for k in range(iterations):
        N = dp.vertex_normals(m, Xc)
        P, Nn, Hh = Xc[idx], N[idx], Hs[idx]  # (n, deg+1, 3)
        x = Xc.copy()
        for _ in range(gn):
            dx = x[:, None, :] - P
            f = np.einsum("nki,nki->nk", Nn, dx) - 0.5 * Hh * np.einsum("nki,nki->nk", dx, dx)
            g = Nn - Hh[..., None] * dx
            JtJ = np.einsum("nk,nki,nkj->nij", wts, g, g) + (eps / h ** 2) * np.eye(3)[None]
            Jtr = np.einsum("nk,nki,nk->ni", wts, g, f) + (eps / h ** 2) * (x - Xc)
            x = x - np.linalg.solve(JtJ, Jtr[..., None])[..., 0]
        # Jacobi updates of neighbouring vertices overshoot together: relax,
        # and cap each move at a tenth of an edge.
        step = relax * (x - Xc)
        sn = np.linalg.norm(step, axis=1)
        step *= np.minimum(1.0, 0.1 * h / np.maximum(sn, 1e-300))[:, None]
        x = Xc + step
        move = float(np.abs(x - Xc).max())
        Xc = x
        out.history.append(move)
        if move < tol:
            break
    out.X, out.iterations, out.seconds = Xc, k + 1, time.time() - t0
    return out


def mean_curvature_poisson(m: ops.Mesh, X, Hs, iterations=60, mu=1.0, tol=1e-10) -> PositionResult:
    """`(mu M + L) x = mu M x_prev - 2 M H n` (`Delta x = 2 H n`, `L` the
    positive cotan Laplacian frozen at the start; `mu` relative)."""
    t0 = time.time()
    L = dp.positive_laplacian(m, X)
    Mv = ops.vertex_areas(m, X)
    mu_ = mu * L.diagonal().mean() / Mv.mean()
    solve = spla.factorized((mu_ * sp.diags(Mv) + L).tocsc())
    out = PositionResult(X.copy(), factorizations=1)
    Xc = X.copy()
    for k in range(iterations):
        N = dp.vertex_normals(m, Xc)
        rhs = mu_ * Mv[:, None] * Xc - 2 * Mv[:, None] * Hs[:, None] * N
        Xn = np.column_stack([solve(rhs[:, c]) for c in range(3)])
        out.backsolves += 3
        move = float(np.abs(Xn - Xc).max())
        Xc = Xn
        out.history.append(move)
        if move < tol:
            break
    out.X, out.iterations, out.seconds = Xc, k + 1, time.time() - t0
    return out


# --------------------------------------------------------------- solver


@dataclass
class SolverConfig:
    outer: int = 8  # projection sweeps (fields rebuilt on the current geometry)
    inner: int = 20  # position iterates per sweep
    route: str = "spoke"  # spoke | envelope | mcp | cnet
    projection: ProjectionConfig = field(default_factory=ProjectionConfig)


class LieShapeSolver:
    """Outer sweeps: project the targets on the current contact elements,
    then positions by the chosen route."""

    def __init__(self, m: ops.Mesh, X0, targets: lsn.EdgeSpheres, cfg: SolverConfig = SolverConfig()):
        self.m, self.cfg, self.tgt = m, cfg, targets
        self.proj = LieShapeProjection(m, cfg.projection)
        self.gather = SpokeGather(m, X0) if cfg.route == "spoke" else None
        ring = dp.rings(m)
        self.pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]])

    def solve(self, X, log=None):
        t0 = time.time()
        Xc = X.copy()
        hist = []
        facts = backs = 0
        tproj = tpos = 0.0
        for k in range(self.cfg.outer):
            N = dp.vertex_normals(self.m, Xc)
            pr = self.proj(Xc, N, self.tgt.H)
            tproj += pr.seconds
            route = self.cfg.route
            if route == "spoke":
                r = self.gather.run(Xc, pr.A, iterations=self.cfg.inner)
            elif route == "envelope":
                r = envelope_landing(self.m, Xc, pr.Hmean, iterations=self.cfg.inner)
            elif route == "mcp":
                r = mean_curvature_poisson(self.m, Xc, pr.Hmean, iterations=self.cfg.inner)
            else:
                d = Xc[self.pairs[:, 1]] - Xc[self.pairs[:, 0]]
                tt = _nrm(d - np.einsum("ij,ij->i", d, N[self.pairs[:, 0]])[:, None] * N[self.pairs[:, 0]])
                kap = np.einsum("pi,pij,pj->p", tt, pr.A[self.pairs[:, 0]], tt)
                c = 0.5 * kap * np.linalg.norm(d, axis=1)
                tg = cn.NetTargets(pr.H[:, 0], pr.H[:, 1], c)
                g = cn.ContactNetGN(self.m, Xc, tg, cn.NetConfig(), self.pairs)
                rg = g.run(Xc.copy(), N, iterations=self.cfg.inner)
                r = PositionResult(rg.X, len(rg.history) - 1, rg.factorizations, rg.backsolves, rg.seconds)
            tpos += r.seconds
            facts += r.factorizations
            backs += r.backsolves
            move = float(np.abs(r.X - Xc).max())
            Xc = r.X
            hist.append((move, int(pr.ok.sum())))
            if log:
                log(k, Xc, move, pr)
        return Xc, dict(history=hist, factorizations=facts, backsolves=backs, t_projection=tproj,
                        t_positions=tpos, seconds=time.time() - t0)
