"""Vertex cyclides + edge spheres, coupled (MATH.md 19.9): the original
sketch of the finite-element cyclidic net.

Shape step, each iterate:
1. Pose the vertex cyclides' pose-free targets on the live mesh (principal
   frames) and couple them by the shape Poisson (`dupin_poisson`).
2. Each vertex predicts its edges' spheres by Euler, linear on the contact
   line: land the vertex on its cyclide, read `k_u, k_v` and the principal
   direction `e_U = dp / dth_v` (the U sphere held fixed traces the U
   curvature line); for an edge at angle `phi` to `e_U`, along
   `k_u cos^2 + k_v sin^2`, across `k_u sin^2 + k_v cos^2`.
3. Glue: each edge's along / across sphere is the blend of its own target and
   its endpoints' predictions, weighted by the cyclide's anisotropy
   `a = |k_u - k_v| / (|k_u| + |k_v|)` (0 at umbilics, where the cyclide
   degenerates and the edge spheres alone drive).

Position step: one constant factor gathering three projective-dynamics
terms: ring landing on each cyclide (blended toward the live spoke by
`1 - a`), hinge wing turns to the glued across sphere, and the glued along
sphere's Rodrigues length blended with the plastic length memory
(`lie_edges`).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_ops as ops
import dupin_poisson as dp
import lie_edges as le
import lie_ops as lo


@dataclass
class NetConfig:
    iterations: int = 300
    w_vertex: float = 1.0  # ring landing rows
    w_target: float = 1.0  # edge spheres: pull of the edge's own target
    w_couple: float = 1.0  # edge spheres: pull of each endpoint cyclide (times a_i)
    w_shape: float = 10.0  # shape Poisson on the posed cyclides
    land_iterations: int = 4
    edge: le.EdgeConfig = field(default_factory=le.EdgeConfig)
    patience: int = 60
    tol: float = 1e-12


@dataclass
class VertexReading:
    ok: np.ndarray  # cyclide usable at this vertex
    p: np.ndarray  # landed contact element
    n: np.ndarray
    ku: np.ndarray
    kv: np.ndarray
    eU: np.ndarray  # unit principal direction of the U family
    aniso: np.ndarray
    conic: list


def principal_direction(EU, EV, th, h=1e-6):
    """`e_U = dp / dth_v`: with the U sphere fixed the point runs along the U
    family's curvature line."""
    pp = lo.dupin_contact(EU, EV, (th[0], th[1] + h))
    pm = lo.dupin_contact(EU, EV, (th[0], th[1] - h))
    if pp is None or pm is None:
        return None
    d = pp[0] - pm[0]
    nd = np.linalg.norm(d)
    return d / nd if nd > 0 else None


def read_vertices(sig: np.ndarray, X: np.ndarray, N: np.ndarray, its: int) -> VertexReading:
    n = len(X)
    out = VertexReading(np.zeros(n, bool), X.copy(), N.copy(), np.zeros(n), np.zeros(n), np.zeros((n, 3)),
                        np.zeros(n), [None] * n)
    for i in range(n):
        try:
            EU, EV = lo.dupin_conic_basis(sig[i])
        except (ValueError, np.linalg.LinAlgError):
            continue
        p, nn, th, ok = lo.dupin_land(sig[i], X[i], n0=N[i], conic=(EU, EV), iterations=its)
        if not ok:
            continue
        e = principal_direction(EU, EV, th)
        if e is None:
            continue
        ku, kv, _, _ = lo.curvatures(sig[i], p, nn)
        if not (np.isfinite(ku) and np.isfinite(kv)):
            continue
        out.ok[i] = True
        out.p[i], out.n[i], out.ku[i], out.kv[i], out.eU[i] = p, nn, ku, kv, e
        out.aniso[i] = abs(ku - kv) / max(abs(ku) + abs(kv), 1e-300)
        out.conic[i] = (EU, EV)
    return out


def euler_predictions(E: np.ndarray, X: np.ndarray, vr: VertexReading, from_first: bool):
    """Along and across normal curvature of each edge `(i, j)` from vertex `i`
    (`from_first`) or `j`; weight `a` (0 where the cyclide is unusable)."""
    i, j = (E[:, 0], E[:, 1]) if from_first else (E[:, 1], E[:, 0])
    t = X[j] - X[i]
    n = vr.n[i]
    t = t - np.einsum("ij,ij->i", t, n)[:, None] * n
    t /= np.maximum(np.linalg.norm(t, axis=1), 1e-300)[:, None]
    c2 = np.einsum("ij,ij->i", t, vr.eU[i]) ** 2
    s2 = 1 - c2
    along = vr.ku[i] * c2 + vr.kv[i] * s2
    across = vr.ku[i] * s2 + vr.kv[i] * c2
    w = np.where(vr.ok[i], vr.aniso[i], 0.0)
    return along, across, w


class NetSolver:
    def __init__(self, m: ops.Mesh, X0: np.ndarray, sigma_body: np.ndarray, tgt: le.EdgeTargets,
                 cfg: NetConfig = NetConfig()):
        self.m, self.cfg, self.sigma_body, self.tgt = m, cfg, sigma_body, tgt
        self.edges = le.EdgeSphereSolver(m, X0, tgt, cfg.edge)
        self.h = self.edges.h
        self.ring = dp.rings(m)
        self.ring2 = dp.two_rings(self.ring)
        self.pairs = np.array([(i, j) for i in range(m.n) for j in self.ring[i]])
        npair = len(self.pairs)
        rows = np.arange(npair)
        sv = np.sqrt(cfg.w_vertex)
        Av = sp.csr_matrix((np.concatenate([-sv * np.ones(npair), sv * np.ones(npair)]),
                            (np.concatenate([rows, rows]), np.concatenate([self.pairs[:, 0], self.pairs[:, 1]]))),
                           shape=(npair, m.n))
        self.Av = Av
        Ae = self.edges.A
        K = Ae.T @ Ae + Av.T @ Av
        self.Mv = ops.vertex_areas(m, X0)
        self.eps = cfg.edge.eps * K.diagonal().mean() / self.Mv.mean()
        K = K + self.eps * sp.diags(self.Mv)
        self.solve = spla.factorized(K.tocsc())
        self.shape = (dp.ShapePoisson(dp.positive_laplacian(m, X0), self.Mv, cfg.w_shape)
                      if cfg.w_shape > 0 else None)
        # Map hinges to edges for the across predictions.
        self.hinge_edge = self.h.e
        self.last = {}

    def glue(self, X, vr):
        E = self.m.E
        a1, c1, w1 = euler_predictions(E, X, vr, True)
        a2, c2, w2 = euler_predictions(E, X, vr, False)
        wt, wc = self.cfg.w_target, self.cfg.w_couple
        k_par = (wt * self.tgt.k_par + wc * (w1 * a1 + w2 * a2)) / (wt + wc * (w1 + w2))
        he = self.hinge_edge
        k_perp = (wt * self.tgt.k_perp + wc * (w1[he] * c1[he] + w2[he] * c2[he])) / (wt + wc * (w1[he] + w2[he]))
        return k_par, k_perp

    def ring_spokes(self, X, sig, vr):
        """Ring landing on each cyclide, blended toward the live spoke by
        `1 - a_i`."""
        d = X[self.pairs[:, 1]] - X[self.pairs[:, 0]]
        its = self.cfg.land_iterations
        N = dp.vertex_normals(self.m, X)
        for k, (i, j) in enumerate(self.pairs):
            if not vr.ok[i] or vr.aniso[i] <= 0:
                continue
            p, _, _, ok = lo.dupin_land(sig[i], X[j], n0=N[j], conic=vr.conic[i], iterations=its)
            if ok:
                a = vr.aniso[i]
                d[k] = (1 - a) * d[k] + a * (p - vr.p[i])
        return d

    def step(self, X):
        cfg = self.cfg
        N = dp.vertex_normals(self.m, X)
        B = dp.body_frames(self.m, X, N, self.ring, "principal", self.ring2)
        sig = dp.to_world(B, self.sigma_body)
        sig = np.stack([dp.canonical_sign(s, X[i], N[i]) for i, s in enumerate(sig)])
        if self.shape is not None:
            sig = self.shape(sig)
        vr = read_vertices(sig, X, N, cfg.land_iterations)
        k_par, k_perp = self.glue(X, vr)
        d_contact, P = self.edges.projections(X, k_par, k_perp)
        d_ring = self.ring_spokes(X, sig, vr)
        rows = np.concatenate([self.edges.sw_contact * d_contact,
                               (self.edges.sw_hinge[:, None, None] * P).reshape(-1, 3)])
        rhs = self.edges.A.T @ rows + self.Av.T @ (np.sqrt(cfg.w_vertex) * d_ring) + self.eps * self.Mv[:, None] * X
        self.last = {"usable": vr.ok.mean(), "aniso": vr.aniso[vr.ok].mean() if vr.ok.any() else 0.0}
        return np.stack([self.solve(rhs[:, c]) for c in range(3)], 1)

    def run(self, X, iterations=None, log=None, keep_best=True) -> le.EdgeResult:
        t0 = time.time()
        out = le.EdgeResult(X.copy())
        E = self.m.E
        self.edges.L_rest = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1)
        best, best_m, since = X.copy(), self.edges.merit(X), 0
        out.history.append((0.0, best_m))
        restart = self.cfg.edge.restart
        for k in range(iterations or self.cfg.iterations):
            if restart > 0 and k > 0 and k % restart == 0:
                self.edges.L_rest = np.linalg.norm(out.X[E[:, 1]] - out.X[E[:, 0]], axis=1)
            Xn = self.step(out.X)
            move = np.abs(Xn - out.X).max()
            out.X = Xn
            mer = self.edges.merit(Xn)
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
