"""Dupin Poisson solver prototype (MATH.md section 19).

Shape state: one Dupin cyclide `sigma_i` per vertex (`lie_ops`), fitted from
its one ring's contact elements and coupled across the ring by the ambient
(Frobenius) Dirichlet energy, which is the scalar Laplacian on each of the
36 entries: one constant factor `(w L + mu M)`, then each vertex projected
back onto the cyclides by the matrix sign function.

Positions (stage 2) chase points and normals only: each vertex's cyclide
predicts its ring (landing), gathered by a constant cotan solve. No CR
operators in the solver; `cr_ops` supplies meshes and the cotan matrix.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_ops as ops
import lie_ops as lo

# ------------------------------------------------------------ mesh helpers


def vertex_normals(m: ops.Mesh, X: np.ndarray) -> np.ndarray:
    a, b, c = X[m.F[:, 0]], X[m.F[:, 1]], X[m.F[:, 2]]
    cr = np.cross(b - a, c - a)
    N = np.zeros_like(X)
    for k in range(3):
        np.add.at(N, m.F[:, k], cr)
    return N / np.linalg.norm(N, axis=1)[:, None]


def rings(m: ops.Mesh) -> list[np.ndarray]:
    nb = [set() for _ in range(m.n)]
    for i, j in m.E:
        nb[i].add(j)
        nb[j].add(i)
    return [np.array(sorted(s)) for s in nb]


def positive_laplacian(m: ops.Mesh, X: np.ndarray, kind: str = "cotan") -> sp.csr_matrix:
    """PSD vertex Laplacian with nonnegative weights (obtuse cotans clipped):
    no negative couplings."""
    w = np.maximum(ops.cotan_weights(m, X), 0.0) if kind == "cotan" else np.ones(len(m.E))
    i, j = m.E[:, 0], m.E[:, 1]
    rows = np.concatenate([i, j, i, j])
    cols = np.concatenate([j, i, i, j])
    vals = np.concatenate([-w, -w, w, w])
    return sp.csr_matrix((vals, (rows, cols)), shape=(m.n, m.n))


def two_rings(ring: list[np.ndarray]) -> list[np.ndarray]:
    out = []
    for i, r in enumerate(ring):
        s = set(r.tolist()) | {i}
        for j in r:
            s |= set(ring[j].tolist())
        out.append(np.array(sorted(s)))
    return out


def tangent_frame(n):
    a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    t1 = np.cross(n, a)
    t1 /= np.linalg.norm(t1)
    return t1, np.cross(n, t1)


@dataclass
class Jet:
    """Height function `h(u, v)` over the tangent plane at `x0` (frame
    t1, t2, n), cubic: coefficients of u, v, u^2/2, uv, v^2/2, u^3, u^2 v,
    u v^2, v^3 (and a constant)."""
    x0: np.ndarray
    t1: np.ndarray
    t2: np.ndarray
    n: np.ndarray
    c: np.ndarray

    @staticmethod
    def basis(u, v):
        return np.stack([np.ones_like(u), u, v, 0.5 * u * u, u * v, 0.5 * v * v,
                         u ** 3, u * u * v, u * v * v, v ** 3], -1)

    @staticmethod
    def grad_basis(u, v):
        z = np.zeros_like(u)
        o = np.ones_like(u)
        du = np.stack([z, o, z, u, v, z, 3 * u * u, 2 * u * v, v * v, z], -1)
        dv = np.stack([z, z, o, z, u, v, z, u * u, 2 * u * v, 3 * v * v], -1)
        return du, dv

    def contact(self, u, v):
        """Point and unit normal of the jet surface at (u, v)."""
        h = self.basis(u, v) @ self.c
        du, dv = self.grad_basis(u, v)
        hu, hv = du @ self.c, dv @ self.c
        p = self.x0 + u[:, None] * self.t1 + v[:, None] * self.t2 + h[:, None] * self.n
        nn = -hu[:, None] * self.t1 - hv[:, None] * self.t2 + self.n
        return p, nn / np.linalg.norm(nn, axis=1)[:, None]

    def shape(self):
        """Normal, principal curvatures (k1 larger |k|) and directions at the
        jet's foot over x0 (u = v = 0)."""
        c = self.c
        g = np.array([c[1], c[2]])
        Hs = np.array([[c[3], c[4]], [c[4], c[5]]])
        # First and second fundamental forms of the graph at the origin.
        I = np.eye(2) + np.outer(g, g)
        w = np.sqrt(1 + g @ g)
        S = np.linalg.solve(I, Hs / w)
        ev, V = np.linalg.eig(S)
        ev, V = np.real(ev), np.real(V)
        o = np.argsort(-np.abs(ev))
        nrm = self.n - g[0] * self.t1 - g[1] * self.t2
        nrm /= np.linalg.norm(nrm)
        dirs = [V[0, k] * self.t1 + V[1, k] * self.t2 + (g @ V[:, k]) * self.n for k in o]
        dirs = [d / np.linalg.norm(d) for d in dirs]
        p = self.x0 + c[0] * self.n
        return p, nrm, ev[o], dirs


def fit_jet(X, ids, i, n0):
    t1, t2 = tangent_frame(n0)
    d = X[ids] - X[i]
    u, v, h = d @ t1, d @ t2, d @ n0
    c = np.linalg.lstsq(Jet.basis(u, v), h, rcond=None)[0]
    return Jet(X[i].copy(), t1, t2, n0.copy(), c)


def jet_contacts(m: ops.Mesh, X: np.ndarray, N: np.ndarray, ring2: list[np.ndarray]):
    """Per vertex: its cubic jet over the two ring, and the ring's points and
    normals read off the jet (consistent contact elements)."""
    jets = [fit_jet(X, ring2[i], i, N[i]) for i in range(m.n)]
    return jets


# ----------------------------------------------------------- shape field


def canonical_sign(sigma, p, n):
    """`sigma` or `-sigma` so that U holds the larger |curvature| at (p, n)."""
    ku, kv, _, _ = lo.curvatures(sigma, p, n)
    return sigma if abs(ku) >= abs(kv) else -sigma


@dataclass
class FieldFit:
    sigma: np.ndarray  # (n, 6, 6)
    gap: np.ndarray  # (n,) eigenvalue ratio of the fit, small when determined
    ok: np.ndarray  # (n,) bool


def fit_field(m: ops.Mesh, X: np.ndarray, N: np.ndarray, ring: list[np.ndarray] | None = None,
              two_ring: bool = False) -> FieldFit:
    ring = rings(m) if ring is None else ring
    sig = np.zeros((m.n, 6, 6))
    gap = np.zeros(m.n)
    ok = np.zeros(m.n, bool)
    for i in range(m.n):
        ids = set([i]) | set(ring[i].tolist())
        if two_ring:
            for j in ring[i]:
                ids |= set(ring[j].tolist())
        ids = sorted(ids)
        s, g, o = lo.fit(X[ids], N[ids], gap_tol=1e-2)
        sig[i] = canonical_sign(s, X[i], N[i])
        gap[i], ok[i] = g, o
    return FieldFit(sig, gap, ok)


class ShapePoisson:
    """`(w L + mu M) m = mu M m*` on each of the 36 entries, factored once."""

    def __init__(self, L: sp.csr_matrix, Mv: np.ndarray, w: float, mu: float = 1.0):
        self.Mv = Mv
        self.mu = mu
        self.solve = spla.factorized((w * L + mu * sp.diags(Mv)).tocsc())

    def __call__(self, sig: np.ndarray, project: bool = True) -> np.ndarray:
        n = len(sig)
        rhs = self.mu * self.Mv[:, None] * sig.reshape(n, 36)
        out = np.stack([self.solve(rhs[:, k]) for k in range(36)], 1).reshape(n, 6, 6)
        if project:
            out = np.stack([lo.normalize(s) for s in out])
        return out


def field_energy(m: ops.Mesh, L_w: np.ndarray, sig: np.ndarray) -> float:
    d = sig[m.E[:, 0]] - sig[m.E[:, 1]]
    return float((L_w * (d * d).sum((1, 2))).sum())


def conjugate_field(A: np.ndarray, sig: np.ndarray) -> np.ndarray:
    Ai = lo.inverse(A)
    return np.einsum("ab,nbc,cd->nad", A, sig, Ai)


def field_error(sig: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Per vertex max-abs distance, up to the family swap."""
    if truth.ndim == 2:
        truth = np.broadcast_to(truth, sig.shape)
    a = np.abs(sig - truth).max((1, 2))
    b = np.abs(sig + truth).max((1, 2))
    return np.minimum(a, b)


# ------------------------------------------------------------- read-out


def quadratic_jet_shape(X, ids, i, n0):
    """Principal curvatures (|k1| >= |k2|) and the k1 direction at x_i from a
    quadratic height fit over `ids` in the frame of `n0`."""
    t1, t2 = tangent_frame(n0)
    d = X[ids] - X[i]
    u, v, h = d @ t1, d @ t2, d @ n0
    B = np.stack([np.ones_like(u), u, v, 0.5 * u * u, u * v, 0.5 * v * v], 1)
    c = np.linalg.lstsq(B, h, rcond=None)[0]
    J = Jet(X[i].copy(), t1, t2, n0.copy(), np.concatenate([c, np.zeros(4)]))
    _, _, k, dirs = J.shape()
    return k[0], k[1], dirs[0]


READOUT = {"mode": "osculating"}


def readout(m: ops.Mesh, X: np.ndarray, N: np.ndarray | None = None, ring2=None, prior=1e-6,
            iterations=6, mode: str | None = None) -> np.ndarray:
    """Per-vertex cyclide of a mesh. `osculating`: the osculating torus of
    the quadratic jet (second order; the third-order freedom fixed by
    convention, axis along the k1 direction). `corrected`: that torus
    corrected on the two ring's contact elements (mesh normals) on the 9-dim
    manifold (`lo.correct`); its third-order part is dominated by the mesh
    normals' error. Sign: U holds the larger |curvature|."""
    mode = READOUT["mode"] if mode is None else mode
    N = vertex_normals(m, X) if N is None else N
    ring2 = two_rings(rings(m)) if ring2 is None else ring2
    out = np.zeros((m.n, 6, 6))
    for i in range(m.n):
        ids = ring2[i]
        k1, k2, e1 = quadratic_jet_shape(X, ids, i, N[i])
        s = lo.osculating_torus(X[i], N[i], k1, k2, e1)
        if mode == "corrected":
            s = lo.correct(s, X[ids], N[ids], [1.0] * len(ids), prior, iterations)
        out[i] = canonical_sign(s, X[i], N[i])
    return out


FRAMES = {"kind": "edge"}


def body_frames(m: ops.Mesh, X: np.ndarray, N: np.ndarray, ring: list[np.ndarray], kind: str | None = None,
                ring2=None) -> np.ndarray:
    """Per vertex the Lie lift of its frame (t1, t2, n) at x_i. `edge`: t1 the
    first ring neighbour's direction projected to the tangent plane.
    `principal`: t1 the k1 principal direction of the quadratic jet, its sign
    toward that neighbour. Same convention on every mesh with the same
    connectivity, so a pose-free shape can be re-posed on another mesh."""
    kind = FRAMES["kind"] if kind is None else kind
    if kind == "principal" and ring2 is None:
        ring2 = two_rings(ring)
    B = np.zeros((m.n, 6, 6))
    for i in range(m.n):
        n = N[i]
        d = X[ring[i][0]] - X[i]
        d = d - (d @ n) * n
        if kind == "principal":
            _, _, e1 = quadratic_jet_shape(X, ring2[i], i, n)
            t1 = e1 - (e1 @ n) * n
            if t1 @ d < 0:
                t1 = -t1
        else:
            t1 = d
        t1 /= np.linalg.norm(t1)
        B[i] = lo.rigid(np.stack([t1, np.cross(n, t1), n], 1), X[i])
    return B


def to_body(B: np.ndarray, sig: np.ndarray) -> np.ndarray:
    return np.stack([lo.transform(lo.inverse(b), s) for b, s in zip(B, sig)])


def to_world(B: np.ndarray, sig: np.ndarray) -> np.ndarray:
    return np.stack([lo.transform(b, s) for b, s in zip(B, sig)])


# --------------------------------------------------------------- stage 2


@dataclass
class SolverConfig:
    iterations: int = 60
    w_shape: float = 10.0  # shape Poisson coupling of the posed cyclides (0: off)
    eps: float = 1e-3  # proximal anchor (translation null space), relative
    weights: str = "cotan"  # spoke weights: uniform | cotan (clipped, rest)
    tension: float = 0.0  # surface tension gamma (gamma L_cot in the position matrix)
    # Tangential memory: neighbours are predicted by ARAP (rest spoke turned by
    # the ring's best rotation) and that prediction is landed, so the rest
    # sets where along the surface a neighbour sits and the cyclide where
    # across it. The rest is reseeded every `restart` iterates (0: never).
    arap_predict: bool = False
    restart: int = 10
    land_iterations: int = 6
    # `land`: spokes from landing (closest point, normal moves only).
    # `lie`: the ring's contact elements pushed by the transvection taking
    # the live read-out cyclide to the target (`(sigma_t sigma_c)^(1/2)`), so
    # stretch comes out of the Lie motion; then landed on the target.
    predict: str = "land"
    land_after_lie: bool = True
    alpha: float = 1.0  # fraction of the geodesic per iterate, g^alpha
    readout: str = "corrected"
    # Weak ARAP energy beside the landing (same spoke matrix, so the factor
    # stays constant): holds the tangential placement landing leaves free.
    # Rest reseeded every `restart` iterates.
    w_arap: float = 0.0
    frames: str = "principal"  # body frame tangent: edge | principal
    # Return the iterate with the smallest landing residual (needs no target);
    # stop after `patience` iterates without a new best (0: never stop).
    keep_best: bool = True
    patience: int = 0
    # Under-relaxation of the position update, X <- X + relax (X_new - X):
    # posing the targets on the live mesh makes this a fixed-point iteration
    # (not descent on a fixed energy), which can overshoot.
    relax: float = 1.0
    tol: float = 1e-10


@dataclass
class SolverResult:
    X: np.ndarray
    history: list = field(default_factory=list)  # (max move, landing residual)
    seconds: float = 0.0


class DupinSolver:
    """Positions from per-vertex target cyclides (pose-free, in body frames).

    Each iterate: pose every target on the live mesh (`B_i sigma_body B_i^-1`),
    optionally couple the posed cyclides by the shape Poisson, land each
    ring on its vertex's cyclide (spokes `p_ij - p_ii`), and gather with one
    constant factor `A^T W A + gamma L + eps M` (A the directed spoke
    incidence)."""

    def __init__(self, m: ops.Mesh, X0: np.ndarray, sigma_body: np.ndarray, cfg: SolverConfig = SolverConfig()):
        self.m, self.cfg, self.sigma_body = m, cfg, sigma_body
        self.ring = rings(m)
        pairs = [(i, j) for i in range(m.n) for j in self.ring[i]]
        self.pairs = np.array(pairs)
        if cfg.weights == "cotan":
            wE = np.maximum(ops.cotan_weights(m, X0), 1e-3)
            index = {(min(a, b), max(a, b)): k for k, (a, b) in enumerate(m.E)}
            w = np.array([wE[index[(min(i, j), max(i, j))]] for i, j in pairs])
        else:
            w = np.ones(len(pairs))
        self.w = w
        rows = np.arange(len(pairs))
        A = sp.csr_matrix((np.concatenate([-np.ones(len(pairs)), np.ones(len(pairs))]),
                           (np.concatenate([rows, rows]), np.concatenate([self.pairs[:, 0], self.pairs[:, 1]]))),
                          shape=(len(pairs), m.n))
        self.A = A
        K = (1.0 + cfg.w_arap) * (A.T @ sp.diags(w) @ A)
        self.Mv = ops.vertex_areas(m, X0)
        scale = K.diagonal().mean() / self.Mv.mean()
        self.eps = cfg.eps * scale
        Kt = K + self.eps * sp.diags(self.Mv)
        if cfg.tension > 0:
            Kt = Kt + cfg.tension * positive_laplacian(m, X0)
        self.solve = spla.factorized(Kt.tocsc())
        self.shape = (ShapePoisson(positive_laplacian(m, X0), self.Mv, cfg.w_shape) if cfg.w_shape > 0 else None)
        self.theta = {}
        self.pair_index = [np.nonzero(self.pairs[:, 0] == i)[0] for i in range(m.n)]
        self.ring2 = two_rings(self.ring)
        self.X_rest = X0.copy()

    def predictions(self, X: np.ndarray, force: bool = False) -> np.ndarray:
        """Per directed pair the ARAP prediction of x_j from x_i: the rest
        spoke turned by the ring's best rotation (Procrustes)."""
        pred = X[self.pairs[:, 1]].copy()
        if not (self.cfg.arap_predict or force):
            return pred
        Xr = self.X_rest
        for i in range(self.m.n):
            k = self.pair_index[i]
            d0 = Xr[self.pairs[k, 1]] - Xr[i]
            d1 = X[self.pairs[k, 1]] - X[i]
            U, _, Vt = np.linalg.svd((d1 * self.w[k, None]).T @ d0)
            R = U @ Vt
            if np.linalg.det(R) < 0:
                U[:, -1] *= -1
                R = U @ Vt
            pred[k] = X[i] + d0 @ R.T
        return pred

    def spokes(self, X: np.ndarray, sig: np.ndarray, N: np.ndarray):
        d = X[self.pairs[:, 1]] - X[self.pairs[:, 0]]
        res = 0.0
        conics = []
        for s in sig:
            try:
                conics.append(lo.dupin_conic_basis(s))
            except (ValueError, np.linalg.LinAlgError):
                conics.append(None)  # not a cyclide: the ring keeps its spokes
        foot = np.zeros((self.m.n, 3))
        its = self.cfg.land_iterations
        # No warm start across iterates: the conic bases come from an
        # eigendecomposition and are not continuous in sigma, so old angles
        # point elsewhere on the surface. Each landing starts from the
        # contact element.
        for i in range(self.m.n):
            if conics[i] is None:
                foot[i] = X[i]
                continue
            p, _, th, ok = lo.dupin_land(sig[i], X[i], n0=N[i], conic=conics[i], iterations=its)
            foot[i] = p if ok else X[i]
        pred = self.predictions(X)
        for k, (i, j) in enumerate(self.pairs):
            if conics[i] is None:
                continue
            p, _, th, ok = lo.dupin_land(sig[i], pred[k], n0=N[j], conic=conics[i], iterations=its)
            if not ok:
                p = pred[k]
            d[k] = p - foot[i]
            res = max(res, np.linalg.norm(p - pred[k]))
        self.failed = sum(c is None for c in conics)
        if self.failed:
            res = np.inf  # an iterate with non-cyclides is never the best
        return d, res

    def lie_spokes(self, X: np.ndarray, sig: np.ndarray, N: np.ndarray):
        cur = readout(self.m, X, N, self.ring2, mode=self.cfg.readout)
        d = np.zeros((len(self.pairs), 3))
        res = 0.0
        its = self.cfg.land_iterations
        a = self.cfg.alpha
        for i in range(self.m.n):
            g = lo.transvection(cur[i], sig[i])
            if not np.all(np.isfinite(g)) or lo.group_defect(g) > 1e-6:
                g = np.eye(6)
            target = sig[i]
            if a != 1.0:
                g = lo.lexp(a * lo.llog(g))
                target = lo.transform(g, cur[i])
            conic = lo.dupin_conic_basis(target) if self.cfg.land_after_lie else None
            moved = {}
            for j in [i] + list(self.ring[i]):
                P, T = lo.contact(X[j], N[j])
                pn = lo.contact_decode(g @ P, g @ T)
                p = X[j] if pn is None else pn[0]
                if conic is not None and pn is not None:
                    q, _, _, ok = lo.dupin_land(target, p, n0=pn[1], conic=conic, iterations=its)
                    if ok:
                        p = q
                moved[j] = p
            for k in self.pair_index[i]:
                j = self.pairs[k, 1]
                d[k] = moved[j] - moved[i]
                res = max(res, np.linalg.norm(moved[j] - moved[i] - (X[j] - X[i])))
        return d, res

    def step(self, X: np.ndarray):
        N = vertex_normals(self.m, X)
        B = body_frames(self.m, X, N, self.ring, self.cfg.frames, self.ring2)
        sig = to_world(B, self.sigma_body)
        if self.shape is not None:
            sig = np.stack([canonical_sign(s, X[i], N[i]) for i, s in enumerate(sig)])
            sig = self.shape(sig)
        if self.cfg.predict == "lie":
            d, res = self.lie_spokes(X, sig, N)
        else:
            d, res = self.spokes(X, sig, N)
        if self.cfg.w_arap > 0:
            pred = self.predictions(X, force=True)
            d = d + self.cfg.w_arap * (pred - X[self.pairs[:, 0]])
        rhs = self.A.T @ (self.w[:, None] * d) + self.eps * self.Mv[:, None] * X
        Xn = np.stack([self.solve(rhs[:, c]) for c in range(3)], 1)
        if self.cfg.relax != 1.0:
            Xn = X + self.cfg.relax * (Xn - X)
        return Xn, res

    def run(self, X: np.ndarray, iterations: int | None = None, log=None) -> SolverResult:
        t = time.time()
        out = SolverResult(X.copy())
        self.X_rest = X.copy()
        best, best_res, since = X.copy(), np.inf, 0
        for k in range(iterations or self.cfg.iterations):
            if self.cfg.restart > 0 and k > 0 and k % self.cfg.restart == 0:
                self.X_rest = out.X.copy()
            Xn, res = self.step(out.X)
            # `res` is the landing residual of the iterate the step started from.
            if res < best_res:
                best, best_res, since = out.X.copy(), res, 0
            else:
                since += 1
            move = np.abs(Xn - out.X).max()
            out.X = Xn
            out.history.append((move, res))
            if log:
                log(k, out.X, move, res)
            if move < self.cfg.tol or (self.cfg.patience and since >= self.cfg.patience):
                break
        if self.cfg.keep_best and np.isfinite(best_res):
            out.X = best
        out.best_landing = best_res
        out.seconds = time.time() - t
        return out


def pose_free_targets(m: ops.Mesh, Xt: np.ndarray, mode: str | None = None, frames: str | None = None) -> np.ndarray:
    Nt = vertex_normals(m, Xt)
    return to_body(body_frames(m, Xt, Nt, rings(m), frames), readout(m, Xt, Nt, mode=mode))


# --------------------------------------------------------------- stage 1


def torus_normals(X, R):
    rho = np.hypot(X[:, 0], X[:, 1])
    c = np.stack([R * X[:, 0] / rho, R * X[:, 1] / rho, np.zeros(len(X))], 1)
    d = X - c
    return d / np.linalg.norm(d, axis=1)[:, None]


def stage1(nu=24, nv=12, R=2.0, r=0.7, seed=0):
    rng = np.random.default_rng(seed)
    X, F, _ = ops.torus(R, r, nu, nv)
    m = ops.build_mesh(F, len(X))
    truth = lo.torus(R, r)
    ring = rings(m)
    Mv = ops.vertex_areas(m, X)
    L = positive_laplacian(m, X)
    Lw = np.maximum(ops.cotan_weights(m, X), 0)
    print(f"torus R {R} r {r}, {nu} x {nv} ({m.n} vertices)")

    t = time.time()
    fa = fit_field(m, X, torus_normals(X, R), ring)
    print(f"  fit, analytic normals:  max err {field_error(fa.sigma, truth).max():.1e}  "
          f"gap max {fa.gap.max():.1e}  ok {fa.ok.mean():.2f}  ({time.time() - t:.2f}s)")
    fm = fit_field(m, X, vertex_normals(m, X), ring)
    e_m = field_error(fm.sigma, truth)
    print(f"  fit, mesh normals:      max err {e_m.max():.1e}  mean {e_m.mean():.1e}  ok {fm.ok.mean():.2f}")

    for w in (0.0, 1.0, 10.0, 100.0):
        sp_ = ShapePoisson(L, Mv, w)
        t = time.time()
        out = sp_(fm.sigma)
        dt = time.time() - t
        e = field_error(out, truth)
        print(f"  poisson w {w:6.1f} on mesh-normal fits: max err {e.max():.1e} mean {e.mean():.1e} "
              f"defect {max(lo.defect(s) for s in out):.1e}  ({dt * 1e3:.1f} ms)")

    exact = np.broadcast_to(truth, (m.n, 6, 6)).copy()
    out = ShapePoisson(L, Mv, 10.0)(exact)
    print(f"  exact field fixed point: change {field_error(out, truth).max():.1e}")

    noisy = np.stack([lo.normalize(truth + 1e-2 * rng.standard_normal((6, 6))) for _ in range(m.n)])
    noisy = np.stack([canonical_sign(s, X[i], torus_normals(X[i:i + 1], R)[0]) for i, s in enumerate(noisy)])
    e0 = field_error(noisy, truth)
    for w in (1.0, 10.0, 100.0):
        out = ShapePoisson(L, Mv, w)(noisy)
        e = field_error(out, truth)
        print(f"  noisy field (1e-2) w {w:5.1f}: err {e0.mean():.1e} -> {e.mean():.1e} (max {e.max():.1e})  "
              f"energy {field_energy(m, Lw, noisy):.1e} -> {field_energy(m, Lw, out):.1e}")

    # Gauge: translate the scene by its size, solve, translate back.
    size = np.ptp(X, axis=0).max()
    for shift in (0.5, 1.0, 3.0):
        tvec = shift * size * np.array([1.0, 0.4, -0.3])
        Tn = lo.translation(tvec)
        moved = conjugate_field(Tn, noisy)
        out_t = conjugate_field(lo.inverse(Tn), ShapePoisson(L, Mv, 10.0)(moved))
        out_0 = ShapePoisson(L, Mv, 10.0)(noisy)
        d = field_error(out_t, out_0)
        print(f"  gauge: shift {shift} x size, result change max {d.max():.1e} "
              f"(entries |sigma| {np.abs(moved).max():.0f} vs {np.abs(noisy).max():.0f})")


def degeneracy_census():
    import closed_targets as ct
    print("degeneracy census (one-ring fits, mesh normals)")
    cases = []
    X, F, _ = ops.torus(2.0, 0.7, 24, 12)
    cases.append(("torus 0.7", X, F))
    Xc, Fc = ct.cube(6)
    cases.append(("ellipsoid 1.8/1.2/0.9", ct.radial_ellipsoid(Xc, (1.8, 1.2, 0.9)), Fc))
    Rs = np.sqrt(24.0 / (4 * np.pi))
    cases.append(("sphere", ct.radial_ellipsoid(Xc, (Rs, Rs, Rs)), Fc))
    cases.append(("cube", Xc, Fc))
    P, Fs = ops.strip_param(13, 9, 1.6, 1.2)
    cases.append(("flat strip", ops.flat_from_param(P), Fs))
    cases.append(("cylinder r 1.5", ops.cylinder_from_param(P, 1.5), Fs))
    for name, X, F in cases:
        m = ops.build_mesh(F, len(X))
        N = vertex_normals(m, X)
        ff = fit_field(m, X, N)
        norms = np.abs(ff.sigma).max((1, 2))
        res = np.array([max(lo.contact_residual(ff.sigma[i], X[j], N[j]) for j in [i] + list(rings(m)[i]))
                        for i in range(m.n)])
        print(f"  {name:22s} ok {ff.ok.mean():.2f}  gap median {np.median(ff.gap):.1e} max {ff.gap.max():.1e}  "
              f"|sigma| median {np.median(norms):.1e} max {norms.max():.1e}  ring residual max {res.max():.1e}")


if __name__ == "__main__":
    import sys
    if "census" in sys.argv[1:]:
        degeneracy_census()
    else:
        stage1()
