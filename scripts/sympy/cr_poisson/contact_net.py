"""Principal contact element net as a solver (MATH.md 19.10).

State: one contact element per vertex, positions `x` and normals `n` both
unknown. Each edge holds one sphere, the common sphere of its endpoint
contact elements, through the discrete Rodrigues relation with twist:

    dn = -k dx_t - tau (nbar x dx),   nbar = normalize(n_i + n_j),

`dx_t` the part of `dx` normal to `nbar`. For unit normals `dn . nbar = 0`
exactly, so `k` (the curvature along the edge: the common sphere) and `tau`
(the twist; 0 on curvature-line edges, where the sphere exists exactly) are
the whole of `dn`. Both are read off the target, pose free.

Rows, per coordinate (no row mixes coordinates; the cross terms are lagged
into the right-hand sides), one constant factor on `(x_c, n_c)`:

- Rodrigues, per edge: `k (x_j - x_i) + (n_j - n_i) = k (dx . nbar) nbar
  - tau (nbar x dx)` (current `dx`, `nbar`). Couples positions and normals
  in the global matrix, so how much of a curvature is turning and how much
  is length is decided globally.
- Tangency, per directed spoke: `x_j - x_i` toward its part tangent to `n_i`
  plus the target's normal offset `c_ij |d| n_i`: positions follow normals.
- Unit normals: `n_i` toward `n_i / |n_i|`.
- Length memory (plastic), per edge; anchor `eps M` on positions.
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

# ------------------------------------------------------------- read-outs


def edge_frames(E, X, N):
    dx = X[E[:, 1]] - X[E[:, 0]]
    nb = N[E[:, 0]] + N[E[:, 1]]
    nb /= np.linalg.norm(nb, axis=1)[:, None]
    dxn = np.einsum("ij,ij->i", dx, nb)
    dxt = dx - dxn[:, None] * nb
    b = np.cross(nb, dx)
    return dx, nb, dxn, dxt, b


def edge_readout(E, X, N):
    """`(k, tau)` per edge from `dn = -k dx_t - tau (nbar x dx)`."""
    dx, nb, _, dxt, b = edge_frames(E, X, N)
    dn = N[E[:, 1]] - N[E[:, 0]]
    q = np.einsum("ij,ij->i", dxt, dxt)
    k = -np.einsum("ij,ij->i", dn, dxt) / q
    tau = -np.einsum("ij,ij->i", dn, b) / q
    return k, tau


def spoke_offsets(pairs, X, N):
    d = X[pairs[:, 1]] - X[pairs[:, 0]]
    return np.einsum("ij,ij->i", N[pairs[:, 0]], d) / np.linalg.norm(d, axis=1)


@dataclass
class NetTargets:
    k: np.ndarray
    tau: np.ndarray
    c: np.ndarray  # per directed spoke


def read_targets(m: ops.Mesh, pairs, X, N=None) -> NetTargets:
    N = dp.vertex_normals(m, X) if N is None else N
    k, tau = edge_readout(m.E, X, N)
    return NetTargets(k, tau, spoke_offsets(pairs, X, N))


# ---------------------------------------------------------------- solver


@dataclass
class NetConfig:
    iterations: int = 400
    w_rodrigues: float = 1.0
    w_tangency: float = 1.0
    w_unit: float = 1.0
    # Length memory: weak, reseeded every iterate (a proximal damping, so the
    # map stays a function of the iterate and Anderson keeps its history).
    w_memory: float = 0.01
    restart: int = 1  # memory reseed interval (0: elastic to the start)
    twist: bool = True
    # Two-way tangency: the unit-normal rows target the normal fitted to the
    # ring's tangency offsets (normals follow positions), `normal_mu` its
    # pull toward the current normal.
    normal_fit: bool = True
    normal_mu: float = 1e-2
    eps: float = 1e-6
    patience: int = 0
    tol: float = 1e-12
    accel: str = "anderson"  # none | anderson | chebyshev
    anderson_m: int = 5
    cheb_delay: int = 30  # plain iterates before Chebyshev starts (rho estimate)
    cheb_window: int = 20
    gn_refresh: int = 5  # ContactNetGN: iterates a factor is kept
    # ContactNetGN: the damping memory weight scales with merit / merit_start.
    memory_fade: bool = True


@dataclass
class NetResult:
    X: np.ndarray
    N: np.ndarray
    history: list = field(default_factory=list)  # (move, merit)
    seconds: float = 0.0
    best_merit: float = np.inf


class ContactNetSolver:
    """`k_ref`: the curvatures built into the factor (default the targets).
    With a reference the factor survives target changes: the Rodrigues rows
    are `k_ref dx + dn = ... - (k - k_ref) dx` with the difference lagged, so
    `set_targets` swaps targets without refactoring (frames of an edit)."""

    def __init__(self, m: ops.Mesh, X0: np.ndarray, tgt: NetTargets, cfg: NetConfig = NetConfig(), pairs=None,
                 k_ref: np.ndarray | None = None):
        self.m, self.cfg, self.tgt = m, cfg, tgt
        self.k_mat = (tgt.k if k_ref is None else k_ref).copy()
        ring = dp.rings(m)
        self.pairs = np.array([(i, j) for i in range(m.n) for j in ring[i]]) if pairs is None else pairs
        self.pair_index = [np.nonzero(self.pairs[:, 0] == i)[0] for i in range(m.n)]
        n, nE, nP = m.n, len(m.E), len(self.pairs)
        E = m.E
        self.h = np.linalg.norm(X0[E[:, 1]] - X0[E[:, 0]], axis=1).mean()
        rows, cols, vals = [], [], []

        def add(r, c, v):
            rows.append(r)
            cols.append(c)
            vals.append(v)

        r = 0
        # Rodrigues: sqrt(w_r) (k (x_j - x_i) + (n_j - n_i)); unknowns x in [0, n), n in [n, 2n).
        sr = np.sqrt(cfg.w_rodrigues)
        for e, (i, j) in enumerate(E):
            add(r, j, sr * self.k_mat[e])
            add(r, i, -sr * self.k_mat[e])
            add(r, n + j, sr)
            add(r, n + i, -sr)
            r += 1
        self.r_tan = r
        st = np.sqrt(cfg.w_tangency) / self.h
        for i, j in self.pairs:
            add(r, j, st)
            add(r, i, -st)
            r += 1
        self.r_unit = r
        su = np.sqrt(cfg.w_unit)
        for i in range(n):
            add(r, n + i, su)
            r += 1
        self.r_mem = r
        sm = np.sqrt(cfg.w_memory) / self.h
        for i, j in E:
            add(r, j, sm)
            add(r, i, -sm)
            r += 1
        self.A = sp.csr_matrix((vals, (rows, cols)), shape=(r, 2 * n))
        K = self.A.T @ self.A
        self.Mv = ops.vertex_areas(m, X0)
        dmean = K.diagonal()[:n].mean()
        self.eps_x = cfg.eps * dmean / self.Mv.mean()
        K = K + sp.diags(np.concatenate([self.eps_x * self.Mv, np.full(n, cfg.eps * dmean)]))
        self.solve = spla.factorized(K.tocsc())
        self.kbar = np.mean(np.abs(tgt.k)) + 1e-300
        self.L_rest = None

    def rhs(self, X, N):
        cfg, m, n = self.cfg, self.m, self.m.n
        E = m.E
        dx, nb, dxn, _, b = edge_frames(E, X, N)
        sr = np.sqrt(cfg.w_rodrigues)
        rod = self.tgt.k[:, None] * dxn[:, None] * nb - (self.tgt.k - self.k_mat)[:, None] * dx
        if cfg.twist:
            rod = rod - self.tgt.tau[:, None] * b
        rod = sr * rod
        # Tangency: spoke toward its tangential part plus the target offset.
        d = X[self.pairs[:, 1]] - X[self.pairs[:, 0]]
        ni = N[self.pairs[:, 0]] / np.linalg.norm(N[self.pairs[:, 0]], axis=1)[:, None]
        dl = np.linalg.norm(d, axis=1)
        dstar = d - np.einsum("ij,ij->i", ni, d)[:, None] * ni + (self.tgt.c * dl)[:, None] * ni
        tan = np.sqrt(cfg.w_tangency) / self.h * dstar
        unit = np.sqrt(cfg.w_unit) * self.fitted_normals(X, N)
        le_ = np.linalg.norm(dx, axis=1)
        mem = np.sqrt(cfg.w_memory) / self.h * (self.L_rest / le_)[:, None] * dx
        rows = np.concatenate([rod, tan, unit, mem])
        b_ = self.A.T @ rows
        b_[:n] += self.eps_x * self.Mv[:, None] * X
        b_[n:] += cfg.eps * self.solve_scale * N
        return b_

    def set_targets(self, tgt: NetTargets):
        """New targets on the same factor (the difference to `k_ref` lagged)."""
        self.tgt = tgt
        self.kbar = np.mean(np.abs(tgt.k)) + 1e-300

    def fitted_normals(self, X, N):
        """Normals follow positions: per vertex the unit normal that best meets
        its ring's tangency offsets `n . d_ij / |d_ij| = c_ij` for the current
        spokes, near the current normal (`min |D n - c|^2 + mu |n - n_i|^2`,
        then normalized). Without it a normal can tilt off its surface and
        drag its spokes along (the director-shear mode)."""
        cfg = self.cfg
        out = N / np.linalg.norm(N, axis=1)[:, None]
        if not cfg.normal_fit:
            return out
        d = X[self.pairs[:, 1]] - X[self.pairs[:, 0]]
        dh = d / np.linalg.norm(d, axis=1)[:, None]
        for i in range(self.m.n):
            k = self.pair_index[i]
            D = dh[k]
            A = D.T @ D + cfg.normal_mu * np.eye(3)
            n = np.linalg.solve(A, D.T @ self.tgt.c[k] + cfg.normal_mu * out[i])
            nn = np.linalg.norm(n)
            if nn > 0:
                out[i] = n / nn
        return out

    def merit(self, X, N):
        k, tau = edge_readout(self.m.E, X, N / np.linalg.norm(N, axis=1)[:, None])
        c = spoke_offsets(self.pairs, X, N / np.linalg.norm(N, axis=1)[:, None])
        return float(np.sqrt(np.mean(((k - self.tgt.k) / self.kbar) ** 2) +
                             np.mean(((tau - self.tgt.tau) / self.kbar) ** 2) +
                             np.mean((c - self.tgt.c) ** 2)))

    def G(self, X, N):
        """One projective-dynamics iterate: the fixed-point map `z -> G(z)`."""
        n = self.m.n
        b = self.rhs(X, N)
        z = np.stack([self.solve(b[:, c]) for c in range(3)], 1)
        self.backsolves += 3
        return z[:n], z[n:]

    def run(self, X, N=None, iterations=None, log=None, keep_best=True) -> NetResult:
        """Plain (`accel = none`), Anderson or Chebyshev iteration of `G`.

        Anderson (Peng et al. 2018, type II, window `anderson_m`): the candidate
        `g_k - dG gamma`, `gamma = argmin |f_k - dF gamma|`, `f = G(z) - z`;
        taken only if its merit is not above the plain iterate's, else the
        plain iterate and the history cleared (their safeguard; the plain
        iterate is computed anyway as `g_k`). Chebyshev (Wang 2015):
        `z_{k+1} = w_{k+1} (G(z_k) - z_{k-1}) + z_{k-1}`, `rho` estimated from
        the plain iterates' merit ratio before it starts, reset on a rise.
        The plastic memory reseed clears both."""
        t0 = time.time()
        cfg = self.cfg
        m = self.m
        E = m.E
        n = m.n
        N = dp.vertex_normals(m, X) if N is None else N.copy()
        self.solve_scale = self.A.T.dot(self.A).diagonal()[:n].mean()
        self.L_rest = np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1)
        self.backsolves = 0
        out = NetResult(X.copy(), N.copy())
        best = (X.copy(), N.copy())
        best_m = self.merit(X, N)
        out.history.append((0.0, best_m))
        since = 0
        pack = lambda A, B: np.concatenate([A.ravel(), B.ravel()])
        unpack = lambda z: (z[:3 * n].reshape(n, 3), z[3 * n:].reshape(n, 3))
        hist_f, hist_g = [], []
        z_prev, omega, cheb_on, rho = None, 1.0, False, None
        mer_prev = best_m
        for it in range(iterations or cfg.iterations):
            if cfg.restart > 0 and it > 0 and it % cfg.restart == 0:
                self.L_rest = np.linalg.norm(out.X[E[:, 1]] - out.X[E[:, 0]], axis=1)
                if cfg.restart > 1:  # the map changed; at 1 it is a function of z (pure damping)
                    hist_f, hist_g = [], []
                    omega = 1.0
            Xg, Ng = self.G(out.X, out.N)
            z = pack(out.X, out.N)
            g = pack(Xg, Ng)
            Xn, Nn = Xg, Ng
            if cfg.accel == "anderson":
                hist_f.append(g - z)
                hist_g.append(g)
                hist_f, hist_g = hist_f[-(cfg.anderson_m + 1):], hist_g[-(cfg.anderson_m + 1):]
                if len(hist_f) > 1:
                    dF = np.stack([hist_f[i + 1] - hist_f[i] for i in range(len(hist_f) - 1)], 1)
                    dG = np.stack([hist_g[i + 1] - hist_g[i] for i in range(len(hist_g) - 1)], 1)
                    gam = np.linalg.lstsq(dF, hist_f[-1], rcond=None)[0]
                    Xa, Na = unpack(g - dG @ gam)
                    if np.all(np.isfinite(Xa)) and self.merit(Xa, Na) <= self.merit(Xg, Ng):
                        Xn, Nn = Xa, Na
                    else:
                        hist_f, hist_g = hist_f[-1:], hist_g[-1:]
            elif cfg.accel == "chebyshev":
                if not cheb_on and it >= cfg.cheb_delay:
                    h = np.array([q[1] for q in out.history[-cfg.cheb_window:]])
                    rho = float(np.clip((h[-1] / h[0]) ** (1.0 / (len(h) - 1)), 0.0, 0.9999))
                    cheb_on, omega = True, 1.0
                if cheb_on and z_prev is not None:
                    omega = 2.0 / (2.0 - rho * rho) if omega == 1.0 else 4.0 / (4.0 - rho * rho * omega)
                    Xn, Nn = unpack(omega * (g - z_prev) + z_prev)
            mer = self.merit(Xn, Nn)
            if cfg.accel == "chebyshev" and cheb_on and mer > mer_prev:
                Xn, Nn, omega = Xg, Ng, 1.0
                mer = self.merit(Xn, Nn)
            z_prev = z
            move = np.abs(Xn - out.X).max()
            out.X, out.N = Xn, Nn
            mer_prev = mer
            out.history.append((move, mer))
            if mer < best_m:
                best, best_m, since = (Xn.copy(), Nn.copy()), mer, 0
            else:
                since += 1
            if log:
                log(it, Xn, move, mer)
            if move < cfg.tol or not np.isfinite(mer) or (cfg.patience and since >= cfg.patience):
                break
        out.best_merit = best_m
        out.backsolves = self.backsolves
        out.factorizations = 1
        if keep_best:
            out.X, out.N = best
        out.seconds = time.time() - t0
        return out


# ------------------------------------------- Gauss-Newton on un-lagged rows


def _rodrigues(V, k, tau, twist):
    xi, xj, ni, nj = V[:, 0:3], V[:, 3:6], V[:, 6:9], V[:, 9:12]
    dx, dn = xj - xi, nj - ni
    nb = ni + nj
    nb = nb / np.linalg.norm(nb, axis=1)[:, None]
    r = k[:, None] * (dx - np.einsum("ij,ij->i", dx, nb)[:, None] * nb) + dn
    if twist:
        r = r + tau[:, None] * np.cross(nb, dx)
    return r


def _tangency(V, c):
    xi, xj, ni = V[:, 0:3], V[:, 3:6], V[:, 6:9]
    d = xj - xi
    nh = ni / np.linalg.norm(ni, axis=1)[:, None]
    return (np.einsum("ij,ij->i", nh, d) / np.linalg.norm(d, axis=1) - c)[:, None]


def _unit(V):
    return (np.linalg.norm(V, axis=1) - 1.0)[:, None]


def _length(V, L):
    return (np.linalg.norm(V[:, 3:6] - V[:, 0:3], axis=1) - L)[:, None]


class ContactNetGN:
    """Levenberg-Marquardt on the contact net's rows without lagging: the
    Rodrigues rows with their `nbar` and twist, scalar tangency
    `n_i . d_ij / |d_ij| - c_ij` (normals and positions both, replacing the
    normal fit), unit length, and length memory reseeded every iterate (a
    proximal damping). Jacobian by vectorized central differences over each
    row's local unknowns; the factor of `J^T J + lambda diag` is kept
    `refresh` iterates with a fresh gradient (the hepworth pattern)."""

    def __init__(self, m, X0, tgt: NetTargets, cfg: NetConfig = NetConfig(), pairs=None,
                 refresh: int | None = None, lam0: float = 1e-4, fd: float = 1e-7):
        refresh = cfg.gn_refresh if refresh is None else refresh
        self.pd = ContactNetSolver(m, X0, tgt, cfg, pairs)  # merit, pairs, scales
        self.m, self.cfg, self.tgt = m, cfg, tgt
        self.pairs, self.h = self.pd.pairs, self.pd.h
        self.refresh, self.lam0, self.fd = refresh, lam0, fd
        self.memory_scale = 1.0

    def _blocks(self, X, N, L):
        m, n, P, E = self.m, self.m.n, self.pairs, self.m.E
        cfg, tgt = self.cfg, self.tgt
        xo = lambda v: 3 * v
        no = lambda v: 3 * n + 3 * v
        sr, st, su, sm = (np.sqrt(cfg.w_rodrigues), np.sqrt(cfg.w_tangency), np.sqrt(cfg.w_unit),
                          np.sqrt(cfg.w_memory * self.memory_scale) / self.h)
        return [
            (sr, lambda V: _rodrigues(V, tgt.k, tgt.tau, cfg.twist),
             np.concatenate([X[E[:, 0]], X[E[:, 1]], N[E[:, 0]], N[E[:, 1]]], 1),
             [xo(E[:, 0]), xo(E[:, 1]), no(E[:, 0]), no(E[:, 1])]),
            (st, lambda V: _tangency(V, tgt.c),
             np.concatenate([X[P[:, 0]], X[P[:, 1]], N[P[:, 0]]], 1),
             [xo(P[:, 0]), xo(P[:, 1]), no(P[:, 0])]),
            (su, _unit, N.copy(), [no(np.arange(n))]),
            (sm, lambda V: _length(V, L), np.concatenate([X[E[:, 0]], X[E[:, 1]]], 1),
             [xo(E[:, 0]), xo(E[:, 1])]),
        ]

    def residual(self, X, N, L):
        return np.concatenate([(w * f(V)).ravel() for w, f, V, _ in self._blocks(X, N, L)])

    def jacobian(self, X, N, L):
        rows, cols, vals = [], [], []
        r0 = 0
        nvar = 6 * self.m.n
        for w, f, V, starts in self._blocks(X, N, L):
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
                             shape=(r0, nvar))

    def run(self, X, N=None, iterations=200, keep_best=True, log=None) -> NetResult:
        t0 = time.time()
        m, n, E = self.m, self.m.n, self.m.E
        N = dp.vertex_normals(m, X) if N is None else N.copy()
        out = NetResult(X.copy(), N.copy())
        out.history.append((0.0, self.pd.merit(X, N)))
        best, best_m = (X.copy(), N.copy()), out.history[0][1]
        lam = self.lam0
        solve, age = None, 0
        out.factorizations = out.backsolves = 0
        mer0 = max(out.history[0][1], 1e-300)
        for it in range(iterations):
            # Fading damping memory: full weight while the model is poor, gone
            # at the end (a proximal term on every step otherwise dominates the
            # soft bending modes and makes convergence linear).
            # Updated only when a refactor is due, so a kept factor and its
            # residual agree.
            if self.cfg.memory_fade and (solve is None or age >= self.refresh):
                self.memory_scale = min(1.0, out.history[-1][1] / mer0)
            L = np.linalg.norm(out.X[E[:, 1]] - out.X[E[:, 0]], axis=1)
            r = self.residual(out.X, out.N, L)
            J = self.jacobian(out.X, out.N, L)
            g = J.T @ r
            cost = 0.5 * r @ r
            if solve is None or age >= self.refresh:
                H = (J.T @ J).tocsc()
                dg = H.diagonal()
                solve = spla.factorized((H + sp.diags(lam * dg + 1e-12 * dg.mean())).tocsc())
                out.factorizations += 1
                age = 0
            ok = False
            for _ in range(10):
                dz = -solve(g)
                out.backsolves += 1
                Xn = out.X + dz[:3 * n].reshape(n, 3)
                Nn = out.N + dz[3 * n:].reshape(n, 3)
                rn = self.residual(Xn, Nn, L)
                if np.all(np.isfinite(rn)) and 0.5 * rn @ rn < cost:
                    ok = True
                    break
                lam *= 4.0
                H = (J.T @ J).tocsc()
                dg = H.diagonal()
                solve = spla.factorized((H + sp.diags(lam * dg + 1e-12 * dg.mean())).tocsc())
                out.factorizations += 1
                age = 0
            if not ok:
                break
            lam = max(lam / 3.0, 1e-12)
            age += 1
            move = np.abs(Xn - out.X).max()
            out.X, out.N = Xn, Nn
            mer = self.pd.merit(Xn, Nn)
            out.history.append((move, mer))
            if mer < best_m:
                best, best_m = (Xn.copy(), Nn.copy()), mer
            if log:
                log(it, Xn, move, mer)
            if move < self.cfg.tol:
                break
        out.best_merit = best_m
        if keep_best:
            out.X, out.N = best
        out.seconds = time.time() - t0
        return out


# ------------------------------------------- Lie sphere rows (bilinear)


def _sphere_on_point(V, k):
    """`<s_ij, P(x_j)> / |d|` with `s_ij = k P(x_i) + T(x_i, n_i)`:
    `(n_i . d - k |d|^2 / 2) / |d|`; `j` on the sphere of normal curvature `k`
    that `i` holds toward it (the tangency offset, exact on any surface)."""
    xi, xj, ni = V[:, 0:3], V[:, 3:6], V[:, 6:9]
    d = xj - xi
    dd = np.einsum("ij,ij->i", d, d)
    return ((np.einsum("ij,ij->i", ni, d) - 0.5 * k * dd) / np.sqrt(dd))[:, None]


def _sphere_on_plane(V, k, target, scale):
    """`<s_ij, T(x_j, n_j)> = (n_i . n_j - 1) - k n_j . d`: `j`'s tangent plane
    touches the sphere too (0 when the contact elements share it: a
    curvature-line edge; `-tau^2 |d|^2 / 2` to leading order otherwise)."""
    xi, xj, ni, nj = V[:, 0:3], V[:, 3:6], V[:, 6:9], V[:, 9:12]
    d = xj - xi
    r = np.einsum("ij,ij->i", ni, nj) - 1.0 - k * np.einsum("ij,ij->i", nj, d)
    return ((r - target) * scale)[:, None]


@dataclass
class SphereTargets:
    k: np.ndarray  # per directed spoke: normal curvature toward j
    plane: np.ndarray  # per directed spoke: <s_ij, T_j> at the target


def read_sphere_targets(m, pairs, X, N=None) -> SphereTargets:
    N = dp.vertex_normals(m, X) if N is None else N
    d = X[pairs[:, 1]] - X[pairs[:, 0]]
    dd = np.einsum("ij,ij->i", d, d)
    k = 2.0 * np.einsum("ij,ij->i", N[pairs[:, 0]], d) / dd
    plane = np.einsum("ij,ij->i", N[pairs[:, 0]], N[pairs[:, 1]]) - 1.0 - k * np.einsum("ij,ij->i", N[pairs[:, 1]], d)
    return SphereTargets(k, plane)


class ContactNetSphereGN(ContactNetGN):
    """The net's Gauss-Newton with Lie sphere rows instead of the Euclidean
    Rodrigues / tangency rows: per directed spoke the sphere
    `s_ij = k_ij P(x_i) + T(x_i, n_i)` and the two bilinear rows
    `<s_ij, P(x_j)> = 0`, `<s_ij, T(x_j, n_j)> = target`; unit normals; the
    fading length memory. `tgt` (Rodrigues targets) is kept for the merit."""

    def __init__(self, m, X0, tgt: NetTargets, sph: SphereTargets, cfg: NetConfig = NetConfig(), pairs=None,
                 w_plane: float = 1.0, **kw):
        super().__init__(m, X0, tgt, cfg, pairs, **kw)
        self.sph, self.w_plane = sph, w_plane
        d = X0[self.pairs[:, 1]] - X0[self.pairs[:, 0]]
        kbar = np.mean(np.abs(sph.k)) + 1e-300
        self.plane_scale = 1.0 / (np.einsum("ij,ij->i", d, d) * kbar)

    def _blocks(self, X, N, L):
        m, n, P, E = self.m, self.m.n, self.pairs, self.m.E
        cfg, sph = self.cfg, self.sph
        xo = lambda v: 3 * v
        no = lambda v: 3 * n + 3 * v
        st, su = np.sqrt(cfg.w_tangency), np.sqrt(cfg.w_unit)
        sp_ = np.sqrt(cfg.w_rodrigues * self.w_plane)
        sm = np.sqrt(cfg.w_memory * self.memory_scale) / self.h
        return [
            (st, lambda V: _sphere_on_point(V, sph.k),
             np.concatenate([X[P[:, 0]], X[P[:, 1]], N[P[:, 0]]], 1),
             [xo(P[:, 0]), xo(P[:, 1]), no(P[:, 0])]),
            (sp_, lambda V: _sphere_on_plane(V, sph.k, sph.plane, self.plane_scale),
             np.concatenate([X[P[:, 0]], X[P[:, 1]], N[P[:, 0]], N[P[:, 1]]], 1),
             [xo(P[:, 0]), xo(P[:, 1]), no(P[:, 0]), no(P[:, 1])]),
            (su, _unit, N.copy(), [no(np.arange(n))]),
            (sm, lambda V: _length(V, L), np.concatenate([X[E[:, 0]], X[E[:, 1]]], 1),
             [xo(E[:, 0]), xo(E[:, 1])]),
        ]


# ------------------------------------------------- stage 2: cyclides on the net


def dupin_residuals(m: ops.Mesh, X, N, prior=1e-6, iterations=6):
    """Per vertex: fit a cyclide (`lo.correct` from the osculating torus) to
    the ring's contact elements from the net state, return the RMS
    membership residual of the ring (0 where the ring lies on one cyclide)."""
    ring = dp.rings(m)
    ring2 = dp.two_rings(ring)
    Nn = N / np.linalg.norm(N, axis=1)[:, None]
    res = np.zeros(m.n)
    for i in range(m.n):
        ids = np.union1d(ring[i], [i])
        k1, k2, e1 = dp.quadratic_jet_shape(X, ring2[i], i, Nn[i])
        s = lo.osculating_torus(X[i], Nn[i], k1, k2, e1)
        s = lo.correct(s, X[ids], Nn[ids], [1.0] * len(ids), prior, iterations)
        res[i] = np.sqrt(np.mean([lo.contact_residual(s, X[j], Nn[j]) ** 2 for j in ids]))
    return res
