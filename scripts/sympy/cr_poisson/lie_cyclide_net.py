"""Lie cyclide net: the Lie edge net plus a Dupin cyclide per vertex, tied by
the Poisson Lie-sphere coupling (MATH.md 19.15).

Layers:
- edges (`lie_edge_net`): along sphere, twist, across sphere per edge;
- vertices: a cyclide `sigma_i` composed from the net, not read off the
  mesh: fitted to its own contact element, its ring's contact elements (the
  net's solved normals, consistent with the edge spheres through Rodrigues)
  and the curvatures of its edge-composed tensor `S_i` (19.9: a cyclide read
  independently off the mesh is bias);
- coupling: the ambient Poisson on the cyclide field, `(w L + mu M) sigma =
  mu M sigma_fit` on the 36 entries with one factor, then the matrix sign
  (19.1: no gauge, equivariant under Lie transformations).

Block Gauss-Seidel, each block with the other frozen:
- B (shape field): fit, couple, read each coupled cyclide's curvature
  spheres and principal direction at its vertex (`lie_net.read_vertices`);
- A (geometry): the edge net's Gauss-Newton with two extra rows per cyclide
  vertex, weighted by its anisotropy: `S_i` toward the cyclide's tensor, and
  the ring's contact elements in contact with the cyclide (third order).

Vertices without a cyclide (umbilic or flat: anisotropy under `min_aniso`,
or the fit fails) keep the edge net's tensor and stay out of the coupling.
Flat surfaces are out of scope (every vertex starts degenerate there).
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
import lie_net as ln
import lie_ops as lo

ETA = lo.ETA

# -------------------------------------------------------------- cyclide fit


def tensor_eig(S, T, B):
    """`(k1, k2, e1)` of the 2-tensor `S` in frame `(T, B)`: `k1` the larger
    magnitude, `e1` its unit direction in space."""
    M = np.array([[S[0], S[1]], [S[1], S[2]]])
    ev, V = np.linalg.eigh(M)
    o = np.argsort(-np.abs(ev))
    v = V[:, o[0]]
    return ev[o[0]], ev[o[1]], v[0] * T + v[1] * B


def anisotropy(k1, k2):
    return abs(k1 - k2) / max(abs(k1) + abs(k2), 1e-300)


def _curv_rows(sig, p, n, k1, k2):
    ku, kv, _, _ = lo.curvatures(sig, p, n)
    a, b = sorted((ku, kv), key=lambda k: -abs(k))
    return np.array([a - k1, b - k2])


def fit_vertex_cyclide(start, x, n, ring_x, ring_n, k1, k2, h, w_self=10.0, w_curv=1.0, prior=1e-6,
                       iterations=4):
    """Damped Gauss-Newton in the 9-dim tangent space (Cayley update, as
    `lie_ops.correct`), in coordinates centred on `x` (conditioning). Rows:
    contact at `(x, n)` (weight `w_self`), contact at each ring element, and
    the cyclide's curvatures at `(x, n)` against `(k1, k2)` (times `h`)."""
    Tm, Ti = lo.translation(-x), lo.translation(x)
    d = Tm @ start @ Ti
    pts = [np.zeros(3)] + [q - x for q in ring_x]
    nrm = [n] + list(ring_n)
    w = [w_self] + [1.0] * len(ring_x)
    frames = [lo.contact_frame(p, q) for p, q in zip(pts, nrm)]
    sw = np.sqrt(w_curv) * h
    for _ in range(iterations):
        basis = lo.tangent_basis(d)
        m = len(basis)
        H = prior * np.eye(m)
        g = np.zeros(m)
        dsig = [2 * B @ d for B in basis]
        for wk, f in zip(w, frames):
            r0 = lo.contact_rows(d, f)
            J = np.stack([lo.contact_rows(ds, f) for ds in dsig], 1)
            H += wk * J.T @ J
            g += wk * J.T @ r0
        r0 = _curv_rows(d, np.zeros(3), n, k1, k2)
        if np.all(np.isfinite(r0)):
            eps = 1e-6
            J = np.stack([(_curv_rows(lo.normalize(lo.cayley(eps * B) @ d @ lo.inverse(lo.cayley(eps * B))),
                                      np.zeros(3), n, k1, k2) - r0) / eps for B in basis], 1)
            if np.all(np.isfinite(J)):
                H += sw ** 2 * J.T @ J
                g += sw ** 2 * J.T @ r0
        xs = np.linalg.solve(H, -g)
        c = lo.cayley(sum(xk * B for xk, B in zip(xs, basis)))
        d = lo.normalize(c @ d @ lo.inverse(c))
    return Ti @ d @ Tm


# ------------------------------------------------------- vectorized contact


def contact_residual_rows(sig_loc, d, n):
    """Per row: the contact element `(d, n)` (local coordinates) against the
    cyclide `sig_loc`: for each family (U: `I - sigma`, U^perp: `I + sigma`)
    the vector `Q (alpha P + beta T)` with `(alpha, beta)` the unit minimizer,
    i.e. how far the contact line is from holding a curvature sphere of that
    family. 12 numbers per row, zero exactly when the element touches."""
    nh = n / np.linalg.norm(n, axis=1)[:, None]
    P = lsn.point_rows(d)
    T = lsn.plane_rows(nh, np.einsum("ij,ij->i", nh, d))
    P = P / np.linalg.norm(P, axis=1)[:, None]
    T = T / np.linalg.norm(T, axis=1)[:, None]
    out = []
    I = np.eye(6)
    for s in (1.0, -1.0):
        Q = I[None] - s * sig_loc
        a, b = np.einsum("rij,rj->ri", Q, P), np.einsum("rij,rj->ri", Q, T)
        g11, g12, g22 = (a * a).sum(1), (a * b).sum(1), (b * b).sum(1)
        tr, det = g11 + g22, g11 * g22 - g12 * g12
        lam = 0.5 * (tr - np.sqrt(np.maximum(tr * tr - 4 * det, 0.0)))
        # Eigenvector of [[g11, g12], [g12, g22]] for lam; the better-conditioned
        # of the two forms, sign fixed so the T coefficient is >= 0.
        v1 = np.column_stack([g12, lam - g11])
        v2 = np.column_stack([lam - g22, g12])
        v = np.where((np.linalg.norm(v1, axis=1) >= np.linalg.norm(v2, axis=1))[:, None], v1, v2)
        v = v / np.maximum(np.linalg.norm(v, axis=1), 1e-300)[:, None]
        v = np.where((v[:, 1] < 0)[:, None], -v, v)
        out.append(v[:, 0:1] * a + v[:, 1:2] * b)
    return np.concatenate(out, 1)


# ------------------------------------------------------------------ solver


@dataclass
class CyclideConfig:
    outer: int = 8  # block Gauss-Seidel sweeps
    inner: int = 8  # Gauss-Newton iterates per geometry block
    min_aniso: float = 0.05  # below: no cyclide (umbilic / flat)
    w_shape: float = 10.0  # Poisson coupling of the cyclide field
    w_tensor: float = 1.0  # S_i toward the coupled cyclide's tensor (times aniso)
    # Ring contact with the coupled cyclide (times aniso). Off: at 1 the
    # geometry block wrecks the tori in the second sweep (MATH.md 19.15).
    w_ring: float = 0.0
    fit_iterations: int = 4
    land_iterations: int = 6
    edge: lsn.EdgeNetConfig = field(default_factory=lsn.EdgeNetConfig)


@dataclass
class CyclideField:
    sig: np.ndarray  # (n, 6, 6) world coordinates
    ok: np.ndarray
    aniso: np.ndarray
    A: np.ndarray  # (n, 3, 3) tensor k_u eU eU^T + k_v eV eV^T at the landed element
    loc: np.ndarray  # (n, 6, 6) sigma translated by -x0 (ring rows)
    x0: np.ndarray


class LieCyclideNet(lsn.LieEdgeNet):
    def __init__(self, m: ops.Mesh, X0, targets: lsn.EdgeSpheres, ccfg: CyclideConfig = CyclideConfig(), N0=None):
        super().__init__(m, X0, targets, ccfg.edge, N0)
        self.ccfg = ccfg
        self.ring = dp.rings(m)
        self.field: CyclideField | None = None
        self.sig_prev = None
        self.Lpos = dp.positive_laplacian(m, X0)
        self.Mv = ops.vertex_areas(m, X0)

    # -- block B
    def shape_block(self, z):
        X, N, H, S = self.unpack(z)
        Nn = lsn._nrm(N)
        T, B, _ = lsn.vertex_frames(N, self.R)
        n = self.m.n
        sig = np.zeros((n, 6, 6))
        ok = np.zeros(n, bool)
        for i in range(n):
            k1, k2, e1 = tensor_eig(S[i], T[i], B[i])
            if anisotropy(k1, k2) < self.ccfg.min_aniso or not np.isfinite(k1 + k2):
                continue
            start = (self.sig_prev[i] if self.sig_prev is not None and self.field is not None and self.field.ok[i]
                     else lo.osculating_torus(X[i], Nn[i], k1, k2, e1))
            try:
                s = fit_vertex_cyclide(start, X[i], Nn[i], X[self.ring[i]], Nn[self.ring[i]], k1, k2, self.h,
                                       iterations=self.ccfg.fit_iterations)
            except (np.linalg.LinAlgError, ValueError):
                continue
            if not np.all(np.isfinite(s)) or lo.defect(s) > 1e-6:
                continue
            sig[i] = dp.canonical_sign(s, X[i], Nn[i])
            ok[i] = True
        self.sig_prev = sig.copy()
        sig = self.couple(sig, ok)
        vr = ln.read_vertices(sig, X, Nn, self.ccfg.land_iterations)
        ok &= vr.ok
        A = np.zeros((n, 3, 3))
        eV = np.cross(vr.n, vr.eU)
        A[:] = (vr.ku[:, None, None] * np.einsum("ni,nj->nij", vr.eU, vr.eU)
                + vr.kv[:, None, None] * np.einsum("ni,nj->nij", eV, eV))
        loc = np.stack([lo.translation(-X[i]) @ sig[i] @ lo.translation(X[i]) if ok[i] else np.eye(6)
                        for i in range(n)])
        self.field = CyclideField(sig, ok, np.where(ok, vr.aniso, 0.0), A, loc, X.copy())
        return ok.sum()

    def couple(self, sig, ok):
        """`(w L_ok + mu M) sigma = mu M sigma_fit` on the cyclide vertices,
        `L_ok` keeping only edges with a cyclide at both ends; then the sign."""
        if self.ccfg.w_shape <= 0 or ok.sum() == 0:
            return sig
        E = self.m.E
        keep = ok[E[:, 0]] & ok[E[:, 1]]
        Lc = self.Lpos.tocoo()
        mask = ok[Lc.row] & ok[Lc.col]
        L = sp.csr_matrix((Lc.data[mask], (Lc.row[mask], Lc.col[mask])), shape=Lc.shape)
        # Remove the dropped edges' diagonal share so L stays a Laplacian.
        d = np.asarray(L.sum(1)).ravel()
        L = L - sp.diags(d)
        if keep.sum() == 0:
            return sig
        n = self.m.n
        solve = spla.factorized((self.ccfg.w_shape * L + sp.diags(self.Mv)).tocsc())
        rhs = self.Mv[:, None] * sig.reshape(n, 36)
        out = np.stack([solve(rhs[:, k]) for k in range(36)], 1).reshape(n, 6, 6)
        for i in np.nonzero(ok)[0]:
            out[i] = lo.normalize(out[i])
        out[~ok] = sig[~ok]
        return out

    # -- block A rows
    def _blocks(self, X, N, H, S, L):
        blocks = super()._blocks(X, N, H, S, L)
        f = self.field
        if f is None or not f.ok.any():
            return blocks
        ccfg, o, R, h = self.ccfg, self.off, self.R, self.h
        vi = np.nonzero(f.ok)[0]
        wv = np.sqrt(f.aniso[vi])[:, None]
        Av = f.A[vi]

        def tensor_rows(V):
            n_, S_ = V[:, 0:3], V[:, 3:6]
            T, B, _ = lsn.vertex_frames(n_, R[vi])
            tAt = np.einsum("ri,rij,rj->r", T, Av, T)
            tAb = np.einsum("ri,rij,rj->r", T, Av, B)
            bAb = np.einsum("ri,rij,rj->r", B, Av, B)
            return wv * h * (S_ - np.column_stack([tAt, tAb, bAb]))

        if ccfg.w_tensor > 0:
            blocks.append((np.sqrt(ccfg.w_tensor), tensor_rows, np.concatenate([N[vi], S[vi]], 1),
                           [o["n"] + 3 * vi, o["S"] + 3 * vi]))
        if ccfg.w_ring <= 0:
            return blocks
        # Ring contact (third order): every directed spoke from a cyclide
        # vertex, plus the vertex itself.
        si = np.concatenate([self.spoke[:, 0], np.arange(self.m.n)])
        sj = np.concatenate([self.spoke[:, 1], np.arange(self.m.n)])
        keep = f.ok[si]
        si, sj = si[keep], sj[keep]
        loc, x0 = f.loc[si], f.x0[si]
        wr = np.sqrt(f.aniso[si])[:, None]

        def ring_rows(V):
            return wr * contact_residual_rows(loc, V[:, 0:3] - x0, V[:, 3:6])

        blocks.append((np.sqrt(ccfg.w_ring), ring_rows, np.concatenate([X[sj], N[sj]], 1),
                       [o["x"] + 3 * sj, o["n"] + 3 * sj]))
        return blocks

    def base_merit(self, z):
        f, self.field = self.field, None
        X = self.unpack(z)[0]
        E = self.m.E
        m_ = self.merit(z, np.linalg.norm(X[E[:, 1]] - X[E[:, 0]], axis=1))
        self.field = f
        return m_

    def solve(self, X, N=None, log=None):
        t0 = time.time()
        z = self.initial(X, N)
        hist = [(0.0, self.base_merit(z))]
        best, best_m = z.copy(), hist[0][1]
        facts = backs = 0
        ncyc = []
        for k in range(self.ccfg.outer):
            ncyc.append(int(self.shape_block(z)))
            Xc, Nc = self.unpack(z)[:2]
            r = self.run(Xc, Nc, iterations=self.ccfg.inner, keep_best=False, z0=z)
            facts += r.factorizations
            backs += r.backsolves
            zn = self.pack(r.X, r.N, r.H, r.S)
            move = float(np.abs(r.X - Xc).max())
            z = zn
            mer = self.base_merit(z)
            hist.append((move, mer))
            if mer < best_m:
                best, best_m = z.copy(), mer
            if log:
                log(k, r.X, move, mer)
        out = lsn.EdgeNetResult(*self.unpack(best))
        out.history, out.best_merit = hist, best_m
        out.factorizations, out.backsolves = facts, backs
        out.cyclides = ncyc
        out.seconds = time.time() - t0
        return out
