"""Goal energy on the CR tangent field: E = 1/2 sum_f A_f |S_f(T) - s_g,f|^2.

`S_f(T)` is the shape operator of Stein et al. 2020 (App. A.2) read entirely
from the tangent field `T` (per edge `[g_par, g_perp]`), including the face
normal and frame, which come from the face's `g_par` (oriented edge vectors
`v_i = s_i l0_i g_par`). With `T = D X` this is exactly the read-out of the
embedding; with `T` free, tilting a face is a change of `T`, so `S` can be
inverted in tangent space.

Solvers:
- `gn_positions`: exact Gauss-Newton on `E(D X)` in positions, Levenberg
  damped by the CR Dirichlet matrix `L = D^T M D` (reference).
- `tangent_loop`: local / global.
    T0 = R(X) G0           continual ARAP: rest tangents rotated per edge
    T* = argmin |T - T0|_M  s.t.  S(T) = s_g     (Newton-Raphson, min-norm)
    X  = L^{-1} D^T M T*    (L prefactored, rebuilt when E rises)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_ops as ops
import cr_loop as cl


# ------------------------------------------------------------ the read-out


@dataclass
class Topology:
    fe: np.ndarray  # (m, 3)
    si: np.ndarray  # (m, 3) orientation of local edge i against E[e]
    l0: np.ndarray  # (nE,) rest lengths
    nE: int


def topology(m: ops.Mesh, X0: np.ndarray) -> Topology:
    si = np.zeros(m.fe.shape)
    for f, t in enumerate(m.F):
        for i in range(3):
            si[f, i] = 1.0 if t[i] < t[(i + 1) % 3] else -1.0
    l0 = np.linalg.norm(X0[m.E[:, 1]] - X0[m.E[:, 0]], axis=1)
    return Topology(fe=m.fe, si=si, l0=l0, nE=len(m.E))


def face_S(tp: Topology, Tf: np.ndarray) -> np.ndarray:
    """`(k11, k12, k22)` per face from its edges' tangents `Tf` (m, 3, 6)."""
    gpar, gperp = Tf[:, :, :3], Tf[:, :, 3:]
    v = (tp.si * tp.l0[tp.fe])[:, :, None] * gpar
    n = (np.cross(v[:, 0], v[:, 1]) + np.cross(v[:, 1], v[:, 2])
         + np.cross(v[:, 2], v[:, 0])) / 3.0
    nn = np.linalg.norm(n, axis=1)
    A = 0.5 * nn
    N = n / nn[:, None]
    t = v[:, 0] - np.sum(v[:, 0] * N, 1)[:, None] * N
    t /= np.linalg.norm(t, axis=1)[:, None]
    b = np.cross(N, t)
    s = np.zeros((len(Tf), 3))
    for i in range(3):
        li = np.linalg.norm(v[:, i], axis=1)
        vp = np.cross(N, v[:, i])
        vp /= np.linalg.norm(vp, axis=1)[:, None]
        c = np.sum(N * gperp[:, i], 1)
        k = ops.K_CR_SHAPE_SIGN * tp.si[:, i] * li / A * c
        tv, bv = np.sum(t * vp, 1), np.sum(b * vp, 1)
        s += k[:, None] * np.stack([tv * tv, tv * bv, bv * bv], 1)
    return s


def S_of(tp: Topology, T: np.ndarray) -> np.ndarray:
    return face_S(tp, T[tp.fe])


def jac_S(tp: Topology, T: np.ndarray, h: float = 1e-7) -> sp.csr_matrix:
    """Sparse `dS / dT` (3 m x 6 nE), central differences per face."""
    Tf = T[tp.fe]
    m = len(Tf)
    blocks = np.zeros((m, 3, 18))
    for d in range(18):
        i, c = divmod(d, 6)
        Tp, Tm = Tf.copy(), Tf.copy()
        Tp[:, i, c] += h
        Tm[:, i, c] -= h
        blocks[:, :, d] = (face_S(tp, Tp) - face_S(tp, Tm)) / (2 * h)
    rows = np.repeat(3 * np.arange(m)[:, None] + np.arange(3)[None, :], 18, 1)
    cols = (6 * tp.fe[:, :, None] + np.arange(6)[None, None, :]).reshape(m, 18)
    cols = np.repeat(cols[:, None, :], 3, 1).reshape(m, 54)
    return sp.csr_matrix((blocks.reshape(m, 54).ravel(),
                          (rows.ravel(), cols.ravel())),
                         shape=(3 * m, 6 * tp.nE))


# --------------------------------------------------------- packing / D, L


def tangents(D0, X: np.ndarray, nE: int) -> np.ndarray:
    G = D0 @ X
    return np.concatenate([G[:nE], G[nE:]], axis=1)


def d_matrix(D0, nE: int, n: int) -> sp.csr_matrix:
    """`vec(T) = Dx vec(X)`, `vec` row-major: T (nE, 6), X (n, 3)."""
    D = D0.tocoo()
    e = np.where(D.row < nE, D.row, D.row - nE)
    w = (D.row >= nE).astype(int)
    rows = np.concatenate([6 * e + 3 * w + a for a in range(3)])
    cols = np.concatenate([3 * D.col + a for a in range(3)])
    vals = np.tile(D.data, 3)
    return sp.csr_matrix((vals, (rows, cols)), shape=(6 * nE, 3 * n))


def arap_tangents(G0: np.ndarray, G: np.ndarray) -> np.ndarray:
    """Continual ARAP: per edge `R_e = polar(G_e G0_e^T)`, returns `R_e G0_e`."""
    out = np.zeros_like(G0)
    for e in range(len(G0)):
        A = np.column_stack([G0[e, :3], G0[e, 3:]])
        B = np.column_stack([G[e, :3], G[e, 3:]])
        U, _, Vt = np.linalg.svd(B @ A.T)
        R = U @ Vt
        if np.linalg.det(R) < 0:
            U[:, -1] *= -1
            R = U @ Vt
        out[e, :3] = R @ G0[e, :3]
        out[e, 3:] = R @ G0[e, 3:]
    return out


def face_arap_tangents(m: ops.Mesh, X0: np.ndarray, X: np.ndarray, Me):
    """Continual ARAP per face: `sum_f R_f D_f X0`, `R_f` the face's Procrustes
    rotation. Equals `D X` whenever every face moves rigidly, so isometric
    bending is a fixed point (the per-edge `R_e G0_e` keeps `|g_perp| = 1`,
    while the CR `g_perp` of a bent edge shortens by about cos(theta / 2))."""
    import face_loops as fl
    R = fl._face_procrustes(m, X0, X)
    F = fl._cr_rotated_frames(m, X0, R, Me)
    nE = len(m.E)
    return np.concatenate([F[:nE], F[nE:]], axis=1)


def _polar3(M: np.ndarray) -> np.ndarray:
    U, _, Vt = np.linalg.svd(M)
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    return R


def vertex_rotations(m: ops.Mesh, X0: np.ndarray, X: np.ndarray, w: np.ndarray):
    """Vertex ARAP (cotan-weighted one-ring spokes): Procrustes `R_v` and the
    right-hand side `b` with `L_cot X = b` (`ops.arap_rhs_and_rotations`)."""
    return ops.arap_rhs_and_rotations(m, X0, X, w)


def vertex_arap_tangents(m: ops.Mesh, Rv: np.ndarray, G0: np.ndarray,
                         push: str = "polar") -> np.ndarray:
    """Tangents live with the vertex rotations; pushed to edge `e = (i, j)` as
    `R_e G0_e` with `R_e = polar((R_i + R_j) / 2)` ("polar"), or the plain
    average `(R_i + R_j) / 2` applied to `G0_e` ("blend")."""
    i, j = m.E[:, 0], m.E[:, 1]
    A = 0.5 * (Rv[i] + Rv[j])
    if push == "polar":
        A = np.array([_polar3(a) for a in A])
    return np.concatenate([np.einsum("eab,eb->ea", A, G0[:, :3]),
                           np.einsum("eab,eb->ea", A, G0[:, 3:])], axis=1)


def vertex_arap_energy(m: ops.Mesh, X0, X, Rv, w) -> float:
    """`1/2 sum_e w_e (|d - R_i d0|^2 + |d - R_j d0|^2) / 2`; gradient
    `L_cot X - b` at the Procrustes `R` (envelope), Hessian `L_cot`."""
    i, j = m.E[:, 0], m.E[:, 1]
    d = X[j] - X[i]
    d0 = X0[j] - X0[i]
    ri = d - np.einsum("eab,eb->ea", Rv[i], d0)
    rj = d - np.einsum("eab,eb->ea", Rv[j], d0)
    return 0.25 * float(np.sum(w * (np.sum(ri * ri, 1) + np.sum(rj * rj, 1))))


# ------------------------------------------------- edge rotations (CR element)


def _polar_batch(C: np.ndarray) -> np.ndarray:
    U, _, Vt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(U @ Vt))
    U[:, :, -1] *= d[:, None]
    return U @ Vt


@dataclass
class EdgeRotRest:
    """Rest data for edge rotations, at the ARAP rest `X0`."""
    X0: np.ndarray
    G0: np.ndarray  # (nE, 6) rest tangents
    wfi: np.ndarray  # (m, 3) half cotan of the angle opposite local edge i
    tails: np.ndarray  # (m, 3)
    heads: np.ndarray  # (m, 3)
    fe: np.ndarray  # (m, 3) edge of local edge i
    EF: sp.csr_matrix  # (nE, m) edge-face incidence
    Lc: sp.csr_matrix  # kron(L_cot, I3)
    pairs: np.ndarray  # (k, 2) edge pairs of a face (diamond neighbors)
    wp: np.ndarray  # (k,) 2 |cot| at the shared vertex (vector Dirichlet block norm)
    wself: np.ndarray  # (nE,) vector Dirichlet diagonal
    Lw: sp.csr_matrix  # kron(L_w, I6), L_w the edge-graph Laplacian on `pairs`


def edge_rot_rest(m: ops.Mesh, X0: np.ndarray, D0) -> EdgeRotRest:
    nE, nF = len(m.E), len(m.F)
    tails = m.F
    heads = np.roll(m.F, -1, axis=1)
    opp = np.roll(m.F, -2, axis=1)
    u = X0[tails] - X0[opp]
    v = X0[heads] - X0[opp]
    wfi = 0.5 * np.sum(u * v, 2) / np.linalg.norm(np.cross(u, v), axis=2)
    EF = sp.csr_matrix((np.ones(3 * nF), (m.fe.ravel(), np.repeat(np.arange(nF), 3))),
                       shape=(nE, nF))
    Lc = sp.kron(ops.cotan_laplacian(m, X0), sp.identity(3)).tocsr()
    Ld = ops.vector_dirichlet(m, X0).tocsr()
    pairs = []
    for f in range(nF):
        for i in range(3):
            pairs.append((m.fe[f, i], m.fe[f, (i + 1) % 3]))
    pairs = np.array(pairs)
    a, b = pairs[:, 0], pairs[:, 1]
    wp = np.hypot(np.asarray(Ld[a, b]).ravel(), np.asarray(Ld[a, b + nE]).ravel())
    wself = Ld.diagonal()[:nE]
    rows = np.concatenate([a, b, a, b])
    cols = np.concatenate([b, a, a, b])
    vals = np.concatenate([-wp, -wp, wp, wp])
    Lw = sp.csr_matrix((vals, (rows, cols)), shape=(nE, nE))
    return EdgeRotRest(X0=X0, G0=tangents(D0, X0, nE), wfi=wfi, tails=tails,
                       heads=heads, fe=m.fe, EF=EF, Lc=Lc, pairs=pairs, wp=wp,
                       wself=wself, Lw=sp.kron(Lw, sp.identity(6)).tocsr())


def cr_arap_rotations(er: EdgeRotRest, X: np.ndarray) -> np.ndarray:
    """CR-ARAP: `R_e` the Procrustes fit over the diamond of `e` (its faces'
    cotan-weighted spokes; the center edge gets its full cotan weight)."""
    d = X[er.heads] - X[er.tails]
    d0 = er.X0[er.heads] - er.X0[er.tails]
    Kf = np.einsum("fi,fia,fib->fab", er.wfi, d, d0).reshape(-1, 9)
    return _polar_batch(np.asarray(er.EF @ Kf).reshape(-1, 3, 3))


def cr_arap_energy_grad(er: EdgeRotRest, X: np.ndarray, R: np.ndarray):
    """`E = sum_e sum_{f ni e} 1/3 sum_{ij in f} w_ij |d - R_e d0|^2`; gradient
    `L_cot X - b` at the Procrustes `R` (envelope), Hessian `L_cot`."""
    d = X[er.heads] - X[er.tails]
    d0 = er.X0[er.heads] - er.X0[er.tails]
    E = 0.0
    for k in range(3):
        Rd0 = np.einsum("fab,fib->fia", R[er.fe[:, k]], d0)
        E += np.sum(er.wfi[:, :, None] * (d - Rd0) ** 2) / 3.0
    Rbar = R[er.fe].mean(1)
    r = er.wfi[:, :, None] * np.einsum("fab,fib->fia", Rbar, d0)
    b = np.zeros_like(X)
    np.add.at(b, er.heads.ravel(), r.reshape(-1, 3))
    np.add.at(b, er.tails.ravel(), -r.reshape(-1, 3))
    return 0.5 * E, er.Lc @ X.ravel() - b.ravel()


def dir_rotations(er: EdgeRotRest, G: np.ndarray, spokes: bool) -> np.ndarray:
    """Vector-Dirichlet-weighted Procrustes on edge frames. `spokes=False`:
    `polar(w_ee G_e G0_e^T + sum_{e'} w_ee' G_e' G0_e'^T)`; `spokes=True`:
    on frame differences `(G_e' - G_e)(G0_e' - G0_e)^T` (the `dir_arap` fit)."""
    def outer(A, B):
        return (np.einsum("ea,eb->eab", A[:, :3], B[:, :3])
                + np.einsum("ea,eb->eab", A[:, 3:], B[:, 3:]))
    a, b = er.pairs[:, 0], er.pairs[:, 1]
    nE = len(G)
    C = np.zeros((nE, 3, 3))
    if spokes:
        P = er.wp[:, None, None] * outer(G[b] - G[a], er.G0[b] - er.G0[a])
        np.add.at(C, a, P)
        np.add.at(C, b, P)
    else:
        C += er.wself[:, None, None] * outer(G, er.G0)
        np.add.at(C, a, er.wp[:, None, None] * outer(G[b], er.G0[b]))
        np.add.at(C, b, er.wp[:, None, None] * outer(G[a], er.G0[a]))
    return _polar_batch(C)


def dir_arap_energy_grad_T(er: EdgeRotRest, G: np.ndarray, R: np.ndarray):
    """`E = 1/4 sum_pairs w (|D - R_a D0|^2 + |D - R_b D0|^2)`, `D = G_b - G_a`;
    gradient in `vec(T)`: `L_w t - b` (envelope), Hessian `L_w` (x I6)."""
    a, b = er.pairs[:, 0], er.pairs[:, 1]
    Dl = (G[b] - G[a]).reshape(-1, 2, 3)
    D0 = (er.G0[b] - er.G0[a]).reshape(-1, 2, 3)
    ra = Dl - np.einsum("kab,kcb->kca", R[a], D0)
    rb = Dl - np.einsum("kab,kcb->kca", R[b], D0)
    E = 0.25 * float(np.sum(er.wp[:, None, None] * (ra ** 2 + rb ** 2)))
    r = 0.5 * er.wp[:, None, None] * np.einsum("kab,kcb->kca", R[a] + R[b], D0)
    bt = np.zeros((len(G), 2, 3))
    np.add.at(bt, b, r)
    np.add.at(bt, a, -r)
    return E, er.Lw @ G.ravel() - bt.reshape(len(G), 6).ravel()


def rotate_tangents(R: np.ndarray, G0: np.ndarray) -> np.ndarray:
    return np.concatenate([np.einsum("eab,eb->ea", R, G0[:, :3]),
                           np.einsum("eab,eb->ea", R, G0[:, 3:])], axis=1)


# --------------------------- intrinsic vector Dirichlet memory (Stein L_dir)


def intrinsic_frames(T: np.ndarray) -> np.ndarray:
    """`U_e = [t, b]^T [g_par, g_perp]` (nE, 2, 2) in each edge's own frame
    `t = g_par / |g_par|`, `n = g_par x g_perp / |.|`, `b = n x t`, which is
    the Stein edge frame (`g_perp` lies along `N x t`). Rotation-invariant;
    an isometric bend changes it only through the `|g_perp|` shortening."""
    gp, gq = T[:, :3], T[:, 3:]
    t = gp / np.linalg.norm(gp, axis=1)[:, None]
    n = np.cross(gp, gq)
    n /= np.linalg.norm(n, axis=1)[:, None]
    b = np.cross(n, t)
    U = np.empty((len(T), 2, 2))
    U[:, 0, 0] = np.sum(t * gp, 1)
    U[:, 1, 0] = np.sum(b * gp, 1)
    U[:, 0, 1] = np.sum(t * gq, 1)
    U[:, 1, 1] = np.sum(b * gq, 1)
    return U


@dataclass
class DirIntrinsic:
    U0: np.ndarray  # (nE, 2, 2) rest intrinsic frames
    K: sp.csr_matrix  # (4 nE)^2: blockdiag over the two columns of L_dir + mu M2


def dir_intrinsic_rest(m: ops.Mesh, X_rest: np.ndarray, G0: np.ndarray,
                       Me: np.ndarray, mu_rel: float = 1e-3) -> DirIntrinsic:
    """`L_dir = ops.vector_dirichlet(m, X_rest)` as written (weights, block
    rotations, orientation signs, diagonal) acting on intrinsic 2-vectors in
    the `[par (nE), perp (nE)]` layout, plus `mu M2` for its null space."""
    Ld = ops.vector_dirichlet(m, X_rest).tocsr()
    M2 = sp.diags(np.concatenate([Me, Me]))
    mu = mu_rel * Ld.diagonal().sum() / M2.diagonal().sum()
    Kc = (Ld + mu * M2).tocsr()
    return DirIntrinsic(U0=intrinsic_frames(G0), K=sp.block_diag([Kc, Kc]).tocsr())


def _dir_stack(dU: np.ndarray) -> np.ndarray:
    """(nE, 2, 2) -> c = [col 0 par, col 0 perp, col 1 par, col 1 perp]."""
    return np.concatenate([dU[:, 0, 0], dU[:, 1, 0], dU[:, 0, 1], dU[:, 1, 1]])


def dir_intrinsic_energy(di: DirIntrinsic, T: np.ndarray) -> float:
    c = _dir_stack(intrinsic_frames(T) - di.U0)
    return 0.5 * float(c @ (di.K @ c))


def dir_intrinsic_jacobian(T: np.ndarray, h: float = 1e-7) -> sp.csr_matrix:
    """`dc / dvec(T)` (4 nE x 6 nE): `U_e` depends on `T_e` only, so the
    Jacobian is one 4 x 6 block per edge (central differences)."""
    nE = len(T)
    blocks = np.zeros((nE, 4, 6))
    for d in range(6):
        Tp, Tm = T.copy(), T.copy()
        Tp[:, d] += h
        Tm[:, d] -= h
        dU = (intrinsic_frames(Tp) - intrinsic_frames(Tm)) / (2 * h)
        blocks[:, :, d] = np.stack([dU[:, 0, 0], dU[:, 1, 0], dU[:, 0, 1], dU[:, 1, 1]], 1)
    e = np.arange(nE)
    rows = (np.arange(4)[None, :] * nE + e[:, None])  # (nE, 4)
    rows = np.repeat(rows[:, :, None], 6, 2)
    cols = np.repeat((6 * e[:, None] + np.arange(6)[None, :])[:, None, :], 4, 1)
    return sp.csr_matrix((blocks.ravel(), (rows.ravel(), cols.ravel())),
                         shape=(4 * nE, 6 * nE))


def dir_intrinsic_energy_grad(di: DirIntrinsic, T: np.ndarray, Dx):
    """`E = 1/2 c^T K c`, `c = stack(U(D X) - U0)`; gradient in `X`."""
    c = _dir_stack(intrinsic_frames(T) - di.U0)
    Jc = dir_intrinsic_jacobian(T) @ Dx
    return 0.5 * float(c @ (di.K @ c)), Jc.T @ (di.K @ c)


def goal_energy(tp, T, s3g, area) -> float:
    return 0.5 * float(np.sum(area[:, None] * (S_of(tp, T) - s3g) ** 2))


@dataclass
class Result:
    X: np.ndarray
    energy: list = field(default_factory=list)
    step: list = field(default_factory=list)
    local_residual: list = field(default_factory=list)
    rebuilds: int = 0


# ------------------------------------------------ reference: GN on positions


@dataclass
class GNConfig:
    iterations: int = 50
    mu: float = 1e-4
    tol: float = 1e-12


def gn_positions(m, X_init, X_rest, s_goal4, cfg: GNConfig = GNConfig()):
    s3g = cl.sym_goal(s_goal4)
    D0, Me = ops.cr_derivative(m, X_rest)
    tp = topology(m, X_rest)
    nE = tp.nE
    _, area = ops.face_normals(m, X_rest)
    Dx = d_matrix(D0, nE, m.n)
    Mt = sp.diags(np.repeat(Me, 6))
    Lx = (Dx.T @ Mt @ Dx).tocsr()
    W = sp.diags(np.repeat(area, 3))
    X = X_init.copy()
    mu = cfg.mu
    E = goal_energy(tp, tangents(D0, X, nE), s3g, area)
    res = Result(X=X)
    for _ in range(cfg.iterations):
        T = tangents(D0, X, nE)
        r = (S_of(tp, T) - s3g).ravel()
        Jx = jac_S(tp, T) @ Dx
        g = Jx.T @ (W @ r)
        H = (Jx.T @ W @ Jx).tocsc()
        while True:
            dx = spla.spsolve((H + mu * Lx + 1e-12 * sp.identity(3 * m.n)).tocsc(), -g)
            Xn = X + dx.reshape(-1, 3)
            En = goal_energy(tp, tangents(D0, Xn, nE), s3g, area)
            if En < E or mu > 1e6:
                break
            mu *= 10
        mu = max(mu / 3, 1e-10)
        step = float(np.max(np.abs(dx)))
        X, E = Xn, En
        res.energy.append(E)
        res.step.append(step)
        if step < cfg.tol:
            break
    res.X = X
    return res


# ---------------------------------------- the loop: ARAP tangents, NR, Poisson


@dataclass
class TangentConfig:
    iterations: int = 200
    nr_iterations: int = 4
    # "arap_vertex": vertex ARAP rotations pushed to edges by polar
    # ((R_i + R_j) / 2) ("_blend": the plain average) | "arap_face":
    # sum_f R_f D_f X0 | "arap": per-edge R_e G0_e | "dx": D X
    base: str = "arap_vertex"
    eps: float = 1e-8
    damping: float = 1e-12
    # PD schedule: `inner` iterates on a frozen system, then rebuild at the
    # current X. "L": D (read-out and L) at the current metric; "L+G0" also
    # reseeds the ARAP rest there; "never" keeps the rest system throughout.
    inner: int = 30
    rebuild: str = "L+G0"
    rebuild_on_ascent: bool = False
    ascent_tol: float = 1e-3
    tol: float = 1e-11


def project_tangents(tp, T0, s3g, Me, cfg: TangentConfig):
    """`argmin |T - T0|_M s.t. S(T) = s_g`: Newton-Raphson, each step the
    M-weighted min-norm solution of the linearized constraint."""
    Minv = sp.diags(1.0 / np.repeat(Me, 6))
    t = T0.ravel().copy()
    for _ in range(cfg.nr_iterations):
        T = t.reshape(-1, 6)
        r = (S_of(tp, T) - s3g).ravel()
        J = jac_S(tp, T)
        K = (J @ Minv @ J.T + cfg.damping * sp.identity(J.shape[0])).tocsc()
        t = t - Minv @ (J.T @ spla.spsolve(K, r))
    T = t.reshape(-1, 6)
    return T, float(np.sqrt(np.mean((S_of(tp, T) - s3g) ** 2)))


def tangent_loop(m, X_init, X_rest, s_goal4, cfg: TangentConfig = TangentConfig()):
    s3g = cl.sym_goal(s_goal4)
    _, area = ops.face_normals(m, X_rest)

    def operators(X_metric, X_tangent_rest):
        D0, Me, M2, solve = cl._factor(m, X_metric, cfg.eps)
        tp = topology(m, X_metric)
        return D0, Me, M2, solve, tp, tangents(D0, X_tangent_rest, tp.nE)

    D0, Me, M2, solve, tp, G0 = operators(X_rest, X_rest)
    X_arap_rest = X_rest
    er = edge_rot_rest(m, X_arap_rest, D0)
    nE = tp.nE
    X = X_init.copy()
    E = goal_energy(tp, tangents(D0, X, nE), s3g, area)
    res = Result(X=X)

    def rebuild(X):
        nonlocal D0, Me, M2, solve, tp, G0, X_arap_rest, er
        D0, Me, M2, solve, tp, G0 = operators(
            X, X if cfg.rebuild == "L+G0" else X_rest)
        if cfg.rebuild == "L+G0":
            X_arap_rest = X
        er = edge_rot_rest(m, X_arap_rest, D0)
        res.rebuilds += 1

    for it in range(cfg.iterations):
        if it > 0 and cfg.rebuild != "never" and it % cfg.inner == 0:
            rebuild(X)
            E = goal_energy(tp, tangents(D0, X, nE), s3g, area)
        G = tangents(D0, X, nE)
        if cfg.base == "arap_face":
            T0 = face_arap_tangents(m, X_arap_rest, X, Me)
        elif cfg.base in ("arap_vertex", "arap_vertex_blend"):
            Rv, _ = vertex_rotations(m, X_arap_rest, X,
                                     ops.cotan_weights(m, X_arap_rest))
            T0 = vertex_arap_tangents(
                m, Rv, G0, "blend" if cfg.base.endswith("blend") else "polar")
        elif cfg.base == "arap_edge":
            T0 = rotate_tangents(cr_arap_rotations(er, X), G0)
        elif cfg.base == "arap_dir":
            T0 = rotate_tangents(dir_rotations(er, G, spokes=False), G0)
        elif cfg.base == "arap":
            T0 = arap_tangents(G0, G)
        else:
            T0 = G
        Ts, lres = project_tangents(tp, T0, s3g, Me, cfg)
        F = np.vstack([Ts[:, :3], Ts[:, 3:]])
        rhs = D0.T @ (M2 @ F) + cfg.eps * X
        Xn = np.column_stack([solve(rhs[:, k]) for k in range(3)])
        En = goal_energy(tp, tangents(D0, Xn, nE), s3g, area)
        if (cfg.rebuild_on_ascent and cfg.rebuild != "never"
                and En > E * (1 + cfg.ascent_tol)):
            rebuild(X)
            En = goal_energy(tp, tangents(D0, Xn, nE), s3g, area)
        step = float(np.max(np.linalg.norm(Xn - X, axis=1)))
        X, E = Xn, En
        res.energy.append(E)
        res.step.append(step)
        res.local_residual.append(lres)
        if step < cfg.tol:
            break
    res.X = X
    return res


# ------------------------- GN on the goal with ARAP tangents, frozen Hessian


@dataclass
class GoalGNConfig:
    iterations: int = 150
    inner: int = 30  # iterates per outer step (frozen system)
    w_arap: float = 1e-2
    # At each outer step: "none" keeps the rest system; "D" re-derives the
    # CR derivative (read-out, L, mass, areas) at the current metric; "D+G0"
    # also moves the ARAP rest to the current X.
    rebuild: str = "D+G0"
    # "vertex_tangent": w/2 |D X - T0|^2_M with vertex ARAP rotations pushed
    # to edges (`push`); "vertex_arap": the vertex ARAP energy itself (exact
    # gradient, Hessian L_cot); "face_tangent": per-face T0.
    # Edge rotations: "edge_arap": the CR-ARAP energy (diamond Procrustes,
    # exact gradient, Hessian L_cot); "edge_tangent": CR-ARAP R_e on G0;
    # "dir_tangent": vector-Dirichlet-weighted frame Procrustes R_e on G0;
    # "dir_arap": ARAP on frame spokes G_e' - G_e with vector Dirichlet
    # weights (exact gradient, Hessian D^T L_w D); "dir_intrinsic": Stein's
    # L_dir on the intrinsic edge frames U(D X) - U0 (GN Hessian J^T K J).
    memory: str = "vertex_tangent"
    push: str = "polar"  # "polar" | "blend"
    # "frozen": S read through the rest D (metric fixed until a rebuild);
    # "live": D recomputed from the current positions in S and its Jacobian.
    readout: str = "frozen"
    # Where the target tensor lives: "face" (per face, face frame);
    # "edge_pair" (each face's reading in each of its edges' frames against
    # the edge's target, 9 rows per face); "edge_avg" (the mean of an edge's
    # two faces' readings against its target, 3 rows per edge).
    site: str = "face"
    # "gn": H = (J D)^T W (J D) + w L; "memory": H = w L + mu_mass M, the
    # curvature gradient only on the right-hand side.
    lhs: str = "gn"
    mu_mass: float = 0.0
    backtracks: int = 12
    eps: float = 1e-10
    tol: float = 1e-9  # inner convergence (max vertex step), triggers the outer step


@dataclass
class _GoalSystem:
    D0: sp.csr_matrix
    Me: np.ndarray
    tp: Topology
    area: np.ndarray
    Dx: sp.csr_matrix
    Mt: sp.dia_matrix
    Lx: sp.csr_matrix
    W: sp.dia_matrix


def _goal_system(m, X) -> _GoalSystem:
    D0, Me = ops.cr_derivative(m, X)
    tp = topology(m, X)
    _, area = ops.face_normals(m, X)
    Dx = d_matrix(D0, tp.nE, m.n)
    Mt = sp.diags(np.repeat(Me, 6))
    return _GoalSystem(D0=D0, Me=Me, tp=tp, area=area, Dx=Dx, Mt=Mt,
                       Lx=(Dx.T @ Mt @ Dx).tocsr(),
                       W=sp.diags(np.repeat(area, 3)))


def goal_gn_loop(m, X_init, X_rest, s_goal4, cfg: GoalGNConfig = GoalGNConfig()):
    """E(X) = 1/2 |S(D X) - s_g|^2_W + memory (`cfg.memory`). Outer step every `inner` iterates: rebuild the system
    (`cfg.rebuild`) and factor `H = (J D)^T W (J D) + w L`. Inner iterates:
    fresh gradient, the frozen factor, Armijo backtracking; when no step
    decreases E, refactor at the current X and retry."""
    # Face site: `s_goal4` (m, 4) per face. Edge sites: (nE, 3) per edge, in
    # each edge's unfolded frame.
    s3g = cl.sym_goal(s_goal4) if cfg.site == "face" else np.asarray(s_goal4)
    reg = cfg.eps * sp.identity(3 * m.n)
    w = cfg.w_arap
    sy = _goal_system(m, X_rest)
    X_arap = X_rest

    def memory_rest():
        wc = ops.cotan_weights(m, X_arap)
        Lc = sp.kron(ops.cotan_laplacian(m, X_arap), sp.identity(3)).tocsr()
        er = edge_rot_rest(m, X_arap, sy.D0)
        return wc, Lc, tangents(sy.D0, X_arap, sy.tp.nE), er

    wc, Lc, G0, er = memory_rest()
    di = dir_intrinsic_rest(m, X_arap, G0, sy.Me) if cfg.memory == "dir_intrinsic" else None

    def memory_hessian():
        """Memory Hessian, scaled so its trace matches `D^T M D` (the
        tangent memories), so `w` means the same for every memory."""
        if cfg.memory in ("vertex_arap", "edge_arap"):
            H = Lc
        elif cfg.memory == "dir_arap":
            H = (sy.Dx.T @ er.Lw @ sy.Dx).tocsr()
        elif cfg.memory == "dir_intrinsic":
            Jc = dir_intrinsic_jacobian(tangents(sy.D0, X_arap, sy.tp.nE)) @ sy.Dx
            H = (Jc.T @ di.K @ Jc).tocsr()
        else:
            return sy.Lx, 1.0
        s = sy.Lx.diagonal().sum() / H.diagonal().sum()
        return s * H, s

    hm = memory_hessian()

    def memory(X, T):
        """Memory energy, its gradient (3 n), for `parts`."""
        s = hm[1]
        if cfg.memory == "vertex_arap":
            Rv, b = vertex_rotations(m, X_arap, X, wc)
            E = vertex_arap_energy(m, X_arap, X, Rv, wc)
            return w * s * E, w * s * (Lc @ X.ravel() - b.ravel())
        if cfg.memory == "edge_arap":
            E, g = cr_arap_energy_grad(er, X, cr_arap_rotations(er, X))
            return w * s * E, w * s * g
        if cfg.memory == "dir_arap":
            E, gt = dir_arap_energy_grad_T(er, T, dir_rotations(er, T, spokes=True))
            return w * s * E, w * s * (sy.Dx.T @ gt)
        if cfg.memory == "dir_intrinsic":
            E, g = dir_intrinsic_energy_grad(di, T, sy.Dx)
            return w * s * E, w * s * g
        if cfg.memory == "vertex_tangent":
            Rv, _ = vertex_rotations(m, X_arap, X, wc)
            T0 = vertex_arap_tangents(m, Rv, G0, cfg.push)
        elif cfg.memory == "edge_tangent":
            T0 = rotate_tangents(cr_arap_rotations(er, X), G0)
        elif cfg.memory == "dir_tangent":
            T0 = rotate_tangents(dir_rotations(er, T, spokes=False), G0)
        else:
            T0 = face_arap_tangents(m, X_arap, X, sy.Me)
        dT = (T - T0).ravel()
        return (0.5 * w * dT @ (sy.Mt @ dT),
                w * (sy.Dx.T @ (sy.Mt @ dT)))

    import face_loops as fl
    import shape_dirichlet as sd
    st = fl.build_stencils(m, sy.D0)
    nF, nE = len(m.F), len(m.E)
    edge_cnt = np.bincount(m.fe.ravel(), minlength=nE).astype(float)

    def Wm():
        """Read-out weights: face area (face site), area / 3 per face-edge
        reading ("edge_pair"), edge mass A_e / 3 ("edge_avg")."""
        if cfg.site == "edge_pair":
            return sp.diags(np.repeat(sy.area / 3.0, 9))
        if cfg.site == "edge_avg":
            return sp.diags(np.repeat(sy.Me, 3))
        return sy.W

    def readout_residual(X, T):
        if cfg.site == "edge_pair":
            R = sd.face_edge_readings(st, X[st.ids])
            return (R - s3g[m.fe]).ravel()
        if cfg.site == "edge_avg":
            return (sd.edge_average_readings(m, st, X) - s3g).ravel()
        if cfg.readout == "live":
            return (fl.stencil_s3_live(st, X[st.ids]) - s3g).ravel()
        return (S_of(sy.tp, T) - s3g).ravel()

    def readout_jacobian(X, T):
        """`dS / dX`. "live": per-face central differences of the live
        read-out over its 6-vertex stencil (D follows the metric). Edge sites
        differentiate the per-face edge readings the same way."""
        if cfg.site == "face" and cfg.readout != "live":
            return jac_S(sy.tp, T) @ sy.Dx
        P = X[st.ids]
        edge_site = cfg.site != "face"
        k_out = 9 if edge_site else 3
        fn = ((lambda Q: sd.face_edge_readings(st, Q).reshape(nF, 9)) if edge_site
              else (lambda Q: fl.stencil_s3_live(st, Q)))
        blocks = np.zeros((nF, k_out, 18))
        h = 1e-6
        for d in range(18):
            k, ax = divmod(d, 3)
            Pp, Pm = P.copy(), P.copy()
            Pp[:, k, ax] += h
            Pm[:, k, ax] -= h
            blocks[:, :, d] = (fn(Pp) - fn(Pm)) / (2 * h)
        blocks *= np.repeat(st.mask, 3, axis=1)[:, None, :]
        cols = (3 * st.ids[:, :, None] + np.arange(3)[None, None, :]).reshape(nF, 18)
        cols = np.repeat(cols[:, None, :], k_out, 1)
        if cfg.site == "edge_avg":
            rows = (3 * m.fe[:, :, None] + np.arange(3)[None, None, :]).reshape(nF, 9)
            vals = blocks / np.repeat(edge_cnt[m.fe], 3, 1)[:, :, None]
            nrow = 3 * nE
        else:
            rows = k_out * np.arange(nF)[:, None] + np.arange(k_out)[None, :]
            vals = blocks
            nrow = k_out * nF
        rows = np.repeat(rows[:, :, None], 18, 2)
        return sp.csr_matrix((vals.ravel(), (rows.ravel(), cols.ravel())),
                             shape=(nrow, 3 * m.n))

    def parts(X):
        T = tangents(sy.D0, X, sy.tp.nE)
        r = readout_residual(X, T)
        Em, gm = memory(X, T)
        return T, r, gm, 0.5 * r @ (Wm() @ r) + Em

    def factor(X):
        Hm = hm[0]
        if cfg.lhs == "memory":
            # Curvature only on the right-hand side: H = w L + mu M (PD-shaped).
            Mv = sp.diags(np.repeat(ops.vertex_areas(m, X_arap), 3))
            return spla.factorized((w * Hm + cfg.mu_mass * Mv + reg).tocsc())
        Jx = readout_jacobian(X, tangents(sy.D0, X, sy.tp.nE))
        return spla.factorized((Jx.T @ Wm() @ Jx + w * Hm + reg).tocsc())

    X = X_init.copy()
    res = Result(X=X)
    res.outer = 0
    res.refactors = 0
    solve = factor(X)
    T, r, gm, E = parts(X)
    since_outer, force_outer, progress = 0, False, True
    for it in range(cfg.iterations):
        if since_outer >= cfg.inner or force_outer:
            if force_outer and not progress:
                break  # converged right after a rebuild: nothing left to update
            since_outer, force_outer, progress = 0, False, False
            if cfg.rebuild != "none":
                sy = _goal_system(m, X)
            if cfg.rebuild == "D+G0":
                X_arap = X
            wc, Lc, G0, er = memory_rest()
            if cfg.memory == "dir_intrinsic":
                di = dir_intrinsic_rest(m, X_arap, G0, sy.Me)
            hm = memory_hessian()
            solve = factor(X)
            T, r, gm, E = parts(X)
            res.outer += 1
        Jx = readout_jacobian(X, T)
        g = Jx.T @ (Wm() @ r) + gm
        accepted = False
        for attempt in range(2):
            dx = -solve(g)
            slope = float(g @ dx)
            a = 1.0
            for _ in range(cfg.backtracks):
                Xn = X + a * dx.reshape(-1, 3)
                Tn, rn, gmn, En = parts(Xn)
                if slope < 0 and En <= E + 1e-4 * a * slope:
                    accepted = True
                    break
                a *= 0.5
            if accepted or attempt == 1:
                break
            solve = factor(X)
            res.refactors += 1
        since_outer += 1
        if accepted:
            step = float(np.max(np.abs(a * dx)))
            progress = progress or step > cfg.tol
            X, T, r, gm, E = Xn, Tn, rn, gmn, En
        else:
            step = 0.0
        res.energy.append(float(0.5 * r @ (Wm() @ r)))
        res.step.append(step)
        if step < cfg.tol:
            # Inner problem converged: the outer step updates the system.
            if cfg.rebuild == "none":
                break
            force_outer = True
    res.X = X
    return res


if __name__ == "__main__":
    import time
    import base_loop as bl
    import face_loops as fl

    m, Xf, Xc, s4 = fl.strip_case()
    D0, Me = ops.cr_derivative(m, Xf)
    tp = topology(m, Xf)
    s_ref = cl.sym_goal(bl.readout(m, D0, Xc))
    print("read-out from tangents vs embedding",
          np.abs(S_of(tp, tangents(D0, Xc, tp.nE)) - s_ref).max())
    T = tangents(D0, Xc, tp.nE)
    J = jac_S(tp, T)
    dT = np.random.default_rng(0).normal(size=T.shape) * 1e-6
    lin = (S_of(tp, T + dT) - S_of(tp, T)).ravel()
    print("jac_S check", np.abs(lin - J @ dT.ravel()).max() / np.abs(lin).max())

    for name, fn in (("gn", gn_positions), ("tangent", tangent_loop)):
        t = time.time()
        r = fn(m, Xf, Xf, s4)
        print(f"{name:8s} its {len(r.step):3d} err {bl.rigid_error(r.X, Xc):.3e} "
              f"E {r.energy[0]:.2e} -> {r.energy[-1]:.2e} t {time.time() - t:.1f}s")
