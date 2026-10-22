"""Phase B loop from MATH.md: X = L^{-1} D^T M F(J^{-1}(S_g)), iterated.

Each iteration:
1. G = D X (edge frames).
2. Continual ARAP per edge: R_e = polar(G_e G0_e^T), T_e = R_e G0_e.
3. Local inverse of J per face: c_f* = K3_f^{-1} s_g,f (perpendicular normal
   components), correction along N_f.
4. Recombine the two faces' corrections per edge (local average, or the
   exact global least-norm solve).
5. X = L^{-1} (D^T M F + eps X), L = D^T M D + eps I, prefactored.
Frozen-Hessian policy: rebuild L (or L and G0) at the current X when E rises.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_ops as ops


@dataclass
class CRLoopConfig:
    iterations: int = 200
    tol: float = 1e-11
    eps: float = 1e-8
    recombine: str = "local"  # "local" | "global"
    lam_global: float = 1e4
    rebuild: str = "on_ascent"  # "on_ascent" | "never"
    rebuild_level: str = "L"  # "L" | "L+G0"


@dataclass
class CRLoopResult:
    X: np.ndarray
    energy: list = field(default_factory=list)
    curvature: list = field(default_factory=list)
    step: list = field(default_factory=list)
    rebuilds: int = 0


def sym_goal(s4: np.ndarray) -> np.ndarray:
    """`(m, 4)` coefficients to `(m, 3)` symmetric `(k11, k12, k22)`."""
    return np.stack([s4[:, 0], 0.5 * (s4[:, 1] + s4[:, 2]), s4[:, 3]], axis=1)


def face_perp_maps(m: ops.Mesh, X: np.ndarray):
    """Per face `K3_f` (3 x 3): perpendicular normal components of its three
    edges (canonical perp column) to `(k11, k12, k22)` in the face frame."""
    N, area = ops.face_normals(m, X)
    T, B, _ = ops.face_frames(m, X)
    K3 = np.zeros((len(m.F), 3, 3))
    for f, t in enumerate(m.F):
        for i in range(3):
            tail, head = t[i], t[(i + 1) % 3]
            vpar = X[head] - X[tail]
            li = np.linalg.norm(vpar)
            vpar /= li
            si = 1.0 if tail < head else -1.0
            vperp = np.cross(N[f], vpar)
            vperp /= np.linalg.norm(vperp)
            Mq = np.outer(vperp, vperp)
            k = ops.K_CR_SHAPE_SIGN * si * li / area[f]
            K3[f, :, i] = k * np.array(
                [T[f] @ Mq @ T[f], T[f] @ Mq @ B[f], B[f] @ Mq @ B[f]])
    return K3, N, area


def arap_frames(G0p, G0q, Gp, Gq):
    """Per edge Procrustes of the current frame onto the rest frame."""
    nE = len(G0p)
    Tp = np.zeros_like(G0p)
    Tq = np.zeros_like(G0q)
    for e in range(nE):
        G0 = np.column_stack([G0p[e], G0q[e]])
        G = np.column_stack([Gp[e], Gq[e]])
        U, _, Vt = np.linalg.svd(G @ G0.T)
        R = U @ Vt
        if np.linalg.det(R) < 0:
            U[:, -1] *= -1
            R = U @ Vt
        Tp[e] = R @ G0p[e]
        Tq[e] = R @ G0q[e]
    return Tp, Tq


def target_frames(m: ops.Mesh, X: np.ndarray, D0, G0p, G0q, s3_goal,
                  cfg: CRLoopConfig, Me: np.ndarray):
    Gp, Gq = ops.edge_G(m, D0, X)
    Tp, Tq = arap_frames(G0p, G0q, Gp, Gq)
    K3, N, area = face_perp_maps(m, X)
    nE = len(m.E)
    Fq = Tq.copy()
    if cfg.recombine == "local":
        acc = np.zeros((nE, 3))
        cnt = np.array([len(fs) for fs in m.edge_faces], dtype=float)
        for f in range(len(m.F)):
            c = np.array([N[f] @ Tq[m.fe[f, i]] for i in range(3)])
            c_star = np.linalg.solve(K3[f], s3_goal[f])
            for i in range(3):
                acc[m.fe[f, i]] += N[f] * (c_star[i] - c[i])
        Fq = Tq + acc / cnt[:, None]
    else:
        # Exact: min sum Me |d_e|^2 + lam sum A_f |K3_f (N_f . (Tq + d)) - s_f|^2.
        rows, cols, vals = [], [], []
        rhs = np.zeros(3 * len(m.F))
        for f in range(len(m.F)):
            c = np.array([N[f] @ Tq[m.fe[f, i]] for i in range(3)])
            rhs[3 * f:3 * f + 3] = s3_goal[f] - K3[f] @ c
            for r in range(3):
                for i in range(3):
                    e = m.fe[f, i]
                    for ax in range(3):
                        rows.append(3 * f + r)
                        cols.append(3 * e + ax)
                        vals.append(K3[f, r, i] * N[f, ax])
        Bm = sp.csr_matrix((vals, (rows, cols)), shape=(3 * len(m.F), 3 * nE))
        W = sp.diags(np.repeat(area, 3))
        H = sp.diags(np.repeat(Me, 3)) + cfg.lam_global * (Bm.T @ W @ Bm)
        d = spla.spsolve(H.tocsc(), cfg.lam_global * (Bm.T @ (W @ rhs)))
        Fq = Tq + d.reshape(nE, 3)
    return Tp, Fq, Gp, Gq


def energy(Me, Gp, Gq, Fp, Fq) -> float:
    return 0.5 * float(np.sum(Me * (np.sum((Gp - Fp) ** 2, 1)
                                     + np.sum((Gq - Fq) ** 2, 1))))


def curvature_residual(m, D0, X, s3_goal) -> float:
    S = ops.face_shape(m, X, D0)
    T, B, _ = ops.face_frames(m, X)
    s = sym_goal(ops.face_coeffs(S, T, B))
    _, area = ops.face_normals(m, X)
    return float(np.sqrt(np.sum(area[:, None] * (s - s3_goal) ** 2)))


def _factor(m, X_rest, eps):
    D0, Me = ops.cr_derivative(m, X_rest)
    M2 = sp.diags(np.concatenate([Me, Me]))
    L = (D0.T @ M2 @ D0 + eps * sp.identity(m.n)).tocsc()
    return D0, Me, M2, spla.factorized(L)


def run(m: ops.Mesh, X_init: np.ndarray, X_rest: np.ndarray,
        s_goal4: np.ndarray, cfg: CRLoopConfig = CRLoopConfig()):
    s3 = sym_goal(s_goal4)
    D0, Me, M2, solve = _factor(m, X_rest, cfg.eps)
    G0p, G0q = ops.edge_G(m, D0, X_rest)
    X = X_init.copy()
    res = CRLoopResult(X=X)
    Tp, Fq, Gp, Gq = target_frames(m, X, D0, G0p, G0q, s3, cfg, Me)
    E = energy(Me, Gp, Gq, Tp, Fq)
    just_rebuilt = False
    for _ in range(cfg.iterations):
        F = np.vstack([Tp, Fq])
        rhs = D0.T @ (M2 @ F) + cfg.eps * X
        X_new = np.column_stack([solve(rhs[:, k]) for k in range(3)])
        step = float(np.max(np.linalg.norm(X_new - X, axis=1)))
        X = X_new
        Tp, Fq, Gp, Gq = target_frames(m, X, D0, G0p, G0q, s3, cfg, Me)
        E_new = energy(Me, Gp, Gq, Tp, Fq)
        if (cfg.rebuild == "on_ascent" and E_new > E * (1 + 1e-9)
                and not just_rebuilt):
            D0, Me, M2, solve = _factor(m, X, cfg.eps)
            if cfg.rebuild_level == "L+G0":
                G0p, G0q = ops.edge_G(m, D0, X)
            Tp, Fq, Gp, Gq = target_frames(m, X, D0, G0p, G0q, s3, cfg, Me)
            E_new = energy(Me, Gp, Gq, Tp, Fq)
            res.rebuilds += 1
            just_rebuilt = True
        else:
            just_rebuilt = False
        E = E_new
        res.energy.append(E)
        res.curvature.append(curvature_residual(m, D0, X, s3))
        res.step.append(step)
        if step < cfg.tol:
            break
    res.X = X
    return res
