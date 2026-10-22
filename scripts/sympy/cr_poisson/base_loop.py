"""Phase B base loop: shape-operator targets through the CR derivative.

One outer iteration:
1. Read the edge derivatives `G = D0 X` of the current positions (two
   vectors per edge; `D0` is intrinsic, so the rest operator stays exact under
   isometric bending).
2. G-step (the inverse map g): the shape read-out is linear in `G` once the
   live face normals and frames are frozen, `s = B(X) g`. Solve
   `min |dg|^2_M + lam |B (g + dg) - s_goal|^2_A` for the target `G* = g + dg`.
3. R-step: per-vertex Procrustes on cotan-weighted spokes against the ARAP rest.
4. X-step (prefactored, constant): `(w_a L_cot + w_d D0^T M D0 + eps Mv) X =
   w_a b(R) + w_d D0^T M G* + eps Mv X`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_ops as ops


@dataclass
class LoopConfig:
    w_arap: float = 1.0
    w_cr: float = 1.0
    eps: float = 1e-3
    lam: float = 1.0
    iterations: int = 200
    tol: float = 1e-10
    relinearize_every: int = 0  # 0: keep the rest D0 / L_cr / L_cot


@dataclass
class LoopResult:
    X: np.ndarray
    energy: list = field(default_factory=list)
    curvature_residual: list = field(default_factory=list)
    step: list = field(default_factory=list)
    factorizations: int = 0


def _system(m, X_rest, cfg):
    D0, Me = ops.cr_derivative(m, X_rest)
    M2 = sp.diags(np.concatenate([Me, Me]))
    L_cr = (D0.T @ M2 @ D0).tocsc()
    w = ops.cotan_weights(m, X_rest)
    L_cot = ops.cotan_laplacian(m, X_rest).tocsc()
    Mv = sp.diags(ops.vertex_areas(m, X_rest))
    A = (cfg.w_arap * L_cot + cfg.w_cr * L_cr + cfg.eps * Mv).tocsc()
    return D0, Me, M2, w, Mv, spla.factorized(A)


def readout(m: ops.Mesh, D0, X: np.ndarray) -> np.ndarray:
    """Per-face shape coefficients `(m, 4)` in each face's own frame."""
    S = ops.face_shape(m, X, D0)
    T, B, _ = ops.face_frames(m, X)
    return ops.face_coeffs(S, T, B)


def run(m: ops.Mesh, X_init: np.ndarray, X_arap_rest: np.ndarray,
        s_goal: np.ndarray, cfg: LoopConfig = LoopConfig()) -> LoopResult:
    X = X_init.copy()
    nE = len(m.E)
    D0, Me, M2, w, Mv, solve = _system(m, X_init, cfg)
    res = LoopResult(X=X)
    res.factorizations = 1
    s_goal_flat = s_goal.ravel()
    for it in range(cfg.iterations):
        if cfg.relinearize_every and it > 0 and it % cfg.relinearize_every == 0:
            D0, Me, M2, w, Mv, solve = _system(m, X, cfg)
            res.factorizations += 1
        _, area = ops.face_normals(m, X)
        Gp, Gq = ops.edge_G(m, D0, X)
        g = ops.pack_G(Gp, Gq)
        Bm = ops.shape_readout_rows(m, X)
        r = s_goal_flat - Bm @ g
        Af = sp.diags(np.repeat(area, 4))
        M6 = sp.diags(np.repeat(Me, 6))
        H = (M6 + cfg.lam * (Bm.T @ Af @ Bm)).tocsc()
        dg = spla.spsolve(H, cfg.lam * (Bm.T @ (Af @ r)))
        Gp_s, Gq_s = ops.unpack_G(g + dg, nE)
        G_star = np.vstack([Gp_s, Gq_s])

        _, b = ops.arap_rhs_and_rotations(m, X_arap_rest, X, w)
        rhs = (cfg.w_arap * b + cfg.w_cr * (D0.T @ (M2 @ G_star))
               + cfg.eps * (Mv @ X))
        X_new = np.column_stack([solve(rhs[:, k]) for k in range(3)])

        step = float(np.max(np.linalg.norm(X_new - X, axis=1)))
        X = X_new
        cres = float(np.sqrt(np.sum(area[:, None] * r.reshape(-1, 4) ** 2)))
        res.curvature_residual.append(cres)
        res.step.append(step)
        res.energy.append(float(r @ (Af @ r)))
        if step < cfg.tol:
            break
    res.X = X
    return res


def rigid_error(X: np.ndarray, Y: np.ndarray) -> float:
    """Max distance after the best rigid alignment of `X` onto `Y`."""
    cx, cy = X.mean(0), Y.mean(0)
    H = (X - cx).T @ (Y - cy)
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T
    return float(np.max(np.linalg.norm((X - cx) @ R.T + cy - Y, axis=1)))
