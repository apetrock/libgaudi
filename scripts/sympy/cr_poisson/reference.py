"""Ground truth for phase B: `least_squares` on positions directly.

Residuals: area-weighted shape coefficients against the goal (live read-out,
same as the loop), edge lengths against the rest (isometry), and a tiny anchor
that fixes the rigid motion.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.optimize import least_squares

import cr_ops as ops
from base_loop import readout


def _sparsity(m: ops.Mesh):
    nF, n, nE = len(m.F), m.n, len(m.E)
    face_verts = []
    for f in range(nF):
        vs = set(m.F[f])
        for e in m.fe[f]:
            for g in m.edge_faces[e]:
                vs.update(m.F[g])
        face_verts.append(sorted(vs))
    rows, cols = [], []
    for f, vs in enumerate(face_verts):
        for r in range(4):
            for v in vs:
                for ax in range(3):
                    rows.append(4 * f + r)
                    cols.append(3 * v + ax)
    base = 4 * nF
    for e, (i, j) in enumerate(m.E):
        for v in (i, j):
            for ax in range(3):
                rows.append(base + e)
                cols.append(3 * v + ax)
    base += nE
    for v in range(n):
        for ax in range(3):
            rows.append(base + 3 * v + ax)
            cols.append(3 * v + ax)
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)),
                         shape=(4 * nF + nE + 3 * n, 3 * n))


def solve(m: ops.Mesh, X_init: np.ndarray, X_rest: np.ndarray,
          s_goal: np.ndarray, lam: float = 1.0, w_iso: float = 1.0,
          anchor: float = 1e-4) -> np.ndarray:
    D0, _ = ops.cr_derivative(m, X_rest)
    l0 = np.linalg.norm(X_rest[m.E[:, 1]] - X_rest[m.E[:, 0]], axis=1)
    _, area0 = ops.face_normals(m, X_rest)
    sa = np.sqrt(lam * np.repeat(area0, 4))

    def fun(z):
        X = z.reshape(-1, 3)
        s = readout(m, D0, X)
        rs = sa * (s - s_goal).ravel()
        l = np.linalg.norm(X[m.E[:, 1]] - X[m.E[:, 0]], axis=1)
        ri = np.sqrt(w_iso) * (l - l0)
        ra = np.sqrt(anchor) * (X - X_init).ravel()
        return np.concatenate([rs, ri, ra])

    out = least_squares(fun, X_init.ravel(), jac_sparsity=_sparsity(m),
                        method="trf", x_scale="jac", max_nfev=400,
                        xtol=1e-12, ftol=1e-12)
    return out.x.reshape(-1, 3)
