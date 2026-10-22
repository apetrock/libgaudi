"""Phase B, MATH.md section 7: the local inverse of J moves faces.

Option 1, stencil projection (PD with a nonlinear local step):
    local:  p_f = argmin |p - X_f|^2  s.t.  s_f(p) = s_g,f, by Newton-Raphson
            (min-norm steps) on the stencil of face f (its 3 vertices and the
            3 opposite vertices, the support of S_f).
    global: min sum_f w_f |C_f X - C_f p_f|^2 + eps |X - X_k|^2, C_f centers
            the stencil (prefactored, topology only).

Option 2, face rotations (curvature as the connection between faces):
    local:  per face c*_f = K3_f^{-1} s_g,f (exact J^{-1}); per interior edge
            the target dihedral theta*_e = asin(c*_+ / |g_perp|) - asin(c*_- / |g_perp|);
            rotations from min sum_e |R_- - R_+ P_e|^2 + mu sum_f A_f |R_f - Rproc_f|^2,
            P_e = Rot(a0_e, theta*_e - theta0_e), then polar.
    global: per-face ARAP Poisson L_cot X = div(R T0) (prefactored).

Both relations used by option 2 hold exactly on the cylinder (checked):
theta_e = asin(c_+/|g|) - asin(c_-/|g|) and R_- = R_+ Rot(a0_e, theta_e), with
`+` the face whose winding runs along the canonical edge E[e, 0] -> E[e, 1].
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import cr_ops as ops
import cr_loop as cl


# ------------------------------------------------------------- stencils


@dataclass
class Stencils:
    ids: np.ndarray  # (m, 6) vertex ids, padded with the face's first vertex
    mask: np.ndarray  # (m, 6) 1 for real stencil vertices
    W: np.ndarray  # (m, 3, 6) perpendicular CR derivative rows per local edge
    si: np.ndarray  # (m, 3) orientation signs


def build_stencils(m: ops.Mesh, D0: sp.csr_matrix) -> Stencils:
    nF, nE = len(m.F), len(m.E)
    ids = np.repeat(m.F[:, :1], 6, axis=1)
    ids[:, :3] = m.F
    mask = np.zeros((nF, 6))
    mask[:, :3] = 1
    W = np.zeros((nF, 3, 6))
    si = np.zeros((nF, 3))
    D0 = D0.tocsr()
    for f, t in enumerate(m.F):
        for i in range(3):
            e = m.fe[f, i]
            si[f, i] = 1.0 if t[i] < t[(i + 1) % 3] else -1.0
            for g in m.edge_faces[e]:
                if g == f:
                    continue
                o = [v for v in m.F[g] if v not in m.E[e]][0]
                ids[f, 3 + i] = o
                mask[f, 3 + i] = 1
        local = {int(v): k for k, v in enumerate(ids[f]) if mask[f, k]}
        for i in range(3):
            row = D0.getrow(nE + m.fe[f, i])
            for v, x in zip(row.indices, row.data):
                W[f, i, local[int(v)]] += x
    return Stencils(ids=ids, mask=mask, W=W, si=si)


def stencil_W_live(st: Stencils, P: np.ndarray) -> np.ndarray:
    """Perpendicular CR derivative rows (m, 3, 6) recomputed from the stencil
    positions (same assembly as `ops.cr_derivative`, M^{-1}-scaled), so the
    read-out follows the current metric. Closed meshes: every edge has its
    opposite vertex in slot `3 + i`."""
    def sq(a, b):
        return np.sum((a - b) ** 2, 1)
    a, b, c = P[:, 0], P[:, 1], P[:, 2]
    Af = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    W = np.zeros((len(P), 3, 6))
    for i in range(3):
        t, h, o = i, (i + 1) % 3, (i + 2) % 3
        Pt, Ph, Po, Po2 = P[:, t], P[:, h], P[:, o], P[:, 3 + i]
        si = st.si[:, i]
        l = np.sqrt(sq(Ph, Pt))
        ejk, eki = sq(Po, Ph), sq(Pt, Po)
        k = -si / (12 * l)
        W[:, i, t] += k * (l * l - (eki - ejk))
        W[:, i, h] += k * (l * l + (eki - ejk))
        W[:, i, o] += si * l / 6
        has = st.mask[:, 3 + i] > 0
        Ag = 0.5 * np.linalg.norm(np.cross(Pt - Ph, Po2 - Ph), axis=1)
        ejk2, eki2 = sq(Po2, Pt), sq(Ph, Po2)  # other face: tail' = h, head' = t
        k2 = si / (12 * l)
        W[:, i, h] += has * k2 * (l * l - (eki2 - ejk2))
        W[:, i, t] += has * k2 * (l * l + (eki2 - ejk2))
        W[:, i, 3 + i] += has * (-si * l / 6)
        Me = (Af + has * Ag) / 3
        W[:, i] /= Me[:, None]
    return W


def stencil_s3_live(st: Stencils, P: np.ndarray) -> np.ndarray:
    """The paper's read-out at the current metric (D recomputed)."""
    live = Stencils(ids=st.ids, mask=st.mask, W=stencil_W_live(st, P), si=st.si)
    return stencil_s3(live, P)


def stencil_s3(st: Stencils, P: np.ndarray) -> np.ndarray:
    """`(k11, k12, k22)` per face from stencil positions `P` (m, 6, 3); equals
    `sym_goal(face_coeffs(face_shape(X)))` with the rest `D0`."""
    a, b, c = P[:, 0], P[:, 1], P[:, 2]
    n = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(n, axis=1)
    N = n / (2 * area[:, None])
    T = (b - a) / np.linalg.norm(b - a, axis=1)[:, None]
    B = np.cross(N, T)
    s = np.zeros((len(P), 3))
    for i in range(3):
        vpar = P[:, (i + 1) % 3] - P[:, i]
        li = np.linalg.norm(vpar, axis=1)
        vperp = np.cross(N, vpar / li[:, None])
        gq = np.einsum("fk,fkd->fd", st.W[:, i], P)
        k = ops.K_CR_SHAPE_SIGN * st.si[:, i] * li / area * np.sum(N * gq, 1)
        tv = np.sum(T * vperp, 1)
        bv = np.sum(B * vperp, 1)
        s += k[:, None] * np.stack([tv * tv, tv * bv, bv * bv], 1)
    return s


# ------------------------------------------------- option 1: stencil projection


@dataclass
class StencilConfig:
    iterations: int = 200
    nr_iterations: int = 5
    eps: float = 1e-3
    rcond: float = 1e-3  # boundary faces: N_f . g_perp = 0 on the boundary edge, J has rank 2
    tol: float = 1e-10


@dataclass
class LoopResult:
    X: np.ndarray
    curvature: list = field(default_factory=list)
    energy: list = field(default_factory=list)
    step: list = field(default_factory=list)
    local_residual: list = field(default_factory=list)


def project_stencils(st: Stencils, P: np.ndarray, s3g: np.ndarray,
                     cfg: StencilConfig) -> tuple[np.ndarray, float]:
    """Newton-Raphson with min-norm (rank-truncated pseudo-inverse) steps on
    each stencil, all faces at once."""
    P = P.copy()
    h = 1e-7
    cols = np.repeat(st.mask, 3, axis=1)  # (m, 18)
    for _ in range(cfg.nr_iterations):
        C = stencil_s3(st, P) - s3g
        J = np.zeros((len(P), 3, 18))
        for d in range(18):
            k, ax = divmod(d, 3)
            Pp = P.copy()
            Pp[:, k, ax] += h
            J[:, :, d] = (stencil_s3(st, Pp) - s3g - C) / h
        J *= cols[:, None, :]
        dp = -(np.linalg.pinv(J, rcond=cfg.rcond) @ C[:, :, None])[:, :, 0]
        P += dp.reshape(-1, 6, 3)
    res = float(np.sqrt(np.mean(np.sum((stencil_s3(st, P) - s3g) ** 2, 1))))
    return P, res


def _stencil_system(m: ops.Mesh, st: Stencils, w: np.ndarray, eps: float):
    rows, cols, vals = [], [], []
    for f in range(len(m.F)):
        idx = st.ids[f][st.mask[f] > 0]
        k = len(idx)
        Cm = np.eye(k) - 1.0 / k
        Q = w[f] * Cm.T @ Cm
        for a in range(k):
            for b in range(k):
                rows.append(idx[a])
                cols.append(idx[b])
                vals.append(Q[a, b])
    A = sp.csr_matrix((vals, (rows, cols)), shape=(m.n, m.n))
    return spla.factorized((A + eps * sp.identity(m.n)).tocsc())


def run_stencil(m: ops.Mesh, X_init: np.ndarray, X_rest: np.ndarray,
                s_goal4: np.ndarray, cfg: StencilConfig = StencilConfig()):
    s3g = cl.sym_goal(s_goal4)
    D0, _ = ops.cr_derivative(m, X_rest)
    st = build_stencils(m, D0)
    _, area = ops.face_normals(m, X_rest)
    w = area / area.mean()
    solve = _stencil_system(m, st, w, cfg.eps)
    X = X_init.copy()
    res = LoopResult(X=X)
    for _ in range(cfg.iterations):
        P, lres = project_stencils(st, X[st.ids], s3g, cfg)
        rhs = cfg.eps * X
        E = 0.0
        for f in range(len(m.F)):
            sel = st.mask[f] > 0
            idx = st.ids[f][sel]
            p = P[f][sel]
            q = X[idx]
            cp = p - p.mean(0)
            E += 0.5 * w[f] * np.sum(((q - q.mean(0)) - cp) ** 2)
            np.add.at(rhs, idx, w[f] * (cp - cp.mean(0)))
        X_new = np.column_stack([solve(rhs[:, k]) for k in range(3)])
        step = float(np.max(np.linalg.norm(X_new - X, axis=1)))
        X = X_new
        res.curvature.append(cl.curvature_residual(m, D0, X, s3g))
        res.energy.append(E)
        res.step.append(step)
        res.local_residual.append(lres)
        if step < cfg.tol:
            break
    res.X = X
    return res


# ---------------------------------------------------- option 2: face rotations


@dataclass
class RotationConfig:
    iterations: int = 200
    mu: float = 1e-2
    eps: float = 1e-8
    tol: float = 1e-10
    global_step: str = "cot"  # "cot": L_cot X = div(R T0) | "cr": D^T M D X = D^T M sum_f R_f D_f X0


def _rot(a: np.ndarray, t: float) -> np.ndarray:
    a = a / np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K


def _polar(M: np.ndarray) -> np.ndarray:
    U, _, Vt = np.linalg.svd(M)
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    return R


def _face_procrustes(m: ops.Mesh, X0: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Per face rotation of the rest triangle (and its normal) onto the current."""
    N0, _ = ops.face_normals(m, X0)
    N, _ = ops.face_normals(m, X)
    R = np.zeros((len(m.F), 3, 3))
    for f, t in enumerate(m.F):
        A = X0[t] - X0[t].mean(0)
        B = X[t] - X[t].mean(0)
        R[f] = _polar(np.vstack([B, N[f]]).T @ np.vstack([A, N0[f]]))
    return R


def _interior_edges(m: ops.Mesh):
    """`(e, f_plus, i_plus, f_minus, i_minus)` per interior edge."""
    out = []
    for e, fs in enumerate(m.edge_faces):
        if len(fs) < 2:
            continue
        side = {}
        for f in fs:
            i = list(m.fe[f]).index(e)
            t = m.F[f]
            side[1 if t[i] < t[(i + 1) % 3] else -1] = (f, i)
        out.append((e, *side[1], *side[-1]))
    return out


def dihedrals(m: ops.Mesh, X: np.ndarray, edges) -> np.ndarray:
    N, _ = ops.face_normals(m, X)
    th = np.zeros(len(edges))
    for k, (e, fp, _, fm, _) in enumerate(edges):
        a = X[m.E[e, 1]] - X[m.E[e, 0]]
        a /= np.linalg.norm(a)
        th[k] = np.arctan2(np.cross(N[fp], N[fm]) @ a, N[fp] @ N[fm])
    return th


def target_dihedrals(m: ops.Mesh, X: np.ndarray, D0, s3g, edges) -> np.ndarray:
    """Exact per-face inverse of J, then the tilt pair per edge."""
    K3, _, _ = cl.face_perp_maps(m, X)
    cstar = np.linalg.solve(K3, s3g[:, :, None])[:, :, 0]
    _, Gq = ops.edge_G(m, D0, X)
    th = np.zeros(len(edges))
    for k, (e, fp, ip, fm, im) in enumerate(edges):
        g = np.linalg.norm(Gq[e])
        th[k] = (np.arcsin(np.clip(cstar[fp, ip] / g, -1, 1))
                 - np.arcsin(np.clip(cstar[fm, im] / g, -1, 1)))
    return th


def sync_rotations(m: ops.Mesh, X0: np.ndarray, edges, dtheta, Rproc, wA, mu):
    """min sum_e |R_- - R_+ P_e|^2 + mu sum_f wA_f |R_f - Rproc_f|^2, per row."""
    nF = len(m.F)
    rows, cols, vals = [], [], []
    for k, (e, fp, _, fm, _) in enumerate(edges):
        a0 = X0[m.E[e, 1]] - X0[m.E[e, 0]]
        Pt = _rot(a0, dtheta[k]).T
        for r in range(3):
            for c in range(3):
                rows.append(3 * k + r)
                cols.append(3 * fp + c)
                vals.append(-Pt[r, c])
            rows.append(3 * k + r)
            cols.append(3 * fm + r)
            vals.append(1.0)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(3 * len(edges), 3 * nF))
    reg = sp.diags(np.repeat(mu * wA, 3))
    H = (A.T @ A + reg).tocsc()
    solve = spla.factorized(H)
    R = np.zeros((nF, 3, 3))
    for r in range(3):
        rhs = reg @ Rproc[:, r, :].ravel()
        R[:, r, :] = solve(rhs).reshape(nF, 3)
    return np.array([_polar(Rf) for Rf in R])


def _face_arap_system(m: ops.Mesh, X0: np.ndarray, eps: float):
    """Per-face ARAP: sum_f sum_i w_fi |x_h - x_t - R_f d0|^2, w = cot / 2."""
    L = ops.cotan_laplacian(m, X0)
    wfi = np.zeros((len(m.F), 3))
    for f, t in enumerate(m.F):
        for i in range(3):
            a, b, o = t[i], t[(i + 1) % 3], t[(i + 2) % 3]
            u, v = X0[a] - X0[o], X0[b] - X0[o]
            wfi[f, i] = 0.5 * np.dot(u, v) / np.linalg.norm(np.cross(u, v))
    return wfi, spla.factorized((L + eps * sp.identity(m.n)).tocsc())


def _face_arap_rhs(m, X0, R, wfi):
    b = np.zeros((m.n, 3))
    for f, t in enumerate(m.F):
        for i in range(3):
            tail, head = t[i], t[(i + 1) % 3]
            r = wfi[f, i] * R[f] @ (X0[head] - X0[tail])
            b[head] += r
            b[tail] -= r
    return b


def _cr_rotated_frames(m: ops.Mesh, X0: np.ndarray, R: np.ndarray, Me):
    """`sum_f R_f D_f X0`, `D_f` the face-`f` part of the CR derivative
    (same stencil weights as `ops.cr_derivative`), scaled by `M^{-1}`."""
    nE = len(m.E)
    _, area = ops.face_normals(m, X0)
    Fr = np.zeros((2 * nE, 3))
    for f, t in enumerate(m.F):
        L = [np.linalg.norm(X0[t[(i + 1) % 3]] - X0[t[i]]) for i in range(3)]
        for i in range(3):
            e = m.fe[f, i]
            tail, head, opp = t[i], t[(i + 1) % 3], t[(i + 2) % 3]
            si = 1.0 if tail < head else -1.0
            li = L[i]
            ejk, eki = L[(i + 1) % 3] ** 2, L[(i + 2) % 3] ** 2
            a = si * area[f] / (3 * li)
            k = -si / (12 * li)
            gp = a * (X0[head] - X0[tail])
            gq = (k * (li * li - (eki - ejk)) * X0[tail]
                  + k * (li * li + (eki - ejk)) * X0[head] + si * li / 6 * X0[opp])
            Fr[e] += R[f] @ gp
            Fr[nE + e] += R[f] @ gq
    return Fr / np.concatenate([Me, Me])[:, None]


def run_rotations(m: ops.Mesh, X_init: np.ndarray, X_rest: np.ndarray,
                  s_goal4: np.ndarray, cfg: RotationConfig = RotationConfig()):
    s3g = cl.sym_goal(s_goal4)
    D0, _ = ops.cr_derivative(m, X_rest)
    edges = _interior_edges(m)
    th0 = dihedrals(m, X_rest, edges)
    _, area = ops.face_normals(m, X_rest)
    wA = area / area.mean()
    if cfg.global_step == "cot":
        wfi, solve = _face_arap_system(m, X_rest, cfg.eps)
    else:
        D0, Me, M2, solve = cl._factor(m, X_rest, cfg.eps)
    X = X_init.copy()
    res = LoopResult(X=X)
    for _ in range(cfg.iterations):
        th = target_dihedrals(m, X, D0, s3g, edges)
        R = sync_rotations(m, X_rest, edges, th - th0,
                           _face_procrustes(m, X_rest, X), wA, cfg.mu)
        if cfg.global_step == "cot":
            rhs = _face_arap_rhs(m, X_rest, R, wfi) + cfg.eps * X
        else:
            rhs = D0.T @ (M2 @ _cr_rotated_frames(m, X_rest, R, Me)) + cfg.eps * X
        X_new = np.column_stack([solve(rhs[:, k]) for k in range(3)])
        step = float(np.max(np.linalg.norm(X_new - X, axis=1)))
        X = X_new
        res.curvature.append(cl.curvature_residual(m, D0, X, s3g))
        res.energy.append(float(np.max(np.abs(dihedrals(m, X, edges) - th))))
        res.step.append(step)
        if step < cfg.tol:
            break
    res.X = X
    return res


# ------------------------------------------------------------------- main


def strip_case(ns=17, nz=6, radius=1.5):
    import base_loop as bl
    P, F = ops.strip_param(ns, nz, 2.0, 1.0)
    Xf = ops.flat_from_param(P)
    Xc = ops.cylinder_from_param(P, radius)
    m = ops.build_mesh(F, len(Xf))
    s4 = bl.readout(m, ops.cr_derivative(m, Xc)[0], Xc)
    return m, Xf, Xc, s4


if __name__ == "__main__":
    import time
    import base_loop as bl

    m, Xf, Xc, s4 = strip_case()
    D0, _ = ops.cr_derivative(m, Xf)
    st = build_stencils(m, D0)
    s_ref = cl.sym_goal(bl.readout(m, D0, Xc))
    print("stencil read-out check", np.abs(stencil_s3(st, Xc[st.ids]) - s_ref).max())

    for name, fn in (("stencil", run_stencil), ("rotations", run_rotations)):
        for X0, label in ((Xc, "idle cyl"), (Xf, "flat->cyl")):
            t = time.time()
            r = fn(m, X0, Xf, s4)
            print(f"{name:9s} {label:9s} its {len(r.step):3d} "
                  f"err {bl.rigid_error(r.X, Xc):.3e} "
                  f"curv {r.curvature[0]:.4f} -> {r.curvature[-1]:.2e} "
                  f"step {r.step[-1]:.1e} t {time.time() - t:.1f}s")
