"""Dirichlet energies for a curvature-tensor field, on edges or on faces.

A symmetric 2-tensor in an orthonormal frame is `(k11, k12, k22)`; its trace
`tau = k11 + k22` is spin 0 and its trace-free part
`z = ((k11 - k22) / 2, k12)` is spin 2. The energy measured is

    E = int |grad z|^2 + 1/2 |grad tau|^2.

Edge element (CR, the vector Dirichlet element's site). Each edge stores its
tensor in its own frame `(t_e, b_e)`, `t_e` the canonical edge direction
(lower vertex id to higher) and `b_e = N x t_e` (unfolded per face, so both
faces agree). Per face, rotate each edge's tensor into the face frame and
interpolate with the CR basis `1 - 2 lambda_opposite`:

    E_f = sum_{a, b in f} K_ab q(M_a h_a, M_b h_b),

`K` the CR scalar stiffness (4 x the P1 stiffness of the opposite vertices),
`M_a` the frame map of edge `a` into the face frame, `q` the inner product
above. At spin 1 (`M_a` a rotation, `q` the dot product) this must reproduce
`ops.vector_dirichlet` (Stein et al. 2020).

Face element (dual graph). Each face stores its tensor in its frame
`(T_f, B_f)`. Per interior edge `e = (f, g)`, read both faces' tensors in the
unfolded edge frame and pay `w_e q(d, d)` for the difference `d`, with
`w_e = |e| / |e*|`: circumcentric (`2 / (cot alpha + cot beta)`) or
barycentric (unfolded centroid distance).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

import cr_ops as ops

# q(c, c') for c = (k11, k12, k22): z . z' + 1/2 tau tau'.
Q = np.array([[0.25 + 0.5, 0.0, -0.25 + 0.5],
              [0.0, 1.0, 0.0],
              [-0.25 + 0.5, 0.0, 0.25 + 0.5]])


def rot2(phi: float) -> np.ndarray:
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c, -s], [s, c]])


def tensor_map(phi: float) -> np.ndarray:
    """(k11, k12, k22) in a frame rotated by `phi` -> the same tensor in the
    unrotated frame (`S = R S' R^T`)."""
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c * c, -2 * c * s, s * s],
                     [c * s, c * c - s * s, -c * s],
                     [s * s, 2 * c * s, c * c]])


def edge_directions(m: ops.Mesh, X: np.ndarray) -> np.ndarray:
    t = X[m.E[:, 1]] - X[m.E[:, 0]]
    return t / np.linalg.norm(t, axis=1)[:, None]


def edge_normals(m: ops.Mesh, X: np.ndarray) -> np.ndarray:
    N, _ = ops.face_normals(m, X)
    n = np.zeros((len(m.E), 3))
    for e, fs in enumerate(m.edge_faces):
        for f in fs:
            n[e] += N[f]
    return n / np.linalg.norm(n, axis=1)[:, None]


def cr_face_stiffness(P: np.ndarray) -> np.ndarray:
    """CR scalar stiffness of one triangle (3 x 3), indexed by local edge
    `a = (P[a], P[a+1])` (opposite vertex `a + 2`)."""
    cot = np.zeros(3)
    for k in range(3):
        u, v = P[(k + 1) % 3] - P[k], P[(k + 2) % 3] - P[k]
        cot[k] = np.dot(u, v) / np.linalg.norm(np.cross(u, v))
    K1 = np.zeros((3, 3))  # P1 stiffness, vertex indexed
    for k in range(3):
        i, j = (k + 1) % 3, (k + 2) % 3
        K1[i, j] = K1[j, i] = -0.5 * cot[k]
    for i in range(3):
        K1[i, i] = -K1[i].sum()
    opp = [(a + 2) % 3 for a in range(3)]
    return 4.0 * K1[np.ix_(opp, opp)]


def _face_edge_angles(m: ops.Mesh, X: np.ndarray, f: int, T, B, t_e):
    """Angle of each local edge's canonical direction in the face frame."""
    return [np.arctan2(t_e[m.fe[f, a]] @ B[f], t_e[m.fe[f, a]] @ T[f]) for a in range(3)]


def edge_dirichlet(m: ops.Mesh, X: np.ndarray, spin: int) -> sp.csr_matrix:
    """Spin 1: (2 nE)^2 in the `[par (nE), perp (nE)]` layout of
    `ops.vector_dirichlet`. Spin 2: (3 nE)^2, `(k11, k12, k22)` per edge,
    edge-major. Energy `c^T L c` (no 1/2)."""
    T, B, _ = ops.face_frames(m, X)
    t_e = edge_directions(m, X)
    nE = len(m.E)
    rows, cols, vals = [], [], []
    for f, tri in enumerate(m.F):
        K = cr_face_stiffness(X[tri])
        phi = _face_edge_angles(m, X, f, T, B, t_e)
        if spin == 1:
            Ms = [rot2(p) for p in phi]
            q = np.eye(2)
        else:
            Ms = [tensor_map(p) for p in phi]
            q = Q
        for a in range(3):
            for b in range(3):
                blk = K[a, b] * Ms[a].T @ q @ Ms[b]
                ea, eb = m.fe[f, a], m.fe[f, b]
                d = blk.shape[0]
                for r in range(d):
                    for c in range(d):
                        if spin == 1:
                            rows.append(ea + r * nE)
                            cols.append(eb + c * nE)
                        else:
                            rows.append(3 * ea + r)
                            cols.append(3 * eb + c)
                        vals.append(blk[r, c])
    n = (2 if spin == 1 else 3) * nE
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def face_dirichlet(m: ops.Mesh, X: np.ndarray, weights: str = "circum") -> sp.csr_matrix:
    """(3 nF)^2 on `(k11, k12, k22)` per face in the face frame. Energy
    `c^T L c`. `weights`: "circum" (2 / (cot alpha + cot beta)) or "bary"."""
    T, B, N = ops.face_frames(m, X)
    t_e = edge_directions(m, X)
    nF = len(m.F)
    rows, cols, vals = [], [], []
    cen = X[m.F].mean(1)
    for e, fs in enumerate(m.edge_faces):
        if len(fs) != 2:
            continue
        f, g = fs
        a, b = X[m.E[e, 0]], X[m.E[e, 1]]
        l = np.linalg.norm(b - a)
        if weights == "circum":
            cots = []
            for h in (f, g):
                o = [v for v in m.F[h] if v not in m.E[e]][0]
                u, v = a - X[o], b - X[o]
                cots.append(np.dot(u, v) / np.linalg.norm(np.cross(u, v)))
            w = 2.0 / (cots[0] + cots[1])
        else:
            # unfolded centroid distance: in-edge offset plus the two heights
            mid = 0.5 * (a + b)
            d = 0.0
            along = 0.0
            for h, sgn in ((f, 1.0), (g, -1.0)):
                r = cen[h] - mid
                along += sgn * (r @ t_e[e])
                d += np.linalg.norm(r - (r @ t_e[e]) * t_e[e])
            w = l / np.hypot(d, along)
        Ms = []
        for h in (f, g):
            phi = np.arctan2(t_e[e] @ B[h], t_e[e] @ T[h])
            # face-frame comps -> edge-frame comps: inverse of tensor_map(phi)
            Ms.append(tensor_map(-phi))
        Mf, Mg = Ms
        blocks = {(f, f): Mf.T @ Q @ Mf, (g, g): Mg.T @ Q @ Mg,
                  (f, g): -Mf.T @ Q @ Mg, (g, f): -Mg.T @ Q @ Mf}
        for (i, j), blk in blocks.items():
            for r in range(3):
                for c in range(3):
                    rows.append(3 * i + r)
                    cols.append(3 * j + c)
                    vals.append(w * blk[r, c])
    return sp.csr_matrix((vals, (rows, cols)), shape=(3 * nF, 3 * nF))


# ------------------------------------------------- edge readings (coupling)


def face_edge_readings(st, P: np.ndarray) -> np.ndarray:
    """(m, 3, 3): each face's live read-out read in each of its edges' frames
    `(t_e, N_f x t_e)` (unfolded; the edge orientation does not matter: both
    `t` and `b` flip together)."""
    import face_loops as fl
    s = fl.stencil_s3_live(st, P)  # face frame (T along P1 - P0, B = N x T)
    a, b, c = P[:, 0], P[:, 1], P[:, 2]
    N = np.cross(b - a, c - a)
    N /= np.linalg.norm(N, axis=1)[:, None]
    T = (b - a) / np.linalg.norm(b - a, axis=1)[:, None]
    B = np.cross(N, T)
    out = np.zeros((len(P), 3, 3))
    for i in range(3):
        v = P[:, (i + 1) % 3] - P[:, i]
        phi = np.arctan2(np.sum(v * B, 1), np.sum(v * T, 1))
        cs, sn = np.cos(phi), np.sin(phi)
        k11, k12, k22 = s[:, 0], s[:, 1], s[:, 2]
        out[:, i, 0] = cs * cs * k11 + 2 * cs * sn * k12 + sn * sn * k22
        out[:, i, 1] = -cs * sn * k11 + (cs * cs - sn * sn) * k12 + cs * sn * k22
        out[:, i, 2] = sn * sn * k11 - 2 * cs * sn * k12 + cs * cs * k22
    return out


def edge_average_readings(m, st, X):
    """(nE, 3): the mean of an edge's faces' readings in its frame."""
    R = face_edge_readings(st, X[st.ids])
    acc = np.zeros((len(m.E), 3))
    cnt = np.zeros(len(m.E))
    for i in range(3):
        np.add.at(acc, m.fe[:, i], R[:, i])
        np.add.at(cnt, m.fe[:, i], 1)
    return acc / cnt[:, None]


# ------------------------------------------------------------- sampling


def sample_edges(m, X, field):
    """Edge-frame comps of a tensor field `field(p, N) -> 3x3 (3D, tangent)`
    at edge midpoints, frame `(t_e, N_e x t_e)`."""
    t = edge_directions(m, X)
    n = edge_normals(m, X)
    b = np.cross(n, t)
    mid = 0.5 * (X[m.E[:, 0]] + X[m.E[:, 1]])
    out = np.zeros((len(m.E), 3))
    for e in range(len(m.E)):
        S = field(mid[e], n[e])
        out[e] = (t[e] @ S @ t[e], t[e] @ S @ b[e], b[e] @ S @ b[e])
    return out


def sample_faces(m, X, field):
    T, B, N = ops.face_frames(m, X)
    cen = X[m.F].mean(1)
    out = np.zeros((len(m.F), 3))
    for f in range(len(m.F)):
        S = field(cen[f], N[f])
        out[f] = (T[f] @ S @ T[f], T[f] @ S @ B[f], B[f] @ S @ B[f])
    return out


def scalar_cr_laplacian(m, X):
    nE = len(m.E)
    rows, cols, vals = [], [], []
    for f, tri in enumerate(m.F):
        K = cr_face_stiffness(X[tri])
        for a in range(3):
            for b in range(3):
                rows.append(m.fe[f, a])
                cols.append(m.fe[f, b])
                vals.append(K[a, b])
    return sp.csr_matrix((vals, (rows, cols)), shape=(nE, nE))


# ---------------------------------------------------------------- checks


def _planar(nx=7, ny=6, jitter=0.25, seed=0):
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.linspace(0, 1, nx), np.linspace(0, 1, ny), indexing="ij")
    P = np.stack([xs.ravel(), ys.ravel(), np.zeros(xs.size)], 1)
    h = 1.0 / (nx - 1)
    inner = (P[:, 0] > 0) & (P[:, 0] < 1) & (P[:, 1] > 0) & (P[:, 1] < 1)
    P[inner, :2] += jitter * h * rng.uniform(-1, 1, (inner.sum(), 2))
    F = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            a, b, c, d = i * ny + j, (i + 1) * ny + j, (i + 1) * ny + j + 1, i * ny + j + 1
            F += [(a, b, c), (a, c, d)]
    return P, np.array(F)


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    # 1. spin-1 construction reproduces Stein's vector Dirichlet matrix.
    for name, (X, F) in (("planar jittered", _planar()),
                         ("torus", ops.torus(2.0, 0.7, 10, 8)[:2])):
        m = ops.build_mesh(F, len(X))
        L1 = edge_dirichlet(m, X, spin=1).toarray()
        Ls = ops.vector_dirichlet(m, X).toarray()
        c = (L1 * Ls).sum() / (L1 * L1).sum()
        print(f"{name:16s} spin 1 vs vector_dirichlet: scale {c:.6f}, "
              f"rel diff {np.abs(c * L1 - Ls).max() / np.abs(Ls).max():.2e}")
    # 2. constant tensor on a planar mesh has zero energy; 3. rotation
    # invariance; 4. trace part equals the scalar Laplacian.
    X, F = _planar()
    m = ops.build_mesh(F, len(X))
    S0 = np.array([[1.0, 0.3, 0.0], [0.3, -0.5, 0.0], [0.0, 0.0, 0.0]])
    const = lambda p, n: S0  # noqa: E731
    he = sample_edges(m, X, const)
    hf = sample_faces(m, X, const)
    Le = edge_dirichlet(m, X, 2)
    Lfc = face_dirichlet(m, X, "circum")
    Lfb = face_dirichlet(m, X, "bary")
    ce, cf = he.ravel(), hf.ravel()
    print(f"constant tensor: edge {ce @ Le @ ce:.2e}  face circum {cf @ Lfc @ cf:.2e}  "
          f"face bary {cf @ Lfb @ cf:.2e}")
    def field(p, n):  # noqa: E306
        a = 1 + 0.5 * np.sin(2 * np.pi * p[1])
        bt = np.pi * p[0]
        z = a * np.array([np.cos(2 * bt), np.sin(2 * bt)])
        tau = 1 + 0.3 * p[0]
        return np.array([[tau / 2 + z[0], z[1], 0], [z[1], tau / 2 - z[0], 0], [0, 0, 0]])
    Rz = lambda th: np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])  # noqa: E731
    th = 0.7
    Xr = X @ Rz(th).T
    fr = lambda p, n: Rz(th) @ field(Rz(-th) @ p, n) @ Rz(th).T  # noqa: E731
    for nm, Lfun, samp in (("edge", lambda Y: edge_dirichlet(m, Y, 2), sample_edges),
                           ("face circum", lambda Y: face_dirichlet(m, Y, "circum"), sample_faces),
                           ("face bary", lambda Y: face_dirichlet(m, Y, "bary"), sample_faces)):
        c0 = samp(m, X, field).ravel()
        c1 = samp(m, Xr, fr).ravel()
        E0, E1 = c0 @ Lfun(X) @ c0, c1 @ Lfun(Xr) @ c1
        print(f"rotation invariance {nm:11s}: {E0:.6f} vs {E1:.6f}")
    trace_only = lambda p, n: np.diag([np.sin(3 * p[0]) + p[1], np.sin(3 * p[0]) + p[1], 0.0])  # noqa: E731
    he = sample_edges(m, X, trace_only)
    tau = he[:, 0] + he[:, 2]
    Ls = scalar_cr_laplacian(m, X)
    ce = he.ravel()
    print(f"trace part: edge {ce @ Le @ ce:.6f} vs 1/2 tau^T L_CR tau {0.5 * tau @ Ls @ tau:.6f}")
