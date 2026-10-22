"""Mesh operators for the CR Poisson shape exploration (phase B, scratch).

Ports of the C++ pieces this exploration builds on:
- CR derivative of a vertex function (Stein et al. 2020, eq. 7), as in
  `hepworth::block::cr_coordinate_gradients`.
- Face shape operator from the CR derivative of the embedding (App. A.2),
  as in `hepworth::block::cr_face_operators`.
- Vector Dirichlet stiffness on CR edge DOFs, as in
  `kusama::laplace::vector_dirichlet::energy`.
- Cotan Laplacian and ARAP spokes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

K_CR_SHAPE_SIGN = -1.0


@dataclass
class Mesh:
    F: np.ndarray  # (m, 3) faces, outward winding
    E: np.ndarray  # (nE, 2) edges, E[e, 0] < E[e, 1]
    fe: np.ndarray  # (m, 3) edge of local edge i = (F[f, i], F[f, (i + 1) % 3])
    edge_faces: list  # incident faces per edge
    n: int


def build_mesh(F: np.ndarray, n: int) -> Mesh:
    F = np.asarray(F, dtype=int)
    index = {}
    E = []
    fe = np.zeros_like(F)
    edge_faces = []
    for f, t in enumerate(F):
        for i in range(3):
            a, b = int(t[i]), int(t[(i + 1) % 3])
            key = (min(a, b), max(a, b))
            if key not in index:
                index[key] = len(E)
                E.append(key)
                edge_faces.append([])
            e = index[key]
            fe[f, i] = e
            edge_faces[e].append(f)
    return Mesh(F=F, E=np.array(E, dtype=int), fe=fe, edge_faces=edge_faces, n=n)


# ---------------------------------------------------------------- meshes


def torus(R: float, r: float, nu: int, nv: int):
    """`asawa::make_torus` layout, outward winding. Returns X, F, phi."""
    X = np.zeros((nu * nv, 3))
    phi = np.zeros(nu * nv)
    for i in range(nu):
        u = 2 * np.pi * i / nu
        for j in range(nv):
            v = 2 * np.pi * j / nv
            rho = R + r * np.cos(v)
            X[i * nv + j] = (rho * np.cos(u), rho * np.sin(u), r * np.sin(v))
            phi[i * nv + j] = v

    def vid(i, j):
        return (i % nu) * nv + (j % nv)

    F = []
    for i in range(nu):
        for j in range(nv):
            F.append((vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)))
            F.append((vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)))
    return X, np.array(F), phi


def strip_param(ns: int, nz: int, width: float, height: float):
    """Grid in (s, z), s across `width`, z across `height`, centered.
    Winding gives +z normals in the plane and outward normals on a cylinder
    built with `cylinder_from_param`."""
    s = np.linspace(-width / 2, width / 2, ns)
    z = np.linspace(-height / 2, height / 2, nz)
    S, Z = np.meshgrid(s, z, indexing="ij")
    P = np.stack([S.ravel(), Z.ravel()], axis=1)

    def vid(i, j):
        return i * nz + j

    F = []
    for i in range(ns - 1):
        for j in range(nz - 1):
            F.append((vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)))
            F.append((vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)))
    return P, np.array(F)


def flat_from_param(P: np.ndarray) -> np.ndarray:
    """Flat strip in the plane y = 0, s along x, z along z. Normals face -y,
    which is the outward side of the matching cylinder at angle 0."""
    return np.stack([P[:, 0], np.zeros(len(P)), P[:, 1]], axis=1)


def cylinder_from_param(P: np.ndarray, radius: float) -> np.ndarray:
    """Isometric wrap of the strip around the z axis, arc length s, centered
    at angle -pi/2 so the flat strip (y = 0 plane, -y normals) is its tangent
    plane at the middle."""
    t = P[:, 0] / radius - np.pi / 2
    return np.stack(
        [radius * np.cos(t), radius * np.sin(t) + radius, P[:, 1]], axis=1
    )


# ---------------------------------------------------------- basic geometry


def face_normals(m: Mesh, X: np.ndarray):
    a, b, c = X[m.F[:, 0]], X[m.F[:, 1]], X[m.F[:, 2]]
    cr = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(cr, axis=1)
    return cr / (2 * area[:, None]), area


def vertex_areas(m: Mesh, X: np.ndarray) -> np.ndarray:
    _, area = face_normals(m, X)
    A = np.zeros(m.n)
    for k in range(3):
        np.add.at(A, m.F[:, k], area / 3)
    return A


def cotan_weights(m: Mesh, X: np.ndarray) -> np.ndarray:
    """Per edge `w_e = (cot a + cot b) / 2` over the opposite angles."""
    w = np.zeros(len(m.E))
    for f, t in enumerate(m.F):
        for i in range(3):
            a, b, o = t[i], t[(i + 1) % 3], t[(i + 2) % 3]
            u, v = X[a] - X[o], X[b] - X[o]
            cot = np.dot(u, v) / np.linalg.norm(np.cross(u, v))
            w[m.fe[f, i]] += 0.5 * cot
    return w


def cotan_laplacian(m: Mesh, X: np.ndarray) -> sp.csr_matrix:
    """PSD `L = sum_e w_e (e_i - e_j)(e_i - e_j)^T`."""
    w = cotan_weights(m, X)
    i, j = m.E[:, 0], m.E[:, 1]
    rows = np.concatenate([i, j, i, j])
    cols = np.concatenate([j, i, i, j])
    vals = np.concatenate([-w, -w, w, w])
    return sp.csr_matrix((vals, (rows, cols)), shape=(m.n, m.n))


# ------------------------------------------------------------ CR derivative


def cr_derivative(m: Mesh, X: np.ndarray):
    """`G = M^{-1} D u` for vertex functions `u`: rows `[0, nE)` are the
    parallel derivative, rows `[nE, 2 nE)` the perpendicular derivative.
    Returns the sparse `D0` (2 nE x n) and the edge mass `A_e / 3` (nE)."""
    nE = len(m.E)
    _, area = face_normals(m, X)
    rows, cols, vals = [], [], []
    Ae = np.zeros(nE)
    for f, t in enumerate(m.F):
        L = [np.linalg.norm(X[t[(i + 1) % 3]] - X[t[i]]) for i in range(3)]
        A = area[f]
        for i in range(3):
            e = m.fe[f, i]
            tail, head, opp = t[i], t[(i + 1) % 3], t[(i + 2) % 3]
            si = 1.0 if tail < head else -1.0
            li = L[i]
            ejk = L[(i + 1) % 3] ** 2
            eki = L[(i + 2) % 3] ** 2
            Ae[e] += A
            a = si * A / (3 * li)
            rows += [e, e]
            cols += [head, tail]
            vals += [a, -a]
            k = -si / (12 * li)
            rows += [nE + e, nE + e, nE + e]
            cols += [tail, head, opp]
            vals += [k * (li * li - (eki - ejk)), k * (li * li + (eki - ejk)),
                     si * li / 6]
    Me = Ae / 3
    D = sp.csr_matrix((vals, (rows, cols)), shape=(2 * nE, m.n))
    Minv = sp.diags(np.concatenate([1 / Me, 1 / Me]))
    return (Minv @ D).tocsr(), Me


def cr_dirichlet(m: Mesh, X: np.ndarray):
    """`D0^T M0 D0` on vertex functions (per coordinate)."""
    D0, Me = cr_derivative(m, X)
    M = sp.diags(np.concatenate([Me, Me]))
    return (D0.T @ M @ D0).tocsr()


# --------------------------------------------------------- shape operator


def face_shape_from_G(m: Mesh, X: np.ndarray, Gpar: np.ndarray,
                      Gperp: np.ndarray) -> np.ndarray:
    """Face operators `(m, 3, 3)` from per-edge derivative vectors, live
    frames and normals (App. A.2 dotted with N, times `K_CR_SHAPE_SIGN`)."""
    N, area = face_normals(m, X)
    S = np.zeros((len(m.F), 3, 3))
    for f, t in enumerate(m.F):
        for i in range(3):
            e = m.fe[f, i]
            tail, head = t[i], t[(i + 1) % 3]
            vpar = X[head] - X[tail]
            li = np.linalg.norm(vpar)
            vpar /= li
            si = 1.0 if tail < head else -1.0
            vperp = np.cross(N[f], vpar)
            vperp /= np.linalg.norm(vperp)
            scale = si * li / area[f]
            cp = scale * N[f] @ Gpar[e]
            cq = scale * N[f] @ Gperp[e]
            S[f] += cp * np.outer(vperp, vpar) + cq * np.outer(vperp, vperp)
    return K_CR_SHAPE_SIGN * S


def edge_G(m: Mesh, D0: sp.csr_matrix, X: np.ndarray):
    nE = len(m.E)
    G = D0 @ X
    return G[:nE], G[nE:]


def face_shape(m: Mesh, X: np.ndarray, D0: sp.csr_matrix | None = None):
    if D0 is None:
        D0, _ = cr_derivative(m, X)
    Gp, Gq = edge_G(m, D0, X)
    return face_shape_from_G(m, X, Gp, Gq)


def face_frames(m: Mesh, X: np.ndarray):
    """`T` along the first local edge, `B = N x T`."""
    N, _ = face_normals(m, X)
    T = X[m.F[:, 1]] - X[m.F[:, 0]]
    T /= np.linalg.norm(T, axis=1)[:, None]
    B = np.cross(N, T)
    return T, B, N


def face_coeffs(S: np.ndarray, T: np.ndarray, B: np.ndarray) -> np.ndarray:
    """`(m, 4)`: `[T S T, T S B, B S T, B S B]` in each face's own frame."""
    out = np.zeros((len(S), 4))
    out[:, 0] = np.einsum("fi,fij,fj->f", T, S, T)
    out[:, 1] = np.einsum("fi,fij,fj->f", T, S, B)
    out[:, 2] = np.einsum("fi,fij,fj->f", B, S, T)
    out[:, 3] = np.einsum("fi,fij,fj->f", B, S, B)
    return out


def shape_readout_rows(m: Mesh, X: np.ndarray):
    """Sparse `B` with `face_coeffs(face_shape_from_G(G)) = B @ g`, where
    `g` stacks `[Gpar (nE x 3), Gperp (nE x 3)]` row-major (6 per edge).
    Linear in `G` with the live normals and frames of `X` frozen."""
    N, area = face_normals(m, X)
    T, Bf, _ = face_frames(m, X)
    nE = len(m.E)
    rows, cols, vals = [], [], []
    for f, t in enumerate(m.F):
        for i in range(3):
            e = m.fe[f, i]
            tail, head = t[i], t[(i + 1) % 3]
            vpar = X[head] - X[tail]
            li = np.linalg.norm(vpar)
            vpar /= li
            si = 1.0 if tail < head else -1.0
            vperp = np.cross(N[f], vpar)
            vperp /= np.linalg.norm(vperp)
            scale = K_CR_SHAPE_SIGN * si * li / area[f]
            Mp = np.outer(vperp, vpar)
            Mq = np.outer(vperp, vperp)
            for which, M in ((0, Mp), (1, Mq)):
                co = np.array([T[f] @ M @ T[f], T[f] @ M @ Bf[f],
                               Bf[f] @ M @ T[f], Bf[f] @ M @ Bf[f]])
                for r in range(4):
                    for ax in range(3):
                        rows.append(4 * f + r)
                        cols.append(6 * e + 3 * which + ax)
                        vals.append(scale * co[r] * N[f, ax])
    return sp.csr_matrix((vals, (rows, cols)), shape=(4 * len(m.F), 6 * nE))


def pack_G(Gp: np.ndarray, Gq: np.ndarray) -> np.ndarray:
    return np.concatenate([Gp, Gq], axis=1).ravel()


def unpack_G(g: np.ndarray, nE: int):
    G = g.reshape(nE, 6)
    return G[:, :3], G[:, 3:]


# --------------------------------------------------------- vector Dirichlet


def vector_dirichlet(m: Mesh, X: np.ndarray) -> sp.csr_matrix:
    """Port of `kusama::laplace::vector_dirichlet::energy` (2 nE x 2 nE)."""
    nE = len(m.E)
    rows, cols, vals = [], [], []

    def sgn(a, b):
        return 1 if a <= b else -1

    for f, t in enumerate(m.F):
        v = [int(t[0]), int(t[1]), int(t[2])]
        lsq = [np.sum((X[v[(k + 1) % 3]] - X[v[k]]) ** 2) for k in range(3)]
        dA = np.linalg.norm(np.cross(X[v[1]] - X[v[0]], X[v[2]] - X[v[0]]))
        for k in range(3):
            eij, ejk, eki = lsq[k], lsq[(k + 1) % 3], lsq[(k + 2) % 3]
            lens = np.sqrt(eij * eki)
            o = sgn(v[k], v[(k + 1) % 3]) * sgn(v[(k + 2) % 3], v[k])
            Ei = m.fe[f, k]
            Ej = m.fe[f, (k + 2) % 3]
            diag = 2.0 / dA * eij
            num = eij - ejk + eki
            Dijki = o * num * num / (2.0 * eij * eki * dA) * lens
            Dperp = -o * num / (eij * eki) * lens
            for r, c, x in (
                (Ei, Ei, diag), (Ei + nE, Ei + nE, diag),
                (Ei, Ej, Dijki), (Ej, Ei, Dijki),
                (Ei + nE, Ej + nE, Dijki), (Ej + nE, Ei + nE, Dijki),
                (Ei, Ej + nE, Dperp), (Ej + nE, Ei, Dperp),
                (Ei + nE, Ej, -Dperp), (Ej, Ei + nE, -Dperp),
            ):
                rows.append(r)
                cols.append(c)
                vals.append(x)
    return sp.csr_matrix((vals, (rows, cols)), shape=(2 * nE, 2 * nE))


# ------------------------------------------------------------------- ARAP


def arap_rhs_and_rotations(m: Mesh, X0: np.ndarray, X: np.ndarray,
                           w: np.ndarray):
    """Per-vertex Procrustes rotations on cotan-weighted spokes and the ARAP
    right-hand side `b` with `L_cot X = b` for those rotations."""
    n = m.n
    C = np.zeros((n, 3, 3))
    for e, (i, j) in enumerate(m.E):
        d0 = X0[j] - X0[i]
        d = X[j] - X[i]
        C[i] += w[e] * np.outer(d, d0)
        C[j] += w[e] * np.outer(d, d0)
    R = np.zeros((n, 3, 3))
    for i in range(n):
        U, _, Vt = np.linalg.svd(C[i])
        Ri = U @ Vt
        if np.linalg.det(Ri) < 0:
            U[:, -1] *= -1
            Ri = U @ Vt
        R[i] = Ri
    b = np.zeros((n, 3))
    for e, (i, j) in enumerate(m.E):
        d0 = X0[j] - X0[i]
        r = 0.5 * w[e] * (R[i] + R[j]) @ d0
        b[j] += r
        b[i] -= r
    return R, b
